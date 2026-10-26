"""Explicit GitHub transport for synthetic evidence, never for production approval.

Uses GITHUB_TOKEN or, only with --git-credential, the repository's existing Git
credential helper. Credentials are kept in memory and never printed or stored.
No endpoint here merges a PR, changes visibility, or enables a real runner.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import urllib.error
import urllib.request
import uuid

from gate import digest, evaluate, merge_decision
from run import ROOT, read, verify, write

REPOSITORY = "KoinoSora/BSpline2Analytic"
ALLOWED_FILES = {"context.json", "policy.json", "registry.json", "evidence.json", "result.json",
                 "lifecycle.json", "production-refusal.json", "measurements.csv", "report.md", "manifest.json"}


class GitHubError(RuntimeError):
    def __init__(self, status, method, path):
        super().__init__(f"GitHub HTTP {status}: {method} {path}")
        self.status = status


class GitHub:
    def __init__(self, use_git_credential=False):
        self.token = os.environ.get("GITHUB_TOKEN")
        if not self.token and use_git_credential:
            origin = subprocess.run(["git", "remote", "get-url", "origin"], check=True,
                                    capture_output=True, text=True).stdout.strip()
            if origin not in (f"https://github.com/{REPOSITORY}.git", f"https://github.com/{REPOSITORY}"):
                raise ValueError("Credential use is restricted to the configured project origin")
            env = dict(os.environ, GIT_TERMINAL_PROMPT="0", GCM_INTERACTIVE="Never")
            filled = subprocess.run(["git", "credential", "fill"],
                                    input="protocol=https\nhost=github.com\n\n",
                                    capture_output=True, text=True, env=env, timeout=30)
            if filled.returncode:
                raise ValueError("Git credential helper did not provide credentials")
            fields = dict(line.split("=", 1) for line in filled.stdout.splitlines() if "=" in line)
            self.token = fields.get("password")
        if not self.token:
            raise ValueError("Provide GITHUB_TOKEN or explicitly use --git-credential")

    def request(self, method, path, body=None):
        if not path.startswith(f"/repos/{REPOSITORY}/") and path != f"/repos/{REPOSITORY}":
            raise ValueError("API path outside target repository")
        raw = None if body is None else json.dumps(body, allow_nan=False).encode()
        request = urllib.request.Request("https://api.github.com" + path, data=raw, method=method,
                    headers={"Authorization": "Bearer " + self.token, "Accept": "application/vnd.github+json",
                             "X-GitHub-Api-Version": "2022-11-28", "Content-Type": "application/json",
                             "User-Agent": "BSpline2Analytic-closed-loop-protocol"})
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                data = response.read()
                return json.loads(data) if data else None
        except urllib.error.HTTPError as error:
            # Do not echo headers, request data or credential helper output.
            raise GitHubError(error.code, method, path) from None

    def repo(self, method, path="", body=None):
        return self.request(method, f"/repos/{REPOSITORY}" + path, body)


def live_context(api, pr_number):
    pr = api.repo("GET", f"/pulls/{pr_number}")
    if pr["state"] != "open" or pr["base"]["ref"] != "main":
        raise ValueError("Expected an open PR targeting main")
    base = api.repo("GET", "/git/ref/heads/main")["object"]["sha"]
    if pr["base"]["sha"] != base:
        raise ValueError("GitHub has not refreshed the PR base; retry after synchronization")
    merge_sha = pr.get("merge_commit_sha")
    if not merge_sha or pr.get("mergeable") is False:
        raise ValueError("PR merge tree is unavailable or conflicted")
    merge = api.repo("GET", f"/git/commits/{merge_sha}")
    parents = {p["sha"] for p in merge["parents"]}
    if not {base, pr["head"]["sha"]}.issubset(parents):
        raise ValueError("Synthetic merge commit does not bind current head and base")
    return dict(repository=REPOSITORY, base_branch="main", head_sha=pr["head"]["sha"],
                base_sha=base, merge_tree_sha=merge["tree"]["sha"], author=pr["user"]["login"],
                trusted_collector=False)


def publish_dry_run(api, bundle, pr_number):
    manifest = verify(bundle)
    if set(manifest["files"]) | {"manifest.json"} != ALLOWED_FILES or manifest.get("synthetic") is not True:
        raise ValueError("Only the fixed synthetic protocol bundle may be uploaded")
    policy = read(ROOT / "policy.json")
    if policy["real_tests_enabled"]:
        raise ValueError("This bridge is limited to the inactive experimental rollout")
    context = live_context(api, pr_number)
    evidence = read(bundle / "evidence.json")
    registry = read(bundle / "registry.json")
    result = evaluate(policy, registry, evidence, context)
    if result["state"] != "SIMULATION_PASS" or result["production_eligible"]:
        raise ValueError("Live PR identity or simulation validation failed")
    if result != read(bundle / "result.json"):
        raise ValueError("Stored result disagrees with recomputation against the live PR")
    for path in evidence["documents"].values():
        if not re.fullmatch(r"docs/(algorithm|requirements)/[A-Za-z0-9_/-]+\.md", path):
            raise ValueError("Unexpected module/impact document path")
        document = api.repo("GET", f"/contents/{path}?ref={context['head_sha']}")
        if document.get("type") != "file" or document.get("size", 0) == 0:
            raise ValueError("Module/impact document is absent from the PR commit")
    branch = policy["evidence_branch"]
    if branch != "validation-evidence":
        raise ValueError("Unexpected evidence branch")
    previous = None
    try:
        previous = api.repo("GET", f"/git/ref/heads/{branch}")["object"]["sha"]
    except GitHubError as error:
        if error.status != 404:
            raise
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex[:12]
    prefix = f"docs/testing/closed-loop/pr-{pr_number}/{run_id}"
    entries = [dict(path=f"{prefix}/{name}", mode="100644", type="blob",
                    content=(bundle / name).read_bytes().decode("utf-8")) for name in sorted(ALLOWED_FILES)]
    tree_input = {"tree": entries}
    if previous:
        tree_input["base_tree"] = api.repo("GET", f"/git/commits/{previous}")["tree"]["sha"]
    tree = api.repo("POST", "/git/trees", tree_input)
    commit = api.repo("POST", "/git/commits", dict(message=f"Record synthetic closed-loop rehearsal for PR #{pr_number}",
                       tree=tree["sha"], parents=[previous] if previous else []))
    # No retry overwrites another publisher: recheck, then use a non-force ref update.
    if previous:
        if api.repo("GET", f"/git/ref/heads/{branch}")["object"]["sha"] != previous:
            raise ValueError("Evidence branch advanced; rerun publication without overwriting it")
        api.repo("PATCH", f"/git/refs/heads/{branch}", {"sha": commit["sha"], "force": False})
    else:
        api.repo("POST", "/git/refs", {"ref": "refs/heads/" + branch, "sha": commit["sha"]})
    remote = api.repo("GET", f"/git/trees/{tree['sha']}?recursive=1")
    if remote.get("truncated"):
        raise ValueError("Cannot verify a truncated remote evidence tree")
    remote_files = {x["path"]: x for x in remote["tree"]}
    # Verify every uploaded Git blob against local bytes, rather than trust a successful POST.
    import hashlib
    for name in ALLOWED_FILES:
        raw = (bundle / name).read_bytes()
        expected_blob = hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()
        if remote_files[f"{prefix}/{name}"]["sha"] != expected_blob:
            raise ValueError("Remote evidence blob mismatch")
    if api.repo("GET", f"/git/ref/heads/{branch}")["object"]["sha"] != commit["sha"]:
        raise ValueError("Evidence ref changed before publication finished")
    current = live_context(api, pr_number)
    if any(current[k] != context[k] for k in ("head_sha", "base_sha", "merge_tree_sha")):
        raise ValueError("PR changed during publication; old evidence remains archived but cannot pass")
    url = f"https://github.com/{REPOSITORY}/blob/{commit['sha']}/{prefix}/report.md"
    status = api.repo("POST", f"/statuses/{context['head_sha']}",
                      dict(state="success", context=policy["simulation_status"], target_url=url,
                           description="Synthetic protocol verified; production gate NOT ACTIVATED"))
    # Deliberately never write closed-loop/gate. A required production gate stays unsatisfied.
    return dict(verified=True, synthetic=True, production_eligible=False, pr=pr_number,
                evidence_sha256=digest(evidence), evidence_commit_sha=commit["sha"],
                report_url=url, status_context=status["context"], status_state=status["state"], **context)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--git-credential", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("inspect")
    cmd = sub.add_parser("context")
    cmd.add_argument("--pr", type=int, required=True)
    cmd.add_argument("--output", type=Path, required=True)
    cmd = sub.add_parser("publish-dry-run")
    cmd.add_argument("--pr", type=int, required=True)
    cmd.add_argument("--bundle", type=Path, required=True)
    cmd.add_argument("--receipt", type=Path, required=True)
    cmd.add_argument("--allow-network-writes", action="store_true", required=True)
    cmd = sub.add_parser("review-status")
    cmd.add_argument("--pr", type=int, required=True)
    cmd.add_argument("--bundle", type=Path, required=True)
    cmd.add_argument("--receipt", type=Path, required=True)
    cmd = sub.add_parser("create-pr")
    cmd.add_argument("--head", required=True)
    cmd.add_argument("--body-file", type=Path, required=True)
    cmd.add_argument("--allow-network-writes", action="store_true", required=True)
    args = parser.parse_args()
    try:
        api = GitHub(args.git_credential)
        if args.command == "inspect":
            repo = api.repo("GET")
            summary = {k: repo.get(k) for k in ("full_name", "visibility", "permissions", "default_branch")}
            try:
                summary["branch_protection"] = api.repo("GET", "/branches/main/protection")
            except GitHubError as error:
                summary["branch_protection_http"] = error.status
            print(json.dumps(summary, ensure_ascii=False, indent=2))
        elif args.command == "context":
            context = live_context(api, args.pr)
            write(args.output, context)
            print(json.dumps(context, ensure_ascii=False))
        elif args.command == "create-pr":
            if not re.fullmatch(r"[A-Za-z0-9_./-]+", args.head) or args.head == "main":
                raise ValueError("Invalid PR source branch")
            result = api.repo("POST", "/pulls", dict(base="main", head=args.head, draft=True,
                              title="Add dormant closed-loop validation agent and protocol rehearsal",
                              body=args.body_file.read_text(encoding="utf-8")))
            print(json.dumps({k: result[k] for k in ("number", "html_url", "state", "draft")}))
        elif args.command == "review-status":
            verify(args.bundle)
            context = live_context(api, args.pr)
            policy = read(ROOT / "policy.json")
            receipt = read(args.receipt)
            reviews = []
            page = 1
            while True:
                batch = api.repo("GET", f"/pulls/{args.pr}/reviews?per_page=100&page={page}")
                reviews.extend(dict(sequence=r["id"], login=r["user"]["login"],
                                    state=r["state"], head_sha=r["commit_id"]) for r in batch)
                if len(batch) < 100:
                    break
                page += 1
            result = read(args.bundle / "result.json")
            state = merge_decision(result, context, receipt, reviews, policy, simulation=True)
            print(json.dumps(dict(state=state, production_eligible=False,
                                  actual_reviews_read=len(reviews), author=context["author"])))
        else:
            receipt = publish_dry_run(api, args.bundle.resolve(), args.pr)
            write(args.receipt, receipt)
            print(json.dumps(receipt, ensure_ascii=False, indent=2))
        return 0
    except (ValueError, KeyError, OSError, GitHubError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
