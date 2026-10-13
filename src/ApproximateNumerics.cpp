#include "ApproximateInternal.h"
#include <numeric>
namespace bspline2analytic::detail
{
bool SolveLeastSquares(std::vector<std::vector<double>> a, std::vector<double> b,
                       std::vector<double> &x)
{
    if (a.empty())
        return false;
    size_t m = a.size(), n = a[0].size();
    if (m < n || !n || b.size() != m)
        return false;
    std::vector<size_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    double first = 0;
    for (size_t k = 0; k < n; ++k)
    {
        size_t pivot = k;
        double best = -1;
        for (size_t j = k; j < n; ++j)
        {
            double norm = 0;
            for (size_t i = k; i < m; ++i)
                norm = std::hypot(norm, a[i][j]);
            if (norm > best)
            {
                best = norm;
                pivot = j;
            }
        }
        if (k == 0)
            first = best;
        if (!(best > std::max(1e-30, first * 1e-12)))
            return false;
        if (pivot != k)
        {
            for (auto &r : a)
                std::swap(r[k], r[pivot]);
            std::swap(perm[k], perm[pivot]);
        }
        double alpha = -std::copysign(best, a[k][k]);
        std::vector<double> v(m - k);
        for (size_t i = k; i < m; ++i)
            v[i - k] = a[i][k];
        v[0] -= alpha;
        double norm = 0;
        for (double z : v)
            norm = std::hypot(norm, z);
        for (auto &z : v)
            z /= norm;
        for (size_t j = k; j < n; ++j)
        {
            double d = 0;
            for (size_t i = k; i < m; ++i)
                d += v[i - k] * a[i][j];
            for (size_t i = k; i < m; ++i)
                a[i][j] -= 2 * v[i - k] * d;
        }
        double d = 0;
        for (size_t i = k; i < m; ++i)
            d += v[i - k] * b[i];
        for (size_t i = k; i < m; ++i)
            b[i] -= 2 * v[i - k] * d;
        a[k][k] = alpha;
    }
    std::vector<double> y(n);
    for (size_t j = n; j-- > 0;)
    {
        double r = b[j];
        for (size_t k = j + 1; k < n; ++k)
            r -= a[j][k] * y[k];
        y[j] = r / a[j][j];
    }
    x.assign(n, 0);
    for (size_t j = 0; j < n; ++j)
        x[perm[j]] = y[j];
    return std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); });
}
} // namespace bspline2analytic::detail
