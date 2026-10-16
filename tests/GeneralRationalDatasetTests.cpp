#include "GeneralRationalDataset.h"
#include <Foundation/init.h>
#include <iostream>
using namespace bspline2analytic::tests;
int main(int argc,char **argv)
{
    sggk::init();try
    {
        if(argc==8&&std::string(argv[1])=="--generate")
            GenerateRationalDataset(argv[2],argv[3],argv[4],argv[5],unsigned(std::stoul(argv[6])),unsigned(std::stoul(argv[7])));
        else if(argc==6&&std::string(argv[1])=="--audit")AuditRationalDataset(argv[2],argv[3],unsigned(std::stoul(argv[4])),unsigned(std::stoul(argv[5])));
        else throw std::runtime_error("Usage: --generate legacy-main legacy-extreme output suite first last | --audit data suite first last");
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
