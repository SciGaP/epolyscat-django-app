#include "basisMatMO.h"
#include "abort.h"
#include "basisMO.h"
#ifdef _USE_HACC_
#include "mo.h"
#endif
#include "printOutput.h"

BasisMatMO::BasisMatMO(std::string Op,const BasisMO* IBas, const BasisMO* JBas)
{
    if(IBas!=JBas)DEVABORT("for now, we cannot handle different BasisMO");

    // translate operator names to Vinay's
    std::string op;
    if(Op=="<<Laplacian>>")op="Laplacian";
    else if(Op=="<<Coulomb>>")op="Coulomb";
    else if(Op=="<<EEInteraction>>")op=Op;
    else DEVABORT("unknown Op for molecular orbitals (BasisMO): "+Op);

    if(op=="dummy"){
        _mat=Eigen::MatrixXcd::Zero(IBas->size(),JBas->size());
        PrintOutput::warning("dummy=0 for <"+IBas->str()+"|"+Op+"|"+JBas->str()+">");
    }
#ifdef _USE_HACC_
    else if(Op=="<<EEInteraction>>")
        _mat=IBas->vinayMO()->matrix_2e();
    else
        // retrieve from Vinay's
        const_cast<mo*>(IBas->vinayMO())->matrix_1e(_mat,op,nullptr);
#endif
}
