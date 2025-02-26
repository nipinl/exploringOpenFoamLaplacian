#include "fvCFD.H" 
//************************************************************// 
int main(int argc, char *argv[]) 
{ 

    #include "setRootCase.H" 
    #include "createTime.H" 
    #include "createMesh.H" 
    #include "createFields.H" 
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
    label Nc = mesh.nCells(); 
    Info<<"Number of Cells =  "<<Nc<<endl; 
    fvScalarMatrix TEqn 
    (fvm::laplacian(DT, T) + su + fvm::Sp(sp, T)); 
    
        Info<<"Diagonal elements:"<<endl<<TEqn.diag()<<endl; 
        Info<<"Upper elements:"<<endl<<TEqn.upper()<<endl; 
        Info<<"Lower elements:"<<endl<<TEqn.lower()<<endl; 
        Info<<"boundaryCoeffs:"<<endl<<TEqn.boundaryCoeffs()<<endl; 
        Info<<"internalCoeffs:"<<endl<<TEqn.internalCoeffs()<<endl; 
    
    Info<<"Size of boundaryField is "<<T.boundaryField().size()<<endl; 
    #include "defineObjectsForWriting.H"
    // Assigning contribution from BC 
    
    forAll(T.boundaryField(), patchI) 
    { 
     const fvPatch &vfp = T.boundaryField()[patchI].patch(); 
     forAll(vfp, faceI) 
     { 
         label cellI = vfp.faceCells()[faceI]; 
	 //Info<<" CellI =  "<<cellI<<endl;


         wdiag[cellI] += TEqn.internalCoeffs()[patchI][faceI]; 
         wsource[cellI] += TEqn.boundaryCoeffs()[patchI][faceI]; 
     } 
    } 
    while (simple.loop(runTime))
    {
    }
    // write to file 
     wdiag.write(); 
     wsource.write(); 
     wlower.write(); 
     wupper.write(); 
    
     solve(TEqn); 
    
     runTime++; 
     runTime.write(); 
    
     Info<< nl; 
     Info<< "hai from myThermalConductionSolver\n" << endl; 
     return 0; 
}
