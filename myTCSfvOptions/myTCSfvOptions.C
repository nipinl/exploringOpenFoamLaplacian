#include "fvCFD.H"
#include "fvOptions.H"
//************************************************************// 
int main(int argc, char *argv[]) 
{ 

    #include "setRootCaseLists.H" 
    #include "createTime.H" 
    #include "createMesh.H" 
    #include "createFields.H" 
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 
    label Nc = mesh.nCells(); 
    Info<<"Number of Cells =  "<<Nc<<endl; 
    fvScalarMatrix TEqn
   		     (
		      //fvm::laplacian(DT, T) + su + fvm::Sp(sp, T) 
		      //above line replaced by the one below,added fvOptions in rhs, followed by constrain & correct are the  change added for fvOptions. Then added createFvOptions.H in the createFields.H
		      fvm::laplacian(DT, T) 
		      == 
		      fvOptions(T)
		     );
   		fvOptions.constrain(TEqn);
	        TEqn.solve();
		fvOptions.correct(T);	
    
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
    // write to file 
     wdiag.write(); 
     wsource.write(); 
     wlower.write(); 
     wupper.write(); 
    
     solve(TEqn); 
    
     runTime++; 
     runTime.write(); 
    
     Info<< nl; 
     Info<< "hai from myTCSfvOptions\n" << endl; 
     return 0; 
}
