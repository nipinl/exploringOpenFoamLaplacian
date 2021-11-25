//to write 
const scalarField& diag = TEqn.diag(); 
IOField<scalar> wdiag 
( 
    IOobject 
    ( 
        "diag", 
        mesh.time().timeName(), 
        mesh, 
        IOobject::NO_READ 
    ), 
    diag 
); 

const scalarField& upper = TEqn.upper(); 
IOField<scalar> wupper 
( 
    IOobject 
    ( 
        "upper", 
        mesh.time().timeName(), 
        mesh, 
        IOobject::NO_READ 
    ), 
    upper 
); 
     
const scalarField& lower = TEqn.lower(); 
IOField<scalar> wlower 
( 
    IOobject 
    ( 
        "lower", 
        mesh.time().timeName(), 
        mesh, 
        IOobject::NO_READ 
    ), 
    lower 
); 
     
const scalarField& source = TEqn.source(); 
IOField<scalar> wsource 
( 
    IOobject 
    ( 
        "source", 
        mesh.time().timeName(), 
        mesh, 
        IOobject::NO_READ 
    ), 
    source 
); 
