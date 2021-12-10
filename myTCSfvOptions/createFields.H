volScalarField T 
( 
   IOobject 
   ( 
      "T", 
      runTime.timeName(), 
      mesh, 
      IOobject::MUST_READ, 
      IOobject::AUTO_WRITE 
   ), 
   mesh 
); 

IOdictionary transportProperties 
( 
   IOobject 
   ( 
      "transportProperties", 
      runTime.constant(), 
      mesh, 
      IOobject::MUST_READ_IF_MODIFIED, 
      IOobject::NO_WRITE 
   ) 
); 

dimensionedScalar DT 
( 
    "DT",
    dimLength*dimLength/dimTime,
	transportProperties 
);

dimensionedScalar su 
( 
    "su",
    dimTemperature/dimTime,
	transportProperties 
);
dimensionedScalar sp 
( 
    "sp",
    dimless/dimTime,
	transportProperties 
);
Info << "DT: " << DT << endl; 
Info << "su: " << su << endl; 
Info << "sp: " << sp << endl; 
