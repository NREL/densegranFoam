IOdictionary ppProperties
(
 IOobject
 (
  "ppProperties",
  runTime.constant(),
  mesh,
  IOobject::MUST_READ_IF_MODIFIED,
  IOobject::NO_WRITE
  )
 );

const scalar alphaMax(readScalar(ppProperties.lookup("alphaMax")));

const Switch packingLimiter(ppProperties.lookup("packingLimiter"));
