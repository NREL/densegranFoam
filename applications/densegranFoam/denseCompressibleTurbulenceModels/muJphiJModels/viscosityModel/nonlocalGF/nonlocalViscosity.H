/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::muJphiJModels::viscosityModels::nonlocal

Description

SourceFiles
    nonlocalViscosity.C

\*---------------------------------------------------------------------------*/

#ifndef nonlocal_H
#define nonlocal_H

#include "viscosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace muJphiJModels
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class nonlocal Declaration
\*---------------------------------------------------------------------------*/

class nonlocal
:
    public viscosityModel
{  

public:

    //- Runtime type information
    TypeName("nonlocal");


    // Constructors

        //- Construct from components
  nonlocal(const dictionary& dict,
	   const volVectorField& U,
	   const surfaceScalarField& phi,
	   const volScalarField& alpha
	   );


    //- Destructor
    virtual ~nonlocal();



  // Member functions

  
  virtual tmp<volScalarField> nu
  (
   const volScalarField& muJ,
   const dimensionedScalar& rho1,
   const volScalarField& srnz,
   const volScalarField& ps,
   const volScalarField& da,
   const volScalarField& g_nlgf
   )const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace muJphiJModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
