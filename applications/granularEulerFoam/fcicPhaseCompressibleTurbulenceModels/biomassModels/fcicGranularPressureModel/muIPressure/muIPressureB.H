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
    Foam::muI

Description

SourceFiles
    muIPressure.C

\*---------------------------------------------------------------------------*/

#ifndef muI_H
#define muI_H

#include "fcicGranularPressureModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
namespace fcicGranularPressureModels
{

/*---------------------------------------------------------------------------*\
                           Class muI Declaration
\*---------------------------------------------------------------------------*/

class muI
:
    public fcicGranularPressureModel
{

  dictionary coeffDict_;
  dimensionedScalar B_phi_;
  //  dimensionedScalar nu_b_;
  //  dimensionedScalar rho_b_;
  dimensionedScalar graind_;
  //  dimensionedScalar delfc_;
  //  dimensionedScalar in0_;
  dimensionedScalar alphaMax_;
  dimensionedScalar rhop_;
  
public:

    //- Runtime type information
    TypeName("muI");


    // Constructors

        //- Construct from components
  muI(
       const dictionary& dict,
       const volVectorField& U,
       const volScalarField& alpha,
       const volScalarField& rho
       );


    //- Destructor
    virtual ~muI();


    // Member Functions


        tmp<volScalarField> granularPressureCoeffPrime
        () const;

  virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularPressureModels
} // End namespace biomassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
