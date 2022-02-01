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
    Foam::bauyer

Description

SourceFiles
    bauyerPressure.C

\*---------------------------------------------------------------------------*/

#ifndef bauyer_H
#define bauyer_H

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
                           Class bauyer Declaration
\*---------------------------------------------------------------------------*/

class bauyer
:
    public fcicGranularPressureModel
{

  dictionary coeffDict_;
  dimensionedScalar graind_;
  dimensionedScalar B_phi;
  dimensionedScalar rho_p;
  dimensionedScalar alphaMax_;
  
public:

    //- Runtime type information
    TypeName("bauyer");


    // Constructors

        //- Construct from components
  bauyer(
       const dictionary& dict,
       const volVectorField& U,
       const volScalarField& alpha,
       const volScalarField& rho
       );


    //- Destructor
    virtual ~bauyer();


    // Member Functions


        tmp<volScalarField> granularPressureCoeff
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
