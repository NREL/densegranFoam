/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::nonlocal

Description
    nonlocal-based dispersed-phase particle physical model.

SourceFiles
    nonlocal.C

\*---------------------------------------------------------------------------*/

#ifndef nonlocal_H
#define nonlocal_H

#include "physicalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace physicalModels
{

/*---------------------------------------------------------------------------*\
                           Class nonlocal Declaration
\*---------------------------------------------------------------------------*/

class nonlocal
:
    public physicalModel
{
    // Private data

        //- The nonlocal physical of the phase
  dimensionedScalar rhoC_;
  dimensionedScalar nuC_;
  dimensionedScalar mu1_;
  dimensionedScalar mu2_;
  dimensionedScalar I0_;
  dimensionedScalar phiM_;
  dimensionedScalar S0_;
  dimensionedScalar aMJ_;
  dimensionedScalar phiRLP_;
  dimensionedScalar phiRCP_;
  dimensionedScalar alphaDeltaMin_;
  dimensionedScalar A_;
  dimensionedScalar t0inv_;
  dimensionedScalar numin_;
  dimensionedScalar numax_;
  volScalarField g_nlgf;

public:

    //- Runtime type information
    TypeName("nonlocal");


    // Constructors

        //- Construct from components
        nonlocal
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~nonlocal();


    // Member Functions

  tmp<volScalarField> calcNu() const;
  tmp<volScalarField> calcGranP() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace physicalModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
