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

\*---------------------------------------------------------------------------*/

#include "constantViscosity.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace physicalModels
{
    defineTypeNameAndDebug(constant, 0);

    addToRunTimeSelectionTable
    (
        physicalModel,
        constant,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicalModels::constant::constant
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    physicalModel(dict, phase),
    nu0_("nu0", dimensionSet(0, 2, -1, 0, 0), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::physicalModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::physicalModels::constant::calcNu() const
{
    return tmp<Foam::volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nu",
                phase_.U().time().timeName(),
                phase_.U().mesh()
            ),
            phase_.U().mesh(),
            nu0_
        )
    );
}
Foam::tmp<Foam::volScalarField> Foam::physicalModels::constant::calcGranP() const
{
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
  volScalarField pS_ = max(phase_.mesh().lookupObject<volScalarField>("p"), smallval_p);
  return scalar(0.0)*pS_;
}

// ************************************************************************* //
