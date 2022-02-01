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

#include "fcicGranularPressureModel.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
    defineTypeNameAndDebug(fcicGranularPressureModel, 0);

    defineRunTimeSelectionTable(fcicGranularPressureModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModel::fcicGranularPressureModel
(
 const dictionary& dict,
 const volVectorField& U,
 const volScalarField& alpha,
 const volScalarField& rho
)
:
  dict_(dict),
  U_(U),
  alpha_(alpha),
  rho_(rho)
{}


// * * * * * * * * * * * * * * Member Function * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::biomassModels::fcicGranularPressureModel::strainRate() const
 {
     return sqrt(2.0)*mag(symm(fvc::grad(U_)));
 }                                  //SNA_edit

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModel::~fcicGranularPressureModel()
{}


// ************************************************************************* //
