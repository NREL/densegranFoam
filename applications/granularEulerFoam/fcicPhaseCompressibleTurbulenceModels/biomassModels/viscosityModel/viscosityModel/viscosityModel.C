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

#include "viscosityModel.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
    defineTypeNameAndDebug(viscosityModel, 0);
    defineRunTimeSelectionTable(viscosityModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::viscosityModel::viscosityModel
(
    const word& name,
    const dictionary& dict,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const volScalarField& alpha,
    const volScalarField& rho,
    const surfaceScalarField& alphaRhoPhi
)
:
    name_(name),
    dict_(dict),
    U_(U),
    phi_(phi),
    alpha_(alpha),
    rho_(rho),
    alphaRhoPhi_(alphaRhoPhi)
{}



 // * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::biomassModels::viscosityModel::strainRate() const
 {
     return sqrt(2.0)*mag(symm(fvc::grad(U_)));
 }                                  //SNA_edit

/*
Foam::tmp<Foam::volSymmTensorField>
Foam::biomassModels::viscosityModel::sigmaDev() const
{
  return rho_*(2.0*nu*symm(fvc::grad(U_)) - (2.0/3.0)*nu*tr(symm(fvc::grad(U_)))*symmTensor::I);
}
*/

/*
bool Foam::biomassModels::viscosityModel::read(const dictionary& dict)
{
    dict_ = dict;

    return true;
    }*/
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::biomassModels::viscosityModel::~viscosityModel()
{}


// ************************************************************************* //

  
  
