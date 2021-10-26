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

#include "nonlocalViscosity.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace muJphiJModels
{
namespace viscosityModels
{
    defineTypeNameAndDebug(nonlocal, 0);
    addToRunTimeSelectionTable(viscosityModel, nonlocal, dictionary);
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::muJphiJModels::viscosityModels::nonlocal::nonlocal
(
 const dictionary& dict,
 const volVectorField& U,
 const surfaceScalarField& phi,
 const volScalarField& alpha
)
:
  viscosityModel(dict,U,phi,alpha)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::muJphiJModels::viscosityModels::nonlocal::~nonlocal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::muJphiJModels::viscosityModels::nonlocal::nu
(
 const volScalarField& muJ,
 const dimensionedScalar& rho1,
 const volScalarField& srnz,
 const volScalarField& ps,
 const volScalarField& da,
 const volScalarField& g_nlgf
 ) const
{ 
  return ps/g_nlgf/rho1;  
}

// ************************************************************************* //
