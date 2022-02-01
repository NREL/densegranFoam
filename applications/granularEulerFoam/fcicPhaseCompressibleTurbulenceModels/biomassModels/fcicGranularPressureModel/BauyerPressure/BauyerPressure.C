/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
0;95;0c   \\    /   O peration     |
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

#include "BauyerPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
namespace fcicGranularPressureModels
{
    defineTypeNameAndDebug(bauyer, 0);

    addToRunTimeSelectionTable
    (
        fcicGranularPressureModel,
        bauyer,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModels::bauyer::~bauyer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::tmp<Foam::volScalarField>
Foam::biomassModels::fcicGranularPressureModels::bauyer::
granularPressureCoeff
() const
{
  
  dimensionedScalar smallp("smallp", dimMass/dimLength/dimTime/dimTime, 1e-6);
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-6);
  // Schneiderbauer et al. (2012) (SAP model) :: S. Schneiderbauer, A. Aigner, and S. Pirker.
  // "A comprehensive frictional-kinetic model for gas–particle flows: Analysis of fluidized and moving bed regimes."
  // Chemical Engineering Science, 80:279–292, 2012.

  volScalarField srnz_
    (
     "srnz_",
     max(strainRate(),smallval_sr)
     );
  volScalarField granP
    (
     "granP",
     rho_p*graind_*graind_*srnz_*srnz_*B_phi*B_phi*pow(alpha_/max(alphaMax_-alpha_,scalar(1e-6)),2)
     );
  
  
  if(U_.time().outputTime())
      {
	srnz_.write();
	granP.write();
      }

  Info<<"passed to this line"<<endl;
  Info<<max(granP)<<endl;
  Info<<min(granP)<<endl;
  return max(granP,smallp);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModels::bauyer::bauyer
(
 const dictionary& dict,
 const volVectorField& U,
 const volScalarField& alpha,
 const volScalarField& rho
 )
  :
  fcicGranularPressureModel(dict, U, alpha, rho),
  coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
  graind_("graind", dimLength, coeffDict_),
  B_phi("b_phi", dimless, coeffDict_),
  rho_p("rhoP", dimMass/(dimLength*dimLength*dimLength), coeffDict_),
  alphaMax_("alphaMax", dimless, coeffDict_)
  //  srnz_(U.mesh().lookupObject<volScalarField>("srnz"))
{}

bool Foam::biomassModels::fcicGranularPressureModels::bauyer::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");
    
    coeffDict_.lookup("graind") >> graind_;
    coeffDict_.lookup("b_phi") >> B_phi;
    coeffDict_.lookup("rhoP") >> rho_p;
    coeffDict_.lookup("alphaMax") >> alphaMax_;

    return true;
}
// ************************************************************************* //
