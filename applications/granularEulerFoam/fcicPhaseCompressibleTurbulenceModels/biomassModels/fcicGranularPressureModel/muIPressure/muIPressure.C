/*---------------------------------------------------------------------------* \
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

#include "muIPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
namespace fcicGranularPressureModels
{
    defineTypeNameAndDebug(muI, 0);

    addToRunTimeSelectionTable
    (
        fcicGranularPressureModel,
        muI,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModels::muI::~muI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::tmp<Foam::volScalarField>
Foam::biomassModels::fcicGranularPressureModels::muI::
granularPressureCoeff
() const
{
  
  dimensionedScalar smallp("smallp", dimMass/dimLength/dimTime/dimTime, 1e-6);
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-6);
  dimensionedScalar smallal("smallal", dimless, 1e-6);
  dimensionedScalar bigP("bigP", dimMass/dimLength/dimTime/dimTime, 1e6);
  // Schneiderbauer et al. (2012) (SAP model) :: S. Schneiderbauer, A. Aigner, and S. Pirker.
  // "A comprehensive frictional-kinetic model for gas–particle flows: Analysis of fluidized and moving bed regimes."
  // Chemical Engineering Science, 80:279–292, 2012.

  volScalarField srnz_
    (
     "srnz_",
     max(strainRate(),smallval_sr)
     );/*
  volScalarField granP
    (
     "granP",
     rho_*pow(((delfc_*graind_*srnz_*alpha_/in0_)/max((alphaMax_-alpha_),scalar(1.0e-4))),2)
     );*/

  /*
  volScalarField granP
    (
     "granP",
     pow(B_phi_*alpha_/max(alphaMax_-alpha_,smallal),2)*nu_b_*srnz_*rho_b_
     );
  */  
  Info << "delfc_ " << delfc_ << endl;
  Info << "in0_ " << in0_ << endl;
  volScalarField granP
    (
     "granP",
     rhop_*pow(delfc_*graind_*srnz_*alpha_/in0_/max(alphaMax_-alpha_,smallal),2)
     );
  /*  volScalarField granPch
    (
     "granPch",
     rhop_*pow(B_phi_*graind_*srnz_*alpha_/max(alphaMax_-alpha_,smallal),2)
     );
  */
  
  
  // Info << "granPch max " << max(granPch) << endl;
  // Info << "granPch min " << min(granPch) << endl;
  Info << "srnz max " << max(srnz_) << endl;
  Info << "srnz min " << min(srnz_) << endl;
  Info << "granP max " << max(granP) << endl;
  Info << "granP min " << min(granP) <<	endl;
  
  /*
  volScalarField granP
    (
     "granP",
     pow(delfc_/in0_*alpha_ / max(alphaMax_-alpha_, scalar(1e-3)), 2)*nu*rho_*srnz_
     );
  */
  
  if(U_.time().outputTime())
      {
	granP.write();
	srnz_.write();
      }
  
  return min(granP,bigP);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModels::muI::muI
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
  delfc_("delfc", dimless, coeffDict_),
  in0_("in0", dimless, coeffDict_),
  //  B_phi_("B_phi", dimless, coeffDict_),
  //nu_b_("nu_b", dimLength*dimLength/dimTime , coeffDict_),
  //rho_b_("rho_b", dimMass/dimLength/dimLength/dimLength, coeffDict_),
  alphaMax_("alphaMax", dimless, coeffDict_),
  rhop_("rho_p", dimMass/dimLength/dimLength/dimLength, coeffDict_)
  //  srnz_(U.mesh().lookupObject<volScalarField>("srnz"))
{}

bool Foam::biomassModels::fcicGranularPressureModels::muI::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");
    /*
    coeffDict_.lookup("B_phi") >> B_phi_;
    coeffDict_.lookup("nu_b") >> nu_b_;
    coeffDict_.lookup("rho_b") >> rho_b_;
    coeffDict_.lookup("alphaMax") >> alphaMax_;
    */
    
    coeffDict_.lookup("graind") >> graind_;
    coeffDict_.lookup("delfc") >> delfc_;
    coeffDict_.lookup("in0") >> in0_;
    //coeffDict_.lookup("B_phi") >> B_phi_;
    coeffDict_.lookup("alphaMax") >> alphaMax_;
    coeffDict_.lookup("rho_p") >> rhop_;
    
    return true;
}
// ************************************************************************* //
