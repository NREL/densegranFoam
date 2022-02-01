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
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"
#include "fvOptions.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace physicalModels
{
    defineTypeNameAndDebug(nonlocal, 0);

    addToRunTimeSelectionTable
    (
        physicalModel,
        nonlocal,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicalModels::nonlocal::nonlocal
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    physicalModel(dict, phase),
    rhoC_("rhoC",dimensionSet(1, -3, 0, 0, 0), dict),
    nuC_("nuC",dimensionSet(0, 2, -1, 0, 0), dict),
    mu1_("mu1",dimensionSet(0, 0, 0, 0, 0), dict),
    mu2_("mu2",dimensionSet(0, 0, 0, 0, 0), dict),
    I0_("I0",dimensionSet(0, 0, 0, 0, 0), dict),
    phiM_("phiM",dimensionSet(0, 0, 0, 0, 0), dict),
    S0_("S0",dimensionSet(0, 0, -1, 0, 0), dict),
    aMJ_("aMJ",dimensionSet(1, -1, -2, 0, 0), dict),
    phiRLP_("phiRLP",dimensionSet(0, 0, 0, 0, 0), dict),
    phiRCP_("phiRCP",dimensionSet(0, 0, 0, 0, 0), dict),
    alphaDeltaMin_("alphaDeltaMin",dimensionSet(0, 0, 0, 0, 0), dict),
    A_("A",dimensionSet(0, 0, 0, 0, 0), dict),
    t0inv_("t0inv",dimensionSet(0, 0, -1, 0, 0), dict),
    numax_("numax",dimensionSet(0, 2, -1, 0, 0), dict),
    numin_("numin",dimensionSet(0, 2, -1, 0, 0), dict),
    g_nlgf
    (
     IOobject
     (
      "g_nlgf",
      phase_.time().timeName(),
      phase_.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     phase_.mesh()
     )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::physicalModels::nonlocal::~nonlocal()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::physicalModels::nonlocal::calcNu() const
{

      // Local references
  const volScalarField& alpha_ = phase_;
  const surfaceScalarField& phi_ = phase_.phi(); 
  const volScalarField& da = phase_.d();
  const dimensionedScalar& rhoS_ = phase_.rho();
  const volVectorField& U = phase_.U();    
  tmp<volTensorField> tgradU(fvc::grad(U));
  const volTensorField& gradU(tgradU());
  volSymmTensorField D(symm(gradU));
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, scalar(1e-6));
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
  
  
  volScalarField srnz("srnz", min(mag(D)-tr(gradU)/scalar(3),S0_));
  const fvPatchList& patches = U.mesh().boundary();

  volScalarField pS_ = max(phase_.mesh().lookupObject<volScalarField>("p"), smallval_p);
  /*
  volScalarField pS_ ("pS_", max(
				 aMJ_*max(alpha_-phiRLP_,scalar(0))/
				 max((phiRCP_-alpha_),alphaDeltaMin_)+
				 scalar(2)*nuC_*rhoC_*srnz/
				 pow(phiPM/max(alpha_-scalar(1),alphaDeltaMin_),2),smallval_p));
  

  volScalarField J1_ = scalar(2)*srnz*nuC_*rhoC_/max(pS_,smallval_p);
  */
  volScalarField muI_ = srnz/g_nlgf;
  
  volScalarField pS_mjpj
    ("pS_mjpj", pS_);
  
  if(U.time().outputTime())
    {
      pS_mjpj.write();
      srnz.write();
    }
  dimensionedScalar b_nlgf("b_nlgf", (mu2_-mu1_)/I0_);
  

  fvScalarMatrix gEqn
    (fvm::ddt(g_nlgf) + fvm::div(phi_, g_nlgf) - fvm::SuSp(fvc::div(phi_), g_nlgf)
     - fvm::laplacian(A_*A_*da*da*t0inv_, g_nlgf) ==
     fvm::Sp(-t0inv_*(mu2_-mu1_)*(mu1_ - muI_)/(mu2_ - muI_), g_nlgf)
     + fvm::Sp(-t0inv_*b_nlgf*muI_* sqrt(rhoS_*da*da/pS_)*srnz/max(muI_,scalar(1e-6)), g_nlgf)
     );
  /*
    fvScalarMatrix gEqn
    (fvm::ddt(g_nlgf) + fvm::div(phi_, g_nlgf) - fvm::SuSp(fvc::div(phi_), g_nlgf)
    - fvm::laplacian(A_*A_*da*da*t0inv_, g_nlgf) ==
    fvm::Sp(-t0inv_*(mu2_-mu1_)*(mu1_*g_nlgf - srnz)/(mu2_*g_nlgf - srnz), g_nlgf)
    + fvm::Sp(-t0inv_*b_nlgf*da*srnz * sqrt(alpha_*rhoS_/pS_), g_nlgf)
    );
  */
  gEqn.relax();
  gEqn.solve();
  //  g_nlgf.max(1e-6);
  //  g_nlgf.min(1e6);
  //  g_nlgf.correctBoundaryConditions();
  
    
  
  return max(min(pS_/g_nlgf,numax_),numin_);
}
Foam::tmp<Foam::volScalarField> Foam::physicalModels::nonlocal::calcGranP() const
{
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
  volScalarField pS_ = max(phase_.mesh().lookupObject<volScalarField>("p"), smallval_p);
  return scalar(0.0)*pS_;
}

// ************************************************************************* //
