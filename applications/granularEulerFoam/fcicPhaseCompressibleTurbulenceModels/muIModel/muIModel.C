/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
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

#include "muIModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::muIModel::muIModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
  eddyViscosity
  <
  RASModel<EddyDiffusivity<fcicPhaseCompressibleTurbulenceModel>>
  >
  (
   type,
   alpha,
   rho,
   U,
   alphaRhoPhi,
   phi,
   phase,
   propertiesName
   ),
  phase_(phase),
  Hmax_("Hmax", dimLength, coeffDict_),
  L_("L", dimLength, coeffDict_),
  muS_("muS", dimless, coeffDict_),
  mu2_("mu2", dimless, coeffDict_),
  alpha0_(readScalar(coeffDict_.lookup("alpha0"))),
  phiMuI_(readScalar(coeffDict_.lookup("phiMuI"))),
  Cs_(readScalar(coeffDict_.lookup("Cs"))),
  ampl_(readScalar(coeffDict_.lookup("pres_ampl"))),
  beta_(readScalar(coeffDict_.lookup("beta"))),
  theta_(readScalar(coeffDict_.lookup("theta"))),
  g0_
  (
   "g0",
   dimensionSet(1, -1, -2, 0, 0),
   coeffDict_.lookup("g0")
   ),
  ppres_(U.db().lookupObject<volScalarField>("p")),
  g_(U.db().lookupObject<uniformDimensionedVectorField>("g")),
  tau
  (
   IOobject
   (
    "tau",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
   rho_*2*nut_*symm(fvc::grad(U_))*scalar(0.0)
   )
{
  //  nut_ == dimensionedScalar("zero", nut_.dimensions(), 0.0);
  if (type == typeName)
      {
        printCoeffs(type);
      }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::muIModel::~muIModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::muIModel::read()
{
    if
    (
        eddyViscosity
    <
        RASModel<EddyDiffusivity<fcicPhaseCompressibleTurbulenceModel>>
    >::read()
    )
    {
        coeffDict().lookup("Hmax") >> Hmax_;
	coeffDict().lookup("L") >> L_;
	coeffDict().lookup("muS") >> muS_;
	coeffDict().lookup("mu2") >> mu2_;
	coeffDict().lookup("alpha0") >> alpha0_;
	coeffDict().lookup("phiMuI") >> phiMuI_;
	coeffDict().lookup("Cs") >> Cs_;
	coeffDict().lookup("pres_ampl") >> ampl_;
	coeffDict().lookup("beta") >> beta_;
	coeffDict().lookup("theta") >> theta_;
	g0_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::muIModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::muIModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::muIModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
	    -dev(tau)/rho_//(sigmaHyp-tr(sigmaHyp)*symmTensor::I)/rho_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::muIModel::pPrime() const
{
  
  tmp<volScalarField> tpPrime
    (
     g0_
     *min
     (
      exp(ampl_*(alpha_ - alpha0_)),
      scalar(1e3)
      )
     );
  
  volScalarField::Boundary& bpPrime =
    tpPrime.ref().boundaryFieldRef();
  
  forAll(bpPrime, patchi)
    {
      if (!bpPrime[patchi].coupled())
        {
	  bpPrime[patchi] == 0;
        }
    }
  
  return tpPrime;
  

  /*  
  const volScalarField& rho = phase_.rho();
  volScalarField K_("K_",ampl_*rho*mag(g_)*Hmax_);

  //  dimensionedScalar small_p("small_p", dimMass/dimLength/dimTime/dimTime, 0.0);
  
  tmp<volScalarField> tpPrime
    (
     K_/alpha0_//K_*(alpha_/alpha0_-scalar(1))
     );
  
    tmp<volScalarField> tpPrime
    (
     tr(tau)/scalar(3.0)
     );
    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;*/
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::muIModel::pPrimef() const
{
  return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::muIModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
	    -dev(tau)//(sigmaHyp-tr(sigmaHyp)*symmTensor::I)
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::muIModel::divDevRhoReff
(
    volVectorField& U
) const
{

  bool debug=true;

  Info<<"constructing divDevRhoReff"<<endl;
  return
    (
     - fvm::laplacian(rho_*nut_, U)
     - fvc::div(dev(tau))
     );
}

void Foam::RASModels::muIModel::correct()
{

  const volScalarField& rho = phase_.rho();
  const volVectorField& U = U_;
  
  tmp<volScalarField> tda(phase_.d());
  const volScalarField& da = tda();
  
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-6);
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, 1e-6);
  
  volScalarField srnz("srnz", max(sqrt(2.0)*mag(symm(fvc::grad(U_))),smallval_sr));
  volSymmTensorField srnzT("srnzT", symm(fvc::grad(U_)));

  Info<<min(da)<<endl;
  Info<<max(da)<<endl;
  volScalarField eta("eta", Cs_*da*Cs_*da*srnz);
  volScalarField kSps("kSps", (nut_/scalar(0.08)/da)*(nut_/scalar(0.08)/da));
  
  tau = rho_*(2*eta*srnzT - 2/3*kSps*dimensioned<symmTensor>("I", dimless, symmTensor::I));

  volScalarField I0("I0", scalar(5/2)*beta_*da/L_/sqrt(phiMuI_*cos(theta_*M_PI/scalar(180))));
  volScalarField zeta("zeta", I0/sqrt(rho*da*da));


  volScalarField K_("K_",ampl_*rho_*mag(g_)*Hmax_);

  //  dimensionedScalar small_p("small_p", dimMass/dimLength/dimTime/dimTime, 0.0);

  volScalarField pp
    (
     "pp",
     max(K_*(alpha_/alpha0_-scalar(1)),smallval_p)
     );
  
  nut_ = (muS_*pp/srnz + (mu2_-muS_)*pp/(srnz+zeta*sqrt(pp+smallval_p)))/rho_;
  
  bool debug=true;
  
  if(debug)
    {
      //      Info<<min(eta)<<"\n";
      //      Info<<max(eta)<<"\n";
      Info<<min(nut_)<<"\n";
      Info<<max(nut_)<<"\n";
    }
}


// ************************************************************************* //
