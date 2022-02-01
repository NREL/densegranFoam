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

#include "hypoplasticModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::hypoplasticModel::hypoplasticModel
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
  alphaMax_(readScalar(coeffDict_.lookup("alphaMax"))),
  phiC_(readScalar(coeffDict_.lookup("phiC"))),
  hs_("hs", dimMass/dimLength/dimTime/dimTime, coeffDict_),
  n_(readScalar(coeffDict_.lookup("n"))),
  ed0_(readScalar(coeffDict_.lookup("ed0"))),
  ec0_(readScalar(coeffDict_.lookup("ec0"))),
  ei0_(readScalar(coeffDict_.lookup("ei0"))),
  alphaHyp_(readScalar(coeffDict_.lookup("alphaHyp"))),
  beta_(readScalar(coeffDict_.lookup("beta"))),
  ThetaC_(readScalar(coeffDict_.lookup("ThetaC"))),

/*
  phiC_("phiC", dimless, coeffDict_),
  hs_("hs", dimMass/dimLength/dimTime/dimTime, coeffDict_),
  n_("n", dimless, coeffDict_),
  ed0_("ed0", dimless, coeffDict_),
  ec0_("ec0", dimless, coeffDict_),
  ei0_("ei0", dimless, coeffDict_),
  alphaHyp_("alphaHyp", dimless, coeffDict_),
  beta_("beta", dimless, coeffDict_),
  ThetaC_("ThetaC", dimless, coeffDict_), */
  sigmaHypVecX
  (
  IOobject
  (
  "sigmaHypVecX",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
   U_.mesh()
    ),
  sigmaHypVecY
  (
  IOobject
  (
  "sigmaHypVecY",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
   U_.mesh()
    ),
  sigmaHypVecZ
  (
  IOobject
  (
  "sigmaHypVecZ",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
    U_.mesh()
    ),
  sigmaHyp
  (
  IOobject
  (
  "sigmaHyp",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
    U_.mesh()
    ),
  g0_
  (
  "g0",
    dimensionSet(1, -1, -2, 0, 0),
    coeffDict_.lookup("g0")
    )
{
    nut_ == dimensionedScalar("zero", nut_.dimensions(), 0.0);

    if (type == typeName)
      {
        printCoeffs(type);
      }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::hypoplasticModel::~hypoplasticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::hypoplasticModel::read()
{
    if
    (
        eddyViscosity
    <
        RASModel<EddyDiffusivity<fcicPhaseCompressibleTurbulenceModel>>
    >::read()
    )
    {
        coeffDict().lookup("alphaMax") >> alphaMax_;
	coeffDict().lookup("phiC") >> phiC_;
	coeffDict().lookup("hs") >> hs_;
	coeffDict().lookup("n") >> n_;
	coeffDict().lookup("ed0") >> ed0_;
	coeffDict().lookup("ec0") >> ec0_;
	coeffDict().lookup("ei0") >> ei0_;
	coeffDict().lookup("alphaHyp") >> alphaHyp_;
	coeffDict().lookup("beta") >> beta_;
	coeffDict().lookup("ThetaC") >> ThetaC_;
        g0_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::hypoplasticModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::hypoplasticModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::hypoplasticModel::R() const
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
	    -dev(sigmaHyp)/rho_//(sigmaHyp-tr(sigmaHyp)*symmTensor::I)/rho_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::hypoplasticModel::pPrime() const
{
  tmp<volScalarField> tpPrime
    (
     tr(sigmaHyp)
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
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::hypoplasticModel::pPrimef() const
{
  return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::hypoplasticModel::devRhoReff() const
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
	    -dev(sigmaHyp)//(sigmaHyp-tr(sigmaHyp)*symmTensor::I)
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::hypoplasticModel::divDevRhoReff
(
    volVectorField& U
) const
{
  Info<<"constructing divDevRhoReff"<<endl;
  return
    (
     - fvm::laplacian(rho_*nut_, U)
     - fvc::div
     (
      //(rho_*nut_)*dev2(T(fvc::grad(U)))
      dev(sigmaHyp)
      )
     );
  //return (- fvm::laplacian(rho_*nut_, U) - fvc::div(dev(sigmaHyp)));//fvc::div(sigmaHyp-tr(sigmaHyp)*symmTensor::I));
}


void Foam::RASModels::hypoplasticModel::correct()
{

  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-6);
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, 1e-6);
  dimensionedScalar limhs("limhs", dimMass/dimLength/dimTime/dimTime, 1.0);
      volScalarField srnz
        (
         "srnz",
         //max(strainRate(), smallval_sr)
         sqrt(2.0)*mag(symm(fvc::grad(U_)))
        );
      volSymmTensorField srnzT
        (
         "srnzT",
	 symm(fvc::grad(U_))
	);
      volTensorField omegaD
        (
         "omegaD",
	 (0.5)*(fvc::grad(U_)-symm(fvc::grad(U_)))
        );

      volScalarField voidRatio
	(
	 "voidRatio",
	 scalar(1)/max(alpha_,VSMALL)-1
	 );

      

      sigmaHyp.component(tensor::XX) = sigmaHypVecX.component(vector::X);
      sigmaHyp.component(tensor::YX) = sigmaHypVecX.component(vector::Y);
      sigmaHyp.component(tensor::ZX) = sigmaHypVecX.component(vector::Z);
      sigmaHyp.component(tensor::XY) = sigmaHypVecY.component(vector::X);
      sigmaHyp.component(tensor::YY) = sigmaHypVecY.component(vector::Y);
      sigmaHyp.component(tensor::ZY) = sigmaHypVecY.component(vector::Z);
      sigmaHyp.component(tensor::XZ) = sigmaHypVecZ.component(vector::X);
      sigmaHyp.component(tensor::YZ) = sigmaHypVecZ.component(vector::Y);
      sigmaHyp.component(tensor::ZZ) = sigmaHypVecZ.component(vector::Z);

      volScalarField trSigmaInv
	("trSigmaInv", 1/max(smallval_p,tr(sigmaHyp)));
      

      //      dimensionedScalar trSigmaInv("trSigmaInv", dimMass/dimLength/dimTime/dimTime,tr(sigmaHyp));
      volSymmTensorField sigmaHypHat
	(
	 "sigmaHypHat",
	 sigmaHyp*trSigmaInv
	 );

      
      volSymmTensorField tauDev
	(
	 "tauDev",
	 dev(sigmaHypHat)
	 );

      vector ii(1,0,0);
      vector jj(0,1,0);
      vector kk(0,0,1);
      

      volVectorField tauDevX
        (
         "tauDevX",
         tauDev.component(tensor::XX)*ii + tauDev.component(tensor::YX)*jj + tauDev.component(tensor::ZX)*kk
         );
      volVectorField tauDevY
	(
         "tauDevY",
         (tauDev.component(tensor::XY)*ii + tauDev.component(tensor::YY)*jj + tauDev.component(tensor::ZY)*kk)
         );
      volVectorField tauDevZ
        (
         "tauDevZ",
         (tauDev.component(tensor::XZ)*ii + tauDev.component(tensor::YZ)*jj + tauDev.component(tensor::ZZ)*kk)
         );

      volVectorField srnzTX
        (
         "srnzTX",
         srnzT.component(tensor::XX)*ii + srnzT.component(tensor::YX)*jj + srnzT.component(tensor::ZX)*kk
         );
      volVectorField srnzTY
        (
         "srnzTY",
	 (srnzT.component(tensor::XY)*ii + srnzT.component(tensor::YY)*jj + srnzT.component(tensor::ZY)*kk)
         );
      volVectorField srnzTZ
        (
         "srnzTZ",
         (srnzT.component(tensor::XZ)*ii + srnzT.component(tensor::YZ)*jj + srnzT.component(tensor::ZZ)*kk)
         );

      volVectorField omegaDX
	(
         "omegaDX",
         (omegaD.component(tensor::XX)*ii + omegaD.component(tensor::YX)*jj + omegaD.component(tensor::ZX)*kk)
         );
      volVectorField omegaDY
        (
         "omegaDY",
         (omegaD.component(tensor::XY)*ii + omegaD.component(tensor::YY)*jj + omegaD.component(tensor::ZY)*kk)
         );
      volVectorField omegaDZ
        (
         "omegaDZ",
         (omegaD.component(tensor::XZ)*ii + omegaD.component(tensor::YZ)*jj + omegaD.component(tensor::ZZ)*kk)
         );

      
      Info<<"ll:1"<<endl;
      scalar Theta_ = ThetaC_*M_PI/scalar(180);
      Info<<"ll:2"<<endl;
      scalar c1_ = sqrt(scalar(3.0)/scalar(8.0))*(scalar(3)-sin(Theta_))/max(sin(Theta_),VSMALL);
      Info<<"ll:3"<<endl;
      scalar c2_ = scalar(3.0)/scalar(8.0)*(3+sin(Theta_))/max(sin(Theta_),VSMALL);
      Info<<"ll:4"<<endl;
      volScalarField ei ("ei", ei0_*exp(-pow(-tr(sigmaHyp)/max(hs_,smallval_p),1-n_)));Info<<"ll:5"<<endl;
      volScalarField ed ("ed", ed0_*exp(-pow(-tr(sigmaHyp)/max(hs_,smallval_p),1-n_)));Info<<"ll:6"<<endl;
      volScalarField ec ("ec", ec0_*exp(-pow(-tr(sigmaHyp)/max(hs_,smallval_p),1-n_)));Info<<"ll:7"<<endl;
      

      scalar hi = 1/max(c1_*c1_,VSMALL)+scalar(1)/scalar(3)*pow((ei0_+ed0_)/max((ec0_-ed0_),VSMALL),alphaHyp_)/max(c1_,VSMALL)/sqrt(scalar(3));Info<<"ll:8"<<endl;
      volScalarField fs_
	(
	 "fs_",
	 min(hs_/n_/max(hi,VSMALL)*pow(ei/max(voidRatio,VSMALL),beta_)*(1+ei)/max(ei,VSMALL)*pow(-tr(sigmaHyp)/max(hs_,smallval_p),1-n_),limhs)
	 );Info<<"max:"<<max(fs_)<<" min: "<<max(fs_)<<endl;Info<<"ll:13"<<endl;
      volScalarField fd_
        (
         "fd_",
	 min(pow((voidRatio-ed)/max((ec-ed),VSMALL),alphaHyp_),scalar(1.0))
         );
      
      volScalarField a1
        (
         "a1",
	 1/max((c1_ + c2_*mag(tauDev)*(1+cos(M_PI-3*Theta_))),VSMALL)
	 );Info<<"ll:11"<<endl;

      volVectorField omXsigX = (omegaDX^sigmaHypVecX) - (sigmaHypVecX^omegaDX);
      volVectorField omYsigY = (omegaDY^sigmaHypVecY) - (sigmaHypVecY^omegaDY);
      volVectorField omZsigZ = (omegaDZ^sigmaHypVecZ) - (sigmaHypVecZ^omegaDZ);

      volVectorField sigmaHypVecXHat
        (
         "sigmaHypVecXHat",
	 sigmaHypVecX*trSigmaInv
	 );
      volVectorField sigmaHypVecYHat
        (
         "sigmaHypVecYHat",
         sigmaHypVecY*trSigmaInv
         );
      volVectorField sigmaHypVecZHat
        (
         "sigmaHypVecZHat",
         sigmaHypVecZ*trSigmaInv
         );
      Info<<"ll:12"<<endl;
            Info<<"max:"<<max(fs_)<<" min: "<<max(fs_)<<endl;Info<<"ll:13"<<endl;
            Info<<"max:"<<max(fd_)<<" min: "<<max(fd_)<<endl;Info<<"ll:14"<<endl;
            Info<<"max:"<<max(a1)<<" min: "<<max(a1)<<endl;Info<<"ll:15"<<endl;
      //Info<<"max:"<<max(trSigmaInv)<<" min: "<<min(trSigmaInv)<<endl;Info<<"ll:16"<<endl;

      
      Info<<"solving for sigmaHypVecXEqn"<<endl;
      fvVectorMatrix sigmaHypVecXEqn
      (fvm::ddt(sigmaHypVecX)
       + fvm::div(phi_, sigmaHypVecX) + fvm::SuSp(-fvc::div(phi_), sigmaHypVecX)
       ==
       omXsigX
       + fs_*a1*a1*srnzTX
       + fvc::Sp(fs_*tr(sigmaHypVecX * srnzTX)*trSigmaInv*trSigmaInv, sigmaHypVecX)
       + fvc::Sp(fs_*fd_*a1*srnz*trSigmaInv, sigmaHypVecX)
       + fs_*fd_*a1/3/sqrt(2.0)*srnz * tauDevX
       );

      Info<<"solving for sigmaHypVecYEqn"<<endl;
      fvVectorMatrix sigmaHypVecYEqn
      (fvm::ddt(sigmaHypVecY)
       + fvm::div(phi_, sigmaHypVecY) + fvm::SuSp(-fvc::div(phi_), sigmaHypVecY)
       ==
       omYsigY
       + fs_*a1*a1*srnzTY
       + fvc::Sp(fs_*tr(sigmaHypVecY * srnzTY)*trSigmaInv*trSigmaInv, sigmaHypVecY)
       + fvc::Sp(fs_*fd_*a1*srnz*trSigmaInv, sigmaHypVecY)
       + fs_*fd_*a1/3/sqrt(2.0)*srnz * tauDevY
       );

      Info<<"solving for sigmaHypVecZEqn"<<endl;
      fvVectorMatrix sigmaHypVecZEqn
      (fvm::ddt(sigmaHypVecZ)
       + fvm::div(phi_, sigmaHypVecZ) + fvm::SuSp(-fvc::div(phi_), sigmaHypVecZ)
       ==
       omZsigZ
       + fs_*a1*a1*srnzTZ
       + fvc::Sp(fs_*tr(sigmaHypVecZ * srnzTZ)*trSigmaInv*trSigmaInv, sigmaHypVecZ)
       + fvc::Sp(fs_*fd_*a1*srnz*trSigmaInv, sigmaHypVecZ)
       + fs_*fd_*a1/3/sqrt(2.0)*srnz * tauDevZ
       );
      
      /*
      Info<<"solving for sigmaHypVecZEqn"<<endl;
      fvVectorMatrix sigmaHypVecZEqn
      (fvm::ddt(sigmaHypVecZ)
       + fvm::div(phi_, sigmaHypVecZ) + fvm::SuSp(-fvc::div(phi_), sigmaHypVecZ)
       ==
       omZsigZ
       + fs_*a1*a1*srnzTZ
       + fs_*sigmaHypVecZHat*tr(sigmaHypVecZHat*srnzTZ)
       + fs_*fd_*a1*sigmaHypVecZHat*srnz
       + fs_*fd_*a1/3/sqrt(2.0)*srnz*tauDevZ
       );
      */

      sigmaHypVecXEqn.relax();
      sigmaHypVecXEqn.solve();
      Info<<"equation 1 solved"<<endl;
      sigmaHypVecYEqn.relax();
      sigmaHypVecYEqn.solve();
      Info<<"equation 2 solved"<<endl;
      sigmaHypVecZEqn.relax();
      sigmaHypVecZEqn.solve();
      Info<<"equations solved"<<endl;
      
      sigmaHyp.component(tensor::XX) = sigmaHypVecX.component(vector::X);
      sigmaHyp.component(tensor::YX) = sigmaHypVecX.component(vector::Y);
      sigmaHyp.component(tensor::ZX) = sigmaHypVecX.component(vector::Z);
      sigmaHyp.component(tensor::XY) = sigmaHypVecY.component(vector::X);
      sigmaHyp.component(tensor::YY) = sigmaHypVecY.component(vector::Y);
      sigmaHyp.component(tensor::ZY) = sigmaHypVecY.component(vector::Z);
      sigmaHyp.component(tensor::XZ) = sigmaHypVecZ.component(vector::X);
      sigmaHyp.component(tensor::YZ) = sigmaHypVecZ.component(vector::Y);
      sigmaHyp.component(tensor::ZZ) = sigmaHypVecZ.component(vector::Z);
      
      bool debug=true;
      
      if(debug)
	{
	  Info<<min(sigmaHypVecX)<<"\n";
	  Info<<max(sigmaHypVecX)<<"\n";
	  Info<<min(sigmaHypVecY)<<"\n";
	  Info<<max(sigmaHypVecY)<<"\n";
	  Info<<min(sigmaHypVecZ)<<"\n";
	  Info<<max(sigmaHypVecZ)<<"\n";
	}
}


// ************************************************************************* //
