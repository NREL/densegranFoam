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
  nu0_("nu0", dimViscosity, coeffDict_),
  n_(readScalar(coeffDict_.lookup("n"))),
  ed0_(readScalar(coeffDict_.lookup("ed0"))),
  ec0_(readScalar(coeffDict_.lookup("ec0"))),
  ei0_(readScalar(coeffDict_.lookup("ei0"))),
  alphaHyp_(readScalar(coeffDict_.lookup("alphaHyp"))),
  beta_(readScalar(coeffDict_.lookup("beta"))),
  ThetaC_(readScalar(coeffDict_.lookup("ThetaC"))),
  pres_(U.db().lookupObject<volScalarField>("p_IR")),
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
   U.mesh()
   )
{
  //  nut_ == dimensionedScalar("zero", nut_.dimensions(), 0.0);
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
	coeffDict().lookup("nu0") >> nu0_;
	coeffDict().lookup("n") >> n_;
	coeffDict().lookup("ed0") >> ed0_;
	coeffDict().lookup("ec0") >> ec0_;
	coeffDict().lookup("ei0") >> ei0_;
	coeffDict().lookup("alphaHyp") >> alphaHyp_;
	coeffDict().lookup("beta") >> beta_;
	coeffDict().lookup("ThetaC") >> ThetaC_;

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
	    - (nut_)*dev(twoSymm(fvc::grad(U_)))//	    -dev(sigmaHyp)/rho_//(sigmaHyp-tr(sigmaHyp)*symmTensor::I)/rho_
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::hypoplasticModel::pPrime() const
{
  tmp<volScalarField> tpPrime
    (
     1/3*tr(sigmaHyp)
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
	    - (rho_*nut_)
	    *dev(twoSymm(fvc::grad(U_)))//-dev(sigmaHyp)//(sigmaHyp-tr(sigmaHyp)*symmTensor::I)
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
  Info<<min(dev(sigmaHyp))<<endl;
  Info<<max(dev(sigmaHyp))<<endl;
  Info<<max(sigmaHyp-(scalar(1)/scalar(3))*tr(sigmaHyp)*symmTensor::I)<<endl;
  Info<<min(fvc::div(dev(sigmaHyp)))<<endl;
  Info<<max(fvc::div(dev(sigmaHyp)))<<endl;
  return
    (
     - fvm::laplacian(rho_*nut_, U)
     - fvc::div
     (
      dev(sigmaHyp)
      )
     );
}

void Foam::RASModels::hypoplasticModel::correct()
{
  volTensorField symmGradVel("symmGradVel", scalar(0.5)*(fvc::grad(U_)+fvc::grad(U_)().T()));
  volTensorField sigmahyp("sigmahyp",rho_*nut_*2*dev(symmGradVel));
  Info<<max(sigmahyp)<<endl;
  Info<<min(sigmahyp)<<endl;
  Info<<"sigmahypfunc: 1"<<endl;
  sigmahyp.replace(0, sigmaHyp.component(symmTensor::XX));
  sigmahyp.replace(1, sigmaHyp.component(symmTensor::XY));
  sigmahyp.replace(2, sigmaHyp.component(symmTensor::XZ));
  sigmahyp.replace(3, sigmaHyp.component(symmTensor::XY));
  sigmahyp.replace(4, sigmaHyp.component(symmTensor::YY));
  sigmahyp.replace(5, sigmaHyp.component(symmTensor::YZ));
  sigmahyp.replace(6, sigmaHyp.component(symmTensor::XZ));
  sigmahyp.replace(7, sigmaHyp.component(symmTensor::YZ));
  sigmahyp.replace(8, sigmaHyp.component(symmTensor::ZZ));
  
  volScalarField sigmahyp11("sigmahyp11", sigmahyp.component(tensor::XX));
  volScalarField sigmahyp12("sigmahyp12", sigmahyp.component(tensor::XY));
  volScalarField sigmahyp13("sigmahyp13", sigmahyp.component(tensor::XZ));
  volScalarField sigmahyp21("sigmahyp21", sigmahyp.component(tensor::YX));
  volScalarField sigmahyp22("sigmahyp22", sigmahyp.component(tensor::YY));
  volScalarField sigmahyp23("sigmahyp23", sigmahyp.component(tensor::YZ));
  volScalarField sigmahyp31("sigmahyp31", sigmahyp.component(tensor::ZX));
  volScalarField sigmahyp32("sigmahyp32", sigmahyp.component(tensor::ZY));
  volScalarField sigmahyp33("sigmahyp33", sigmahyp.component(tensor::ZZ));
  
  Info<<"sigmahypfunc: 2"<<endl;
  volTensorField phidivsigmahyp("phidivsigmahyp", fvc::div(phi_,sigmahyp));
  Info<<"sigmahypfunc: 3"<<endl;
  volScalarField phidivsigmahyp11("phidivsigmahyp11", phidivsigmahyp.component(tensor::XX));
  volScalarField phidivsigmahyp12("phidivsigmahyp12", phidivsigmahyp.component(tensor::XY));
  volScalarField phidivsigmahyp13("phidivsigmahyp13", phidivsigmahyp.component(tensor::XZ));
  volScalarField phidivsigmahyp21("phidivsigmahyp21", phidivsigmahyp.component(tensor::YX));
  volScalarField phidivsigmahyp22("phidivsigmahyp22", phidivsigmahyp.component(tensor::YY));
  volScalarField phidivsigmahyp23("phidivsigmahyp23", phidivsigmahyp.component(tensor::YZ));
  volScalarField phidivsigmahyp31("phidivsigmahyp31", phidivsigmahyp.component(tensor::ZX));
  volScalarField phidivsigmahyp32("phidivsigmahyp32", phidivsigmahyp.component(tensor::ZY));
  volScalarField phidivsigmahyp33("phidivsigmahyp33", phidivsigmahyp.component(tensor::ZZ));
  Info<<"sigmahypfunc: 4"<<endl;
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, scalar(1e-6));
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
  dimensionedScalar limhs("limhs", dimMass/dimLength/dimTime/dimTime, scalar(1.0));
  scalar smallVal = scalar(1e-6);
  Info<<"sigmahypfunc: 5"<<endl;
  volScalarField srnz("srnz", sqrt(2.0)*mag(symm(fvc::grad(U_))));
  Info<<"sigmahypfunc: 6"<<endl;
  volTensorField srnzT("srnzT", symmGradVel);
  Info<<max(srnzT)<<endl;
  Info<<min(srnzT)<<endl;
  Info<<"sigmahypfunc: 7"<<endl;
  volScalarField srnzT11("srnzT11", srnzT.component(tensor::XX));
  volScalarField srnzT12("srnzT12", srnzT.component(tensor::XY));
  volScalarField srnzT13("srnzT13", srnzT.component(tensor::XZ));
  volScalarField srnzT21("srnzT21", srnzT.component(tensor::YX));
  volScalarField srnzT22("srnzT22", srnzT.component(tensor::YY));
  volScalarField srnzT23("srnzT23", srnzT.component(tensor::YZ));
  volScalarField srnzT31("srnzT31", srnzT.component(tensor::ZX));
  volScalarField srnzT32("srnzT32", srnzT.component(tensor::ZY));
  volScalarField srnzT33("srnzT33", srnzT.component(tensor::ZZ));
  Info<<"sigmahypfunc: 8"<<endl;
  volTensorField omegaD("omegaD", fvc::grad(U_)-symmGradVel);
  Info<<"sigmahypfunc: 9"<<endl;
  volScalarField omegaD11("omegaD11", omegaD.component(tensor::XX));
  volScalarField omegaD12("omegaD12", omegaD.component(tensor::XY));
  volScalarField omegaD13("omegaD13", omegaD.component(tensor::XZ));
  volScalarField omegaD21("omegaD21", omegaD.component(tensor::YX));
  volScalarField omegaD22("omegaD22", omegaD.component(tensor::YY));
  volScalarField omegaD23("omegaD23", omegaD.component(tensor::YZ));
  volScalarField omegaD31("omegaD31", omegaD.component(tensor::ZX));
  volScalarField omegaD32("omegaD32", omegaD.component(tensor::ZY));
  volScalarField omegaD33("omegaD33", omegaD.component(tensor::ZZ));
  Info<<"sigmahypfunc: 10"<<endl;
  volScalarField trsigmahyp("trsigmahyp", tr(sigmahyp));
  Info<<"sigmahypfunc: 11"<<endl;
  volScalarField trsigmahypinv("trsigmahypinv", scalar(1.0)/max(smallval_p,trsigmahyp));
  Info<<"sigmahypfunc: 12"<<endl;
  volTensorField sigmahypn("sigmahypn", sigmahyp*trsigmahypinv);
  Info<<"sigmahypfunc: 13"<<endl;
  volTensorField sigmahypndev("sigmahypndev", dev(sigmahypn));
  Info<<"sigmahypfunc: 14"<<endl;
  volScalarField sigmahypndev11("sigmahypndev11", sigmahypndev.component(tensor::XX));
  volScalarField sigmahypndev12("sigmahypndev12", sigmahypndev.component(tensor::XY));
  volScalarField sigmahypndev13("sigmahypndev13", sigmahypndev.component(tensor::XZ));
  volScalarField sigmahypndev21("sigmahypndev21", sigmahypndev.component(tensor::YX));
  volScalarField sigmahypndev22("sigmahypndev22", sigmahypndev.component(tensor::YY));
  volScalarField sigmahypndev23("sigmahypndev23", sigmahypndev.component(tensor::YZ));
  volScalarField sigmahypndev31("sigmahypndev31", sigmahypndev.component(tensor::ZX));
  volScalarField sigmahypndev32("sigmahypndev32", sigmahypndev.component(tensor::ZY));
  volScalarField sigmahypndev33("sigmahypndev33", sigmahypndev.component(tensor::ZZ));

  
  
  Info<<"sigmahypfunc: 15"<<endl;
  scalar Theta_ = ThetaC_*M_PI/scalar(180);
  Info<<"sigmahypfunc: 16"<<endl;
  scalar c1_ = sqrt(scalar(3.0)/scalar(8.0))*(scalar(3)-sin(Theta_))/max(sin(Theta_),smallVal);
  Info<<"sigmahypfunc: 17"<<endl;
  scalar c2_ = scalar(3.0)/scalar(8.0)*(scalar(3.0)+sin(Theta_))/max(sin(Theta_),smallVal);
  Info<<"sigmahypfunc: 18"<<endl;

  //  volScalarField powtrsigma1n ("powtrsigma1n", pow(-trsigmahyp/hs_,1-n_));
  Info<<"sigpow1"<<endl;
  volScalarField magtrSig ("magtrSig", sqrt(trsigmahyp*trsigmahyp));
  volScalarField powtrsigma ("powtrsigma", pow(-trsigmahyp/hs_,n_));
  Info<<"sigpow2"<<endl;
  volScalarField powtrsigma1n ("powtrsigma1n", pow(-trsigmahyp/hs_,1-n_));
  Info<<"sigpow3"<<endl;
  /*
  forAll(powtrsigma, ii){
    Info<<powtrsigma[ii]<<"::"<<n_<<"::"<<ii<<endl;
    powtrsigma[ii] = pow(powtrsigma[ii],n_);
    Info<<powtrsigma1n[ii]<<"::"<<1-n_<<"::"<<ii<<endl;
    powtrsigma1n[ii] = pow(powtrsigma1n[ii],1-n_);
  }
  */
  
  volScalarField ei ("ei", ei0_*exp(-powtrsigma));
  Info<<"sigmahypfunc: 19"<<endl;
  volScalarField ed ("ed", ed0_*exp(-powtrsigma));
  Info<<"sigmahypfunc: 20"<<endl;
  volScalarField ec ("ec", ec0_*exp(-powtrsigma));
  Info<<"sigmahypfunc: 21"<<endl;
  
  volScalarField voidRatio("voidRatio", max(min(scalar(1.0)/max(alpha_,smallVal)-scalar(1.0),ec),ed));
  Info<<"sigmahypfunc: 22"<<endl;

  scalar hi = scalar(1.0)/max(c1_*c1_,smallVal)+scalar(1.0)/scalar(3.0)*pow((ei0_+ed0_)/max((ec0_-ed0_),smallVal),alphaHyp_)
    /max(c1_,smallVal)/sqrt(scalar(3.0));
Info<<"sigmahypfunc: 23"<<endl;
 
 volScalarField fs_
    ("fs_",
     hs_/max(n_,scalar(0.1))/max(hi,smallVal)*pow(ei/max(voidRatio,smallVal),beta_)*(scalar(1.0)+ei)/max(ei,smallVal)*powtrsigma1n
     );
  Info<<"sigmahypfunc: 24"<<endl;
  Info<<min((voidRatio-ed)/max((ec-ed)),scalar(1e-6))<<endl;
  volScalarField fd_("fd_", pow((voidRatio-ed)/max((ec-ed),scalar(1e-6)),alphaHyp_));
  Info<<"sigmahypfunc: 25"<<endl;
  volScalarField a1("a1", scalar(1.0)/(c1_ + c2_*mag(sigmahypndev)*(scalar(1.0)+cos(M_PI-scalar(3.0)*Theta_))));
  Info<<"sigmahypfunc: 26"<<endl;
  
  Info<<"constracting Dsigmahyp"<<endl;

  Info<<min(-phidivsigmahyp11)<<endl;
  Info<<min(- omegaD21*sigmahyp12 - omegaD31*sigmahyp13 + omegaD12*sigmahyp21 + omegaD13*sigmahyp31)<<endl;
  Info<<min((a1*fs_*(scalar(6.0)*a1*srnzT11 + sqrt(scalar(2.0))*fd_*sigmahypndev11*srnz))/scalar(6.0))<<endl;
  Info<<min(fs_*sigmahyp11*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 +
                                                            sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 +
                                                            sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33)*
	    pow(trsigmahypinv,2))<<endl;
  Info<<min(pow(trsigmahypinv,2))<<endl;
  Info<<min((a1*fd_*fs_*sigmahyp11*srnz)*trsigmahypinv)<<endl;
  volScalarField Dsigmahyp11("Dsigmahyp11", -phidivsigmahyp11 - omegaD21*sigmahyp12 - omegaD31*sigmahyp13 + omegaD12*sigmahyp21 +
			     omegaD13*sigmahyp31 + (a1*fs_*(scalar(6.0)*a1*srnzT11 + sqrt(scalar(2.0))*fd_*sigmahypndev11*srnz))/
			     scalar(6.0) + (fs_*sigmahyp11*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 +
							    sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 +
							    sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33))*
			     pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp11*srnz)*trsigmahypinv);

  Info<<"constracting Dsigmahyp12"<<endl;
  
  volScalarField Dsigmahyp12("Dsigmahyp12", -phidivsigmahyp12 - omegaD32*sigmahyp13 + omegaD12*(-sigmahyp11 + sigmahyp22) + omegaD13*sigmahyp32 + (a1*fs_*(scalar(6.0)*a1*srnzT12 + sqrt(scalar(2.0))*fd_*sigmahypndev12*srnz))/scalar(6.0) + sigmahyp12*(omegaD11 - omegaD22 + (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

  Info<<"constracting Dsigmahyp13"<<endl;
  
  volScalarField Dsigmahyp13("Dsigmahyp13", -phidivsigmahyp13 - omegaD23*sigmahyp12 + (omegaD11 - omegaD33)*sigmahyp13 + omegaD12*sigmahyp23 + omegaD13*(-sigmahyp11 + sigmahyp33) + (a1*fs_*(scalar(6.0)*a1*srnzT13 + sqrt(scalar(2.0))*fd_*sigmahypndev13*srnz))/scalar(6.0) + (fs_*sigmahyp13*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2));

  Info<<"constracting Dsigmahyp21"<<endl;
  
  volScalarField Dsigmahyp21("Dsigmahyp21", -phidivsigmahyp21 + omegaD21*(sigmahyp11 - sigmahyp22) - omegaD31*sigmahyp23 + omegaD23*sigmahyp31 + (a1*fs_*(scalar(6.0)*a1*srnzT21 + sqrt(scalar(2.0))*fd_*sigmahypndev21*srnz))/scalar(6.0) + sigmahyp21*(-omegaD11 + omegaD22 + (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));
  
  Info<<"constracting Dsigmahyp22"<<endl;
  
  volScalarField Dsigmahyp22("Dsigmahyp22", -phidivsigmahyp22 + omegaD21*sigmahyp12 - omegaD12*sigmahyp21 - omegaD32*sigmahyp23 + omegaD23*sigmahyp32 + (a1*fs_*(scalar(6.0)*a1*srnzT22 + sqrt(scalar(2.0))*fd_*sigmahypndev22*srnz))/scalar(6.0) + (fs_*sigmahyp22* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33))/ pow(trsigmahyp,2) + (a1*fd_*fs_*sigmahyp22*srnz)*trsigmahypinv);
  
  Info<<"constracting Dsigmahyp23"<<endl;
  
  volScalarField Dsigmahyp23("Dsigmahyp23", -phidivsigmahyp23 + omegaD21*sigmahyp13 - omegaD13*sigmahyp21 + (omegaD22 - omegaD33)*sigmahyp23 + omegaD23*(-sigmahyp22 + sigmahyp33) + (a1*fs_*(scalar(6.0)*a1*srnzT23 + sqrt(scalar(2.0))*fd_*sigmahypndev23*srnz))/scalar(6.0) + (fs_*sigmahyp23*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33))*pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp23*srnz)*trsigmahypinv);
  
  
  Info<<"constracting Dsigmahyp31"<<endl;
  
  volScalarField Dsigmahyp31("Dsigmahyp31", -phidivsigmahyp31 + omegaD32*sigmahyp21 + (-omegaD11 + omegaD33)*sigmahyp31 -
			     omegaD21*sigmahyp32 + omegaD31*(sigmahyp11 - sigmahyp33) +
			     (a1*fs_*(scalar(6.0)*a1*srnzT31 + sqrt(scalar(2.0))*fd_*sigmahypndev31*srnz))/scalar(6.0) +
			     (fs_*sigmahyp31* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 +
					       sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 +
					       sigmahyp33*srnzT33 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2));
  
  Info<<"constracting Dsigmahyp32"<<endl;
 
  volScalarField Dsigmahyp32("Dsigmahyp32", -phidivsigmahyp32 + omegaD31*sigmahyp12 - omegaD12*sigmahyp31 +
			     (-omegaD22 + omegaD33)*sigmahyp32 + omegaD32*(sigmahyp22 - sigmahyp33) +
			     (a1*fs_*(scalar(6.0)*a1*srnzT32 + sqrt(scalar(2.0))*fd_*sigmahypndev32*srnz))/scalar(6.0) +
			     (fs_*sigmahyp32* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 + sigmahyp12*srnzT21 +
					       sigmahyp22*srnzT22 + sigmahyp32*srnzT23 + sigmahyp13*srnzT31 + sigmahyp23*srnzT32 +
					       sigmahyp33*srnzT33))/ pow(trsigmahyp,2) + (a1*fd_*fs_*sigmahyp32*srnz)*trsigmahypinv);
 
  Info<<"constracting Dsigmahyp33"<<endl;
  
  volScalarField Dsigmahyp33("Dsigmahyp33", -phidivsigmahyp33 + omegaD31*sigmahyp13 + omegaD32*sigmahyp23 - omegaD13*sigmahyp31 -
			     omegaD23*sigmahyp32 + (a1*fs_*(scalar(6.0)*a1*srnzT33 + sqrt(scalar(2.0))*fd_*sigmahypndev33*srnz))/
			     scalar(6.0) + (fs_*sigmahyp33* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp31*srnzT13 +
							     sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + sigmahyp32*srnzT23 +
							     sigmahyp13*srnzT31 + sigmahyp23*srnzT32 + sigmahyp33*srnzT33))/
			     pow(trsigmahyp,2) + (a1*fd_*fs_*sigmahyp33*srnz)*trsigmahypinv);
  
  
  
 Info<< min(Dsigmahyp11)<<endl;
 Info<< min(Dsigmahyp12)<<endl;
 Info<< min(Dsigmahyp13)<<endl;
 Info<< min(Dsigmahyp21)<<endl;
 Info<< min(Dsigmahyp22)<<endl;
 Info<< min(Dsigmahyp23)<<endl;
 Info<< min(Dsigmahyp31)<<endl;
 Info<< min(Dsigmahyp32)<<endl;
 Info<< min(Dsigmahyp33)<<endl;

 Info<< max(Dsigmahyp11)<<endl;
 Info<< max(Dsigmahyp12)<<endl;
 Info<< max(Dsigmahyp13)<<endl;
 Info<< max(Dsigmahyp21)<<endl;
 Info<< max(Dsigmahyp22)<<endl;
 Info<< max(Dsigmahyp23)<<endl;
 Info<< max(Dsigmahyp31)<<endl;
 Info<< max(Dsigmahyp32)<<endl;
 Info<< max(Dsigmahyp33)<<endl;


 dimensionedScalar zeroSig ("zeroSig", sigmaHyp.dimensions(), scalar(0.0));

 dimensionedScalar dt = U_.mesh().time().deltaT();
 sigmaHyp.replace(0, (sigmaHyp.component(symmTensor::XX)+Dsigmahyp11*dt));
 sigmaHyp.replace(1, (sigmaHyp.component(symmTensor::XY)+Dsigmahyp12*dt));
 sigmaHyp.replace(2, (sigmaHyp.component(symmTensor::XZ)+Dsigmahyp13*dt));
 sigmaHyp.replace(3, (sigmaHyp.component(symmTensor::YY)+Dsigmahyp22*dt));
 sigmaHyp.replace(4, (sigmaHyp.component(symmTensor::YZ)+Dsigmahyp23*dt));
 sigmaHyp.replace(5, (sigmaHyp.component(symmTensor::ZZ)+Dsigmahyp33*dt));

 Info<<"sigmaHyp constructed!!"<<endl;

 
  bool debug = true;
  if (debug)
    {
      Info<<min(sigmaHyp.component(symmTensor::XX))<<endl;
      Info<<max(sigmaHyp.component(symmTensor::XX))<<endl;
      Info<<min(sigmaHyp.component(symmTensor::XY))<<endl;
      Info<<max(sigmaHyp.component(symmTensor::XY))<<endl;
      Info<<min(sigmaHyp.component(symmTensor::XZ))<<endl;
      Info<<max(sigmaHyp.component(symmTensor::XZ))<<endl;
      Info<<min(sigmaHyp.component(symmTensor::YY))<<endl;
      Info<<max(sigmaHyp.component(symmTensor::YY))<<endl;
      Info<<min(sigmaHyp.component(symmTensor::YZ))<<endl;
      Info<<max(sigmaHyp.component(symmTensor::YZ))<<endl;
      Info<<min(sigmaHyp.component(symmTensor::ZZ))<<endl;
      Info<<max(sigmaHyp.component(symmTensor::ZZ))<<endl;
    }

  Info<<"sigmaHyp displayed!!"<<endl;

  //  forAll(nut_,ii){
  nut_ = min(mag(dev(sigmaHyp))/max(srnz,smallval_sr)/rho_, nu0_);
    //  }
}


// ************************************************************************* //
