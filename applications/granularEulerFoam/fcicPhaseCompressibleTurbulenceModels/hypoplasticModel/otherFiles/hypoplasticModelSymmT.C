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
  Info<<min(dev(sigmaHyp))<<endl;
  Info<<max(dev(sigmaHyp))<<endl;
  Info<<max(sigmaHyp-(1/3)*tr(sigmaHyp)*symmTensor::I)<<endl;
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
  volSymmTensorField sigmahyp("sigmahyp", rho_*nut_*2*dev(symm(fvc::grad(U_))));
  Info<<"sigmahypfunc: 1"<<endl;
  sigmahyp.replace(0, sigmaHyp.component(tensor::XX));
  sigmahyp.replace(1, sigmaHyp.component(tensor::XY));
  sigmahyp.replace(2, sigmaHyp.component(tensor::XZ));
  sigmahyp.replace(3, sigmaHyp.component(tensor::YX));
  sigmahyp.replace(4, sigmaHyp.component(tensor::YY));
  sigmahyp.replace(5, sigmaHyp.component(tensor::YZ));
  volScalarField sigmahyp11("sigmahyp11", sigmaHyp.component(tensor::XX));
  volScalarField sigmahyp12("sigmahyp12", sigmaHyp.component(tensor::XY));
  volScalarField sigmahyp13("sigmahyp13", sigmaHyp.component(tensor::XZ));
  volScalarField sigmahyp21("sigmahyp21", sigmaHyp.component(tensor::YX));
  volScalarField sigmahyp22("sigmahyp22", sigmaHyp.component(tensor::YY));
  volScalarField sigmahyp23("sigmahyp23", sigmaHyp.component(tensor::YZ));
  Info<<"sigmahypfunc: 2"<<endl;
  volSymmTensorField phidivsigmahyp("phidivsigmahyp", fvc::div(phi_,sigmahyp));
  Info<<"sigmahypfunc: 3"<<endl;
  volScalarField phidivsigmahyp11("phidivsigmahyp11", phidivsigmahyp.component(tensor::XX));
  volScalarField phidivsigmahyp12("phidivsigmahyp12", phidivsigmahyp.component(tensor::XY));
  volScalarField phidivsigmahyp13("phidivsigmahyp13", phidivsigmahyp.component(tensor::XZ));
  volScalarField phidivsigmahyp21("phidivsigmahyp21", phidivsigmahyp.component(tensor::YX));
  volScalarField phidivsigmahyp22("phidivsigmahyp22", phidivsigmahyp.component(tensor::YY));
  volScalarField phidivsigmahyp23("phidivsigmahyp23", phidivsigmahyp.component(tensor::YZ));
  Info<<"sigmahypfunc: 4"<<endl;
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-6);
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, 1e-6);
  dimensionedScalar limhs("limhs", dimMass/dimLength/dimTime/dimTime, 1.0);
  scalar smallVal = scalar(1e-6);
  Info<<"sigmahypfunc: 5"<<endl;
  volScalarField srnz("srnz", sqrt(2.0)*mag(symm(fvc::grad(U_))));
  Info<<"sigmahypfunc: 6"<<endl;
  volSymmTensorField srnzT("srnzT", symm(fvc::grad(U_)));
  volTensorField srnzTen("srnzTen", scalar(0.5)*(fvc::grad(U_)+fvc::grad(U_)().T()));
  Info<<max(srnzTen)<<endl;
  Info<<"sigmahypfunc: 7"<<endl;
  volScalarField srnzT11("srnzT11", srnzT.component(tensor::XX));
  volScalarField srnzT12("srnzT12", srnzT.component(tensor::XY));
  volScalarField srnzT13("srnzT13", srnzT.component(tensor::XZ));
  volScalarField srnzT21("srnzT21", srnzT.component(tensor::YX));
  volScalarField srnzT22("srnzT22", srnzT.component(tensor::YY));
  volScalarField srnzT23("srnzT23", srnzT.component(tensor::YZ));
  Info<<"sigmahypfunc: 8"<<endl;
  volTensorField omegaD("omegaD", fvc::grad(U_)-symm(fvc::grad(U_)));
  Info<<"sigmahypfunc: 9"<<endl;
  Info<<max(omegaD)<<endl;
  Info<<min(omegaD)<<endl;
  volScalarField omegaD11("omegaD11", omegaD.component(tensor::XX));
  volScalarField omegaD12("omegaD12", omegaD.component(tensor::XY));
  volScalarField omegaD13("omegaD13", omegaD.component(tensor::XZ));
  volScalarField omegaD21("omegaD21", omegaD.component(tensor::YX));
  volScalarField omegaD22("omegaD22", omegaD.component(tensor::YY));
  volScalarField omegaD23("omegaD23", omegaD.component(tensor::YZ));
  Info<<"sigmahypfunc: 10"<<endl;
  volScalarField trsigmahyp("trsigmahyp", tr(sigmahyp));
  Info<<"sigmahypfunc: 11"<<endl;
  volScalarField trsigmahypinv("trsigmahypinv", 1/max(smallval_p,tr(sigmahyp)));
  Info<<"sigmahypfunc: 12"<<endl;
  volSymmTensorField sigmahypn("sigmahypn", sigmahyp*trsigmahypinv);
  Info<<"sigmahypfunc: 13"<<endl;
  volSymmTensorField sigmahypndev("sigmahypndev", dev(sigmahypn));
  Info<<"sigmahypfunc: 14"<<endl;
  volScalarField sigmahypndev11("sigmahypndev11", sigmahypndev.component(tensor::XX));
  volScalarField sigmahypndev12("sigmahypndev12", sigmahypndev.component(tensor::XY));
  volScalarField sigmahypndev13("sigmahypndev13", sigmahypndev.component(tensor::XZ));
  volScalarField sigmahypndev21("sigmahypndev21", sigmahypndev.component(tensor::YX));
  volScalarField sigmahypndev22("sigmahypndev22", sigmahypndev.component(tensor::YY));
  volScalarField sigmahypndev23("sigmahypndev23", sigmahypndev.component(tensor::YZ));

  
  
  Info<<"sigmahypfunc: 15"<<endl;
  scalar Theta_ = ThetaC_*M_PI/scalar(180);
  Info<<"sigmahypfunc: 16"<<endl;
  scalar c1_ = sqrt(scalar(3.0)/scalar(8.0))*(scalar(3)-sin(Theta_))/max(sin(Theta_),smallVal);
  Info<<"sigmahypfunc: 17"<<endl;
  scalar c2_ = scalar(3.0)/scalar(8.0)*(3+sin(Theta_))/max(sin(Theta_),smallVal);
  Info<<"sigmahypfunc: 18"<<endl;

  //  volScalarField powtrsigma1n ("powtrsigma1n", pow(-trsigmahyp/hs_,1-n_));
  Info<<"sigpow1"<<endl;
  volScalarField powtrsigma ("powtrsigma", pow(-3*pres_/hs_,n_));
  Info<<"sigpow2"<<endl;
  volScalarField powtrsigma1n ("powtrsigma1n", pow(-3*pres_/hs_,1-n_));
  Info<<"sigpow3"<<endl;
  
  
  volScalarField ei ("ei", ei0_*exp(-powtrsigma));
  Info<<"sigmahypfunc: 19"<<endl;
  volScalarField ed ("ed", ed0_*exp(-powtrsigma));
  Info<<"sigmahypfunc: 20"<<endl;
  volScalarField ec ("ec", ec0_*exp(-powtrsigma));
  Info<<"sigmahypfunc: 21"<<endl;
  
  volScalarField voidRatio("voidRatio", max(min(scalar(1.0)/max(alpha_,smallVal)-scalar(1.0),ec),ed));
  Info<<"sigmahypfunc: 22"<<endl;

  scalar hi = 1/max(c1_*c1_,smallVal)+scalar(1)/scalar(3)*pow((ei0_+ed0_)/max((ec0_-ed0_),smallVal),alphaHyp_)
    /max(c1_,smallVal)/sqrt(scalar(3));
Info<<"sigmahypfunc: 23"<<endl;
 
 volScalarField fs_
    ("fs_",
     hs_/max(n_,scalar(1e-1))/max(hi,smallVal)*pow(ei/max(voidRatio,smallVal),beta_)*(1+ei)/max(ei,smallVal)*powtrsigma1n
     );
  Info<<"sigmahypfunc: 24"<<endl;
  volScalarField fd_("fd_", pow((voidRatio-ed)/(ec-ed),alphaHyp_));
  Info<<"sigmahypfunc: 25"<<endl;
  volScalarField a1("a1", 1/(c1_ + c2_*mag(sigmahypndev)*(1+cos(M_PI-3*Theta_))));
  Info<<"sigmahypfunc: 26"<<endl;
  
  Info<<"constracting Dsigmahyp"<<endl;

  volScalarField Dsigmahyp11("Dsigmahyp11", -phidivsigmahyp11 +
                             (a1*fs_*(scalar(6.0)*a1*srnzT11 + sqrt(scalar(2.0))*fd_*sigmahypndev11*srnz))/scalar(6.0) +
                             (fs_*sigmahyp11*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))*
                             pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp11*srnz)*trsigmahypinv);

  Info<<"constracting Dsigmahyp12"<<endl;
  
  volScalarField Dsigmahyp12("Dsigmahyp12",  -phidivsigmahyp12 +
                             (a1*fs_*(scalar(6.0)*a1*srnzT12 + sqrt(scalar(2.0))*fd_*sigmahypndev12*srnz))/scalar(6.0) +
                             sigmahyp12*(omegaD11 - omegaD22 +
                                         (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 +
                                               sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

  Info<<"constracting Dsigmahyp13"<<endl;
  
  volScalarField Dsigmahyp13("Dsigmahyp13", -phidivsigmahyp13 +
                             (a1*fs_*(scalar(6.0)*a1*srnzT13 + sqrt(scalar(2.0))*fd_*sigmahypndev13*srnz))/scalar(6.0) +
                             (fs_*sigmahyp13*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 +
                                              a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2));

  Info<<"constracting Dsigmahyp21"<<endl;
  
  volScalarField Dsigmahyp21("Dsigmahyp21", -phidivsigmahyp21 +
                             (a1*fs_*(scalar(6.0)*a1*srnzT21 + sqrt(scalar(2.0))*fd_*sigmahypndev21*srnz))/scalar(6.0) +
                             sigmahyp21*(-omegaD11 + omegaD22 +
                                         (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 +
                                               sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

  Info<<"constracting Dsigmahyp22"<<endl;
  
 volScalarField Dsigmahyp22("Dsigmahyp22", -phidivsigmahyp22 
                            + (a1*fs_*(scalar(6.0)*a1*srnzT22 + sqrt(scalar(2.0))*fd_*sigmahypndev22*srnz))/scalar(6.0)
                            + (fs_*sigmahyp22* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))
                            *pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp22*srnz)*trsigmahypinv);

 Info<<"constracting Dsigmahyp23"<<endl;
 
 volScalarField Dsigmahyp23("Dsigmahyp23", -phidivsigmahyp23 +
                            (a1*fs_*(scalar(6.0)*a1*srnzT23 + sqrt(scalar(2.0))*fd_*sigmahypndev23*srnz))/scalar(6.0) +
                            (fs_*sigmahyp23*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 +
                                             sigmahyp22*srnzT22))*pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp23*srnz)*trsigmahypinv);

 

  

  /*

  
  volScalarField Dsigmahyp11("Dsigmahyp11", -phidivsigmahyp11 - omegaD21*sigmahyp12 + omegaD12*sigmahyp21 +
			     (a1*fs_*(scalar(6.0)*a1*srnzT11 + sqrt(scalar(2.0))*fd_*sigmahypndev11*srnz))/scalar(6.0) +
			     (fs_*sigmahyp11*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))*
			     pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp11*srnz)*trsigmahypinv);
  
  volScalarField Dsigmahyp12("Dsigmahyp12",  -phidivsigmahyp12 + omegaD12*(-sigmahyp11 + sigmahyp22) +
			     (a1*fs_*(scalar(6.0)*a1*srnzT12 + sqrt(scalar(2.0))*fd_*sigmahypndev12*srnz))/scalar(6.0) +
			     sigmahyp12*(omegaD11 - omegaD22 +
					 (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 +
					       sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

  volScalarField Dsigmahyp13("Dsigmahyp13", -phidivsigmahyp13 - omegaD23*sigmahyp12 + omegaD11*sigmahyp13 + omegaD12*sigmahyp23 -
			     omegaD13*sigmahyp11 +
			     (a1*fs_*(scalar(6.0)*a1*srnzT13 + sqrt(scalar(2.0))*fd_*sigmahypndev13*srnz))/scalar(6.0) +
			     (fs_*sigmahyp13*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 +
					      a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2));

  volScalarField Dsigmahyp21("Dsigmahyp21", -phidivsigmahyp21 + omegaD21*(sigmahyp11 - sigmahyp22) +
			     (a1*fs_*(scalar(6.0)*a1*srnzT21 + sqrt(scalar(2.0))*fd_*sigmahypndev21*srnz))/scalar(6.0) +
			     sigmahyp21*(-omegaD11 + omegaD22 +
					 (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 +
					       sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

 volScalarField Dsigmahyp22("Dsigmahyp22", -phidivsigmahyp22 + omegaD21*sigmahyp12 - omegaD12*sigmahyp21
			    + (a1*fs_*(scalar(6.0)*a1*srnzT22 + sqrt(scalar(2.0))*fd_*sigmahypndev22*srnz))/scalar(6.0)
			    + (fs_*sigmahyp22* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))
			    *pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp22*srnz)*trsigmahypinv);

 volScalarField Dsigmahyp23("Dsigmahyp23", -phidivsigmahyp23 + omegaD21*sigmahyp13 - omegaD13*sigmahyp21 +
			    omegaD22*sigmahyp23 - omegaD23*sigmahyp22 +
			    (a1*fs_*(scalar(6.0)*a1*srnzT23 + sqrt(scalar(2.0))*fd_*sigmahypndev23*srnz))/scalar(6.0) +
			    (fs_*sigmahyp23*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 +
					     sigmahyp22*srnzT22))*pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp23*srnz)*trsigmahypinv);

  */

 
 Info<< min(Dsigmahyp11)<<endl;
 Info<< min(Dsigmahyp12)<<endl;
 Info<< min(Dsigmahyp13)<<endl;
 Info<< min(Dsigmahyp21)<<endl;
 Info<< min(Dsigmahyp22)<<endl;
 Info<< min(Dsigmahyp23)<<endl;

 Info<< max(Dsigmahyp11)<<endl;
 Info<< max(Dsigmahyp12)<<endl;
 Info<< max(Dsigmahyp13)<<endl;
 Info<< max(Dsigmahyp21)<<endl;
 Info<< max(Dsigmahyp22)<<endl;
 Info<< max(Dsigmahyp23)<<endl;

 

 dimensionedScalar dt = U_.mesh().time().deltaT();
 sigmaHyp.replace(0, sigmaHyp.component(tensor::XX)+Dsigmahyp11*dt);
 sigmaHyp.replace(1, sigmaHyp.component(tensor::XY)+Dsigmahyp12*dt);
 sigmaHyp.replace(2, sigmaHyp.component(tensor::XZ)+Dsigmahyp13*dt);
 sigmaHyp.replace(3, sigmaHyp.component(tensor::YX)+Dsigmahyp21*dt);
 sigmaHyp.replace(4, sigmaHyp.component(tensor::YY)+Dsigmahyp22*dt);
 sigmaHyp.replace(5, sigmaHyp.component(tensor::YZ)+Dsigmahyp23*dt);

 Info<<"sigmaHyp constructed!!"<<endl;

 
  bool debug = true;
  if (debug)
    {
      Info<<min(sigmaHyp.component(tensor::XX))<<endl;
      Info<<max(sigmaHyp.component(tensor::XX))<<endl;
      Info<<min(sigmaHyp.component(tensor::XY))<<endl;
      Info<<max(sigmaHyp.component(tensor::XY))<<endl;
    }

  Info<<"sigmaHyp displayed!!"<<endl;
}


// ************************************************************************* //
