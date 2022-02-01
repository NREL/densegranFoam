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
  sigmaHyp
  (
   IOobject
   (
    "sigmaHyp",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
   rho_*nut_*2*dev(symm(fvc::grad(U)))*scalar(0.0)
   )
{
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
      dev(sigmaHyp)
      )
     );
}

Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::hypoplasticModel::RK4solver2D(volSymmTensorField sigmahyp) const
{
  dimensionedScalar rk_h = U_.mesh().time().deltaT();

  volScalarField sh11("sh11", sigmahyp.component(tensor::XX));
  volScalarField sh12("sh12", sigmahyp.component(tensor::XY));
  volScalarField sh13("sh13", sigmahyp.component(tensor::XZ));
  volScalarField sh21("sh21", sigmahyp.component(tensor::YX));
  volScalarField sh22("sh22", sigmahyp.component(tensor::YY));
  volScalarField sh23("sh23", sigmahyp.component(tensor::YZ));

  sigmahyp = this->sigmahypfunc(sh11,sh12,sh13,sh21,sh22,sh23);
  volScalarField k1_11("k1_11", sigmahyp.component(tensor::XX));
  volScalarField k1_12("k1_12", sigmahyp.component(tensor::XY));
  volScalarField k1_13("k1_13", sigmahyp.component(tensor::XZ));
  volScalarField k1_21("k1_21", sigmahyp.component(tensor::YX));
  volScalarField k1_22("k1_22", sigmahyp.component(tensor::YY));
  volScalarField k1_23("k1_23", sigmahyp.component(tensor::YZ));
  sigmahyp = this->sigmahypfunc(sh11+rk_h*k1_11/2,sh12+rk_h*k1_12/2,sh13+rk_h*k1_13/2,sh21+rk_h*k1_21/2,sh22+rk_h*k1_22/2,sh23+rk_h*k1_23/2);
  volScalarField k2_11("k2_11", sigmahyp.component(tensor::XX));
  volScalarField k2_12("k2_12", sigmahyp.component(tensor::XY));
  volScalarField k2_13("k2_13", sigmahyp.component(tensor::XZ));
  volScalarField k2_21("k2_21", sigmahyp.component(tensor::YX));
  volScalarField k2_22("k2_22", sigmahyp.component(tensor::YY));
  volScalarField k2_23("k2_23", sigmahyp.component(tensor::YZ));
  sigmahyp = this->sigmahypfunc(sh11+rk_h*k2_11/2,sh12+rk_h*k2_12/2,sh13+rk_h*k2_13/2,sh21+rk_h*k2_21/2,sh22+rk_h*k2_22/2,sh23+rk_h*k2_23/2);
  volScalarField k3_11("k3_11", sigmahyp.component(tensor::XX));
  volScalarField k3_12("k3_12", sigmahyp.component(tensor::XY));
  volScalarField k3_13("k3_13", sigmahyp.component(tensor::XZ));
  volScalarField k3_21("k3_21", sigmahyp.component(tensor::YX));
  volScalarField k3_22("k3_22", sigmahyp.component(tensor::YY));
  volScalarField k3_23("k3_23", sigmahyp.component(tensor::YZ));
  sigmahyp = this->sigmahypfunc(sh11+rk_h*k3_11,sh12+rk_h*k3_12,sh13+rk_h*k3_13,sh21+rk_h*k3_21,sh22+rk_h*k3_22,sh23+rk_h*k3_23);
  volScalarField k4_11("k4_11", sigmahyp.component(tensor::XX));
  volScalarField k4_12("k4_12", sigmahyp.component(tensor::XY));
  volScalarField k4_13("k4_13", sigmahyp.component(tensor::XZ));
  volScalarField k4_21("k4_21", sigmahyp.component(tensor::YX));
  volScalarField k4_22("k4_22", sigmahyp.component(tensor::YY));
  volScalarField k4_23("k4_23", sigmahyp.component(tensor::YZ));
  sh11 = sh11 + 1/6*rk_h*(k1_11+2*k2_11+2*k3_11+k4_11);
  sh12 = sh12 + 1/6*rk_h*(k1_12+2*k2_12+2*k3_12+k4_12);
  sh13 = sh13 + 1/6*rk_h*(k1_13+2*k2_13+2*k3_13+k4_13);
  sh21 = sh21 + 1/6*rk_h*(k1_21+2*k2_21+2*k3_21+k4_21);
  sh22 = sh22 + 1/6*rk_h*(k1_22+2*k2_22+2*k3_22+k4_22);
  sh23 = sh23 + 1/6*rk_h*(k1_23+2*k2_23+2*k3_23+k4_23);

  sigmahyp.replace(0, sh11);
  sigmahyp.replace(1, sh12);
  sigmahyp.replace(2, sh13);
  sigmahyp.replace(3, sh21);
  sigmahyp.replace(4, sh22);
  sigmahyp.replace(5, sh23);
  
  return sigmahyp;
}

Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::hypoplasticModel::sigmahypfunc(volScalarField sigmahyp11, volScalarField sigmahyp12, volScalarField sigmahyp13,
						  volScalarField sigmahyp21, volScalarField sigmahyp22, volScalarField sigmahyp23) const
{
  volSymmTensorField sigmahyp("sigmahyp", rho_*nut_*2*dev(symm(fvc::grad(U_)))*scalar(0.0));
  sigmahyp.replace(0, sigmahyp11);
  sigmahyp.replace(1, sigmahyp12);
  sigmahyp.replace(2, sigmahyp13);
  sigmahyp.replace(3, sigmahyp21);
  sigmahyp.replace(4, sigmahyp22);
  sigmahyp.replace(5, sigmahyp23);
  volSymmTensorField phidivsigmahyp("phidivsigmahyp", fvc::div(phi_, sigmahyp));
  volScalarField phidivsigmahyp11("phidivsigmahyp11", phidivsigmahyp.component(tensor::XX));
  volScalarField phidivsigmahyp12("phidivsigmahyp12", phidivsigmahyp.component(tensor::XY));
  volScalarField phidivsigmahyp13("phidivsigmahyp13", phidivsigmahyp.component(tensor::XZ));
  volScalarField phidivsigmahyp21("phidivsigmahyp21", phidivsigmahyp.component(tensor::YX));
  volScalarField phidivsigmahyp22("phidivsigmahyp22", phidivsigmahyp.component(tensor::YY));
  volScalarField phidivsigmahyp23("phidivsigmahyp23", phidivsigmahyp.component(tensor::YZ));
  
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-6);
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, 1e-6);
  dimensionedScalar limhs("limhs", dimMass/dimLength/dimTime/dimTime, 1.0);
  scalar smallVal = scalar(1e-6);
  volScalarField srnz("srnz", sqrt(2.0)*mag(symm(fvc::grad(U_))));
  volSymmTensorField srnzT("srnzT", symm(fvc::grad(U_)));
  volScalarField srnzT11("srnzT11", srnzT.component(tensor::XX));
  volScalarField srnzT12("srnzT12", srnzT.component(tensor::XY));
  volScalarField srnzT13("srnzT13", srnzT.component(tensor::XZ));
  volScalarField srnzT21("srnzT21", srnzT.component(tensor::YX));
  volScalarField srnzT22("srnzT22", srnzT.component(tensor::YY));
  volScalarField srnzT23("srnzT23", srnzT.component(tensor::YZ));
  volTensorField omegaD("omegaD", (0.5)*(fvc::grad(U_)-symm(fvc::grad(U_))));
  volScalarField omegaD11("omegaD11", omegaD.component(tensor::XX));
  volScalarField omegaD12("omegaD12", omegaD.component(tensor::XY));
  volScalarField omegaD13("omegaD13", omegaD.component(tensor::XZ));
  volScalarField omegaD21("omegaD21", omegaD.component(tensor::YX));
  volScalarField omegaD22("omegaD22", omegaD.component(tensor::YY));
  volScalarField omegaD23("omegaD23", omegaD.component(tensor::YZ));
  volScalarField trsigmahyp("trsigmahyp", tr(sigmahyp));
  volScalarField trsigmahypinv("trsigmahypinv", 1/max(smallval_p,tr(sigmahyp)));
  volSymmTensorField sigmahypn("sigmahypn", sigmahyp*trsigmahypinv);
  volSymmTensorField sigmahypndev("sigmahypndev", dev(sigmahypn));
  volScalarField sigmahypndev11("sigmahypndev11", sigmahypndev.component(tensor::XX));
  volScalarField sigmahypndev12("sigmahypndev12", sigmahypndev.component(tensor::XY));
  volScalarField sigmahypndev13("sigmahypndev13", sigmahypndev.component(tensor::XZ));
  volScalarField sigmahypndev21("sigmahypndev21", sigmahypndev.component(tensor::YX));
  volScalarField sigmahypndev22("sigmahypndev22", sigmahypndev.component(tensor::YY));
  volScalarField sigmahypndev23("sigmahypndev23", sigmahypndev.component(tensor::YZ));
  scalar Theta_ = ThetaC_*M_PI/scalar(180);
  scalar c1_ = sqrt(scalar(3.0)/scalar(8.0))*(scalar(3)-sin(Theta_))/max(sin(Theta_),smallVal);
  scalar c2_ = scalar(3.0)/scalar(8.0)*(3+sin(Theta_))/max(sin(Theta_),smallVal);
  volScalarField ei ("ei", ei0_*exp(-pow(-tr(sigmahyp)/max(hs_,smallval_p),1-n_)));
  volScalarField ed ("ed", ed0_*exp(-pow(-tr(sigmahyp)/max(hs_,smallval_p),1-n_)));
  volScalarField ec ("ec", ec0_*exp(-pow(-tr(sigmahyp)/max(hs_,smallval_p),1-n_)));
  volScalarField voidRatio("voidRatio", max(min(scalar(1.0)/max(alpha_,smallVal)-scalar(1.0),ec),ed));
  
  scalar hi = 1/max(c1_*c1_,smallVal)+scalar(1)/scalar(3)*pow((ei0_+ed0_)/max((ec0_-ed0_),smallVal),alphaHyp_)
    /max(c1_,smallVal)/sqrt(scalar(3));
  volScalarField fs_
    ("fs_",
     hs_/n_/max(hi,smallVal)*pow(ei/max(voidRatio,smallVal),beta_)*(1+ei)/max(ei,smallVal)*pow(-tr(sigmahyp)/max(hs_,smallval_p),1-n_)
     );
  volScalarField fd_("fd_", pow((voidRatio-ed)/(ec-ed),alphaHyp_));
  volScalarField a1("a1", 1/max((c1_ + c2_*mag(sigmahypndev)*(1+cos(M_PI-3*Theta_))),smallVal));

  Info<<"constracting Dsigmahyp"<<endl;

  volScalarField Dsigmahyp11("Dsigmahyp11", -phidivsigmahyp11 - omegaD21*sigmahyp12 + omegaD12*sigmahyp21 +
			     (a1*fs_*(scalar(6.0)*a1*srnzT11 + sqrt(scalar(2.0))*fd_*sigmahypndev11*srnz))/scalar(6.0) +
			     (fs_*sigmahyp11*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))*
			     pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp11*srnz)*trsigmahypinv);
  
 volScalarField Dsigmahyp12("Dsigmahyp12",-phidivsigmahyp12 + omegaD12*(-sigmahyp11 + sigmahyp22) + (a1*fs_*(scalar(6.0)*a1*srnzT12 + sqrt(scalar(2.0))*fd_*sigmahypndev12*srnz))/scalar(6.0) + sigmahyp12*(omegaD11 - omegaD22 + (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

 volScalarField Dsigmahyp13("Dsigmahyp13",-phidivsigmahyp13 - omegaD23*sigmahyp12 + omegaD11*sigmahyp13 + omegaD12*sigmahyp23 - omegaD13*sigmahyp11 + (a1*fs_*(scalar(6.0)*a1*srnzT13 + sqrt(scalar(2.0))*fd_*sigmahypndev13*srnz))/scalar(6.0) + (fs_*sigmahyp13*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2));

 volScalarField Dsigmahyp21("Dsigmahyp21",-phidivsigmahyp21 + omegaD21*(sigmahyp11 - sigmahyp22) + (a1*fs_*(scalar(6.0)*a1*srnzT21 + sqrt(scalar(2.0))*fd_*sigmahypndev21*srnz))/scalar(6.0) + sigmahyp21*(-omegaD11 + omegaD22 + (fs_*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22 + a1*fd_*srnz*trsigmahyp))*pow(trsigmahypinv,2)));

 volScalarField Dsigmahyp22("Dsigmahyp22",-phidivsigmahyp22 + omegaD21*sigmahyp12 - omegaD12*sigmahyp21
			    + (a1*fs_*(scalar(6.0)*a1*srnzT22 + sqrt(scalar(2.0))*fd_*sigmahypndev22*srnz))/scalar(6.0)
			    + (fs_*sigmahyp22* (sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))
			    *pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp22*srnz)*trsigmahypinv);

 volScalarField Dsigmahyp23("Dsigmahyp23",-phidivsigmahyp23 + omegaD21*sigmahyp13 - omegaD13*sigmahyp21 + omegaD22*sigmahyp23 - omegaD23*sigmahyp22 + (a1*fs_*(scalar(6.0)*a1*srnzT23 + sqrt(scalar(2.0))*fd_*sigmahypndev23*srnz))/scalar(6.0) + (fs_*sigmahyp23*(sigmahyp11*srnzT11 + sigmahyp21*srnzT12 + sigmahyp12*srnzT21 + sigmahyp22*srnzT22))*pow(trsigmahypinv,2) + (a1*fd_*fs_*sigmahyp23*srnz)*trsigmahypinv);

 sigmahyp.replace(0, Dsigmahyp11);
 sigmahyp.replace(1, Dsigmahyp12);
 sigmahyp.replace(2, Dsigmahyp13);
 sigmahyp.replace(3, Dsigmahyp21);
 sigmahyp.replace(4, Dsigmahyp22);
 sigmahyp.replace(5, Dsigmahyp23);

 Info<<"Returning sigmaHyp"<<endl;
  return sigmahyp;
    
}

void Foam::RASModels::hypoplasticModel::correct()
{
  Info<<"Constructing sigmaHyp:"<<endl;
  sigmaHyp == this->RK4solver2D(sigmaHyp);
  Info<<"Constructed sigmaHyp"<<endl;


  bool debug = true;
  if (debug)
    {
      Info<<min(sigmaHyp.component(tensor::XX))<<endl;
      Info<<max(sigmaHyp.component(tensor::XX))<<endl;
    }
}


// ************************************************************************* //
