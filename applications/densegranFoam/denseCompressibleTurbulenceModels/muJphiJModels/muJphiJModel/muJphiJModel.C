
/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "muJphiJModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::muJphiJModel::muJphiJModel
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
        RASModel<EddyDiffusivity<denseCompressibleTurbulenceModel>>
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

    viscosityModel_
    (
        muJphiJModels::viscosityModel::New
        (
	 coeffDict_,
	 U,
	 phi,
	 alpha
        )
    ),
    conductivityModel_
    (
        muJphiJModels::conductivityModel::New
        (
            coeffDict_
        )
    ),
    radialModel_
    (
        muJphiJModels::radialModel::New
        (
            coeffDict_
        )
    ),
    effectivePressureModel_
    (
        muJphiJModels::effectivePressureModel::New
        (
            coeffDict_
        )
    ),
    granularPressureModel_
    (
        muJphiJModels::granularPressureModel::New
        (
            coeffDict_
        )
    ),
    frictionalStressModel_
    (
        muJphiJModels::frictionalStressModel::New
        (
            coeffDict_
        )
    ),

    e_("e", dimless, coeffDict_),
    alphaMax_("alphaMax", dimless, coeffDict_),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict_
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffDict_
    ),
    rhoC_("rhoC", dimensionSet(1,-3,0,0,0), coeffDict_),
    nuC_("nuC", dimensionSet(0,2,-1,0,0), coeffDict_),
    mu1_("mu1", dimensionSet(0,0,0,0,0), coeffDict_),
    mu2_("mu2", dimensionSet(0,0,0,0,0), coeffDict_),
    J0_("J0", dimensionSet(0,0,0,0,0), coeffDict_),
    rhoS_("rhoS", dimensionSet(1,-3,0,0,0), coeffDict_),
    phiM_("phiM", dimensionSet(0,0,0,0,0), coeffDict_),
    S0_("S0", dimensionSet(0,0,-1,0,0), coeffDict_),
    A_("A", dimless, coeffDict_),
    t0inv_("t0inv", dimless/dimTime, coeffDict_),
    relaxPa_("relaxPa", dimensionSet(0, 0, 0, 0, 0, 0, 0), coeffDict_),
    g_nlgf
    (
     IOobject
     (
      "g_nlgf",
      U_.time().timeName(),
      U_.mesh(),
      IOobject::READ_IF_PRESENT,
      IOobject::AUTO_WRITE
      ),
     U_.mesh()
     ),

    pG_
    (
     IOobject
     (
      "pG",
      U_.time().timeName(),
      U_.mesh(),
      IOobject::MUST_READ,
      IOobject::AUTO_WRITE
      ),
     U_.mesh()
     ),
    
    maxNut_
    (
        "maxNut",
        dimensionSet(0,2,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("maxNut",1000)
    ),

    minNut_
    (
        "minNut",
        dimensionSet(0,2,-1,0,0),
        coeffDict_.lookupOrDefault<scalar>("minNut",1e-6)
    ),

    p(U_.mesh().lookupObject<volScalarField>("p")),
    
    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    
    nuFric_
    (
        IOobject
        (
            IOobject::groupName("nuFric", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    )
{
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::muJphiJModel::~muJphiJModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::muJphiJModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<denseCompressibleTurbulenceModel>>
        >::read()
    )
    {
        e_.readIfPresent(coeffDict());
        alphaMax_.readIfPresent(coeffDict());
        alphaMinFriction_.readIfPresent(coeffDict());

	coeffDict().lookup("rhoC") >> rhoC_;
	coeffDict().lookup("nuC") >> nuC_;
	coeffDict().lookup("mu1") >> mu1_;
	coeffDict().lookup("mu2") >> mu2_;
	coeffDict().lookup("J0") >> J0_;
	coeffDict().lookup("rhoS") >> rhoS_;
        coeffDict().lookup("phiM") >> phiM_;
	coeffDict().lookup("S0") >> S0_;
	coeffDict().lookup("relaxPa") >> relaxPa_;
	A_.readIfPresent(coeffDict());
	t0inv_.readIfPresent(coeffDict());

        viscosityModel_->read();
        conductivityModel_->read();
        radialModel_->read();
	effectivePressureModel_->read();
        granularPressureModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::muJphiJModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::muJphiJModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::muJphiJModel::R() const
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
          - (nut_)*dev(twoSymm(fvc::grad(U_)))
          - (lambda_*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::muJphiJModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();

    tmp<volScalarField> tpPrime
    (
     effectivePressureModel_->effectivePressurePrime
      (
       phase_,
       phiM_,
       nuC_,
       rhoC_,
       S0_
       )
     /*Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            alpha_,
            radialModel_->g0(alpha_, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(alpha_, alphaMinFriction_, alphaMax_),
            rho,
            e_
        )
     +  frictionalStressModel_->frictionalPressurePrime
        (
            phase_,
            alphaMinFriction_,
            alphaMax_
	    )*/
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
    
    volScalarField tpPrimeO("tpPrimeO",  
     effectivePressureModel_->effectivePressurePrime
      (
       phase_,
       phiM_,
       nuC_,
       rhoC_,
       S0_
       )
     );
    if(U_.time().outputTime())
      {tpPrimeO.write();}

    Info<<max(tpPrimeO)<<endl;
    Info<<min(tpPrimeO)<<endl;

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::muJphiJModel::pPrimef() const
{
    return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::muJphiJModel::devRhoReff() const
{
  Info<<"devRhoReff\n";
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
	    *dev(twoSymm(fvc::grad(U_)))
	    //          - ((rho_*lambda_)*fvc::div(phi_))*symmTensor::I
        )
    );
    Info<<"devRhoReff end\n";
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::muJphiJModel::divDevRhoReff
(
    volVectorField& U
) const
{
  Info<<"divDevRhoReff\n";
    return
    (
      - fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            (rho_*nut_)*dev2(T(fvc::grad(U)))
	    //          + ((rho_*lambda_)*fvc::div(phi_))
	    //           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
    Info<<"divDevRhoReff\n";
}


void Foam::RASModels::muJphiJModel::correct()
{
    // Local references
    volScalarField alpha(max(alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;
    const volVectorField& Uc_ =
        refCast<const twoPhaseSystem>(phase_.fluid()).otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));
    
    dimensionedScalar tau_inv_min("tau_inv_min", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1e-12);
    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));
    dimensionedScalar smallval_sr("smallvel_sr", dimless/dimTime, scalar(1e-6));
    dimensionedScalar smallval_p("smallvel_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
    
    Info<<"cf.1\n";
    volScalarField srnz("srnz", min(mag(D),S0_));//-tr(gradU)/scalar(3),S0_));
    Info<<"cf.2\n";
    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);
    Info<<"cf.3\n";
    // Frictional pressure
    volScalarField pf
      (
       frictionalStressModel_->frictionalPressure
       (
        phase_,
        alphaMinFriction_,
        alphaMax_
        )
       );
    //    volScalarField phiG_ = alpha_;//phiM_/(scalar(1)+sqrt(J1_));
    const fvMesh& mesh = U_.mesh();

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff
      (
       granularPressureModel_->granularPressureCoeff
       (
        alpha,
        gs0_,
        rho,
        e_
        )
       );
    
    
    volScalarField pE_ = effectivePressureModel_->effectivePressure
      (
       phase_,
       phiM_,
       nuC_,
       rhoC_,
       S0_
       );



    volScalarField tau_inv_par(relaxPa_*alpha_*srnz);
    tau_inv_par.max(tau_inv_min);

    
    fvScalarMatrix pgEqn
    (
         fvm::ddt(pG_)
         + fvm::div(phi_, pG_)
         - fvm::Sp(fvc::div(phi_), pG_)
        ==
        tau_inv_par*(pE_)
        -fvm::Sp(tau_inv_par, pG_)
    );
    pgEqn.relax();
    pgEqn.solve();

    pG_.max(smallval_p);

    pG_.correctBoundaryConditions();

    
    //    volScalarField pS_ = max(mesh.lookupObject<volScalarField>("p")-p, smallval_p);
    volScalarField pS_ = pG_;


    
    //    pS_ = max(mesh.lookupObject<volScalarField>("p")-p+pS_, smallval_p);
    
    //    volScalarField pS_ = max(aMJ_*(alpha_-phiRLP_)/max((phiRCP_-alpha_),scalar(1e-6))+scalar(2)*nuC_*rhoC_*srnz/pow(phiPM_/max(alpha_-scalar(1),scalar(1e-6)),2),smallval_p);
    volScalarField J1_("J1_",  scalar(2)*srnz*nuC_*rhoC_/max(pS_,smallval_p));
    volScalarField muJ_("muJ_",  mu1_+(mu2_-mu1_)*J1_/max(J1_+J0_,scalar(1e-6))+J1_+scalar(5)/scalar(2)*phiM_*sqrt(J1_));
    
    // Particle viscosity (Table 3.2, p.47)
    Info<<max(pS_)<<endl;
    Info<<max(srnz)<<endl;
    Info<<max(muJ_)<<endl;

    //------------------------------------------------------//

    dimensionedScalar smallval_g("smallval_g", dimless/dimTime, 1e-6);
    dimensionedScalar max_g("max_g", dimless/dimTime, 1e6);
    dimensionedScalar b_nlgf = (mu2_-mu1_)/J0_;
    volScalarField t0inv = sqrt(pS_/rho)/da;
        
    fvScalarMatrix gEqn
      (fvm::ddt(g_nlgf) + fvm::div(phi_, g_nlgf) - fvm::SuSp(fvc::div(phi_), g_nlgf)
       - fvm::laplacian(A_*A_*da*da*t0inv_, g_nlgf) ==
       fvm::Sp(-t0inv_*(mu2_-mu1_)*(mu1_ - muJ_)/(mu2_ - muJ_), g_nlgf)
       + fvm::Sp(-t0inv_*b_nlgf*muJ_* sqrt(rho*da*da/pS_)*srnz/max(muJ_,scalar(1e-6)), g_nlgf)
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
    g_nlgf.max(1e-6);
    g_nlgf.min(1e6);
    g_nlgf.correctBoundaryConditions();
    
    //------------------------------------------------------//

    nut_ = viscosityModel_->nu(muJ_, rhoS_, srnz, pS_, da, g_nlgf);
    nut_.correctBoundaryConditions();
    Info<<"cf.11\n";
    volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));
    
    // Stress tensor, Definitions, Table 3.1, p. 43
    volSymmTensorField tau
      (
       rho*(2.0*nut_*D + ( - (2.0/3.0)*nut_)*tr(D)*I)
       );

    /*
    const fvMesh& mesh = U_.mesh();
    const volScalarField& pM = mesh.lookupObject<volScalarField>("p");
    */
    volScalarField pS_mjpj
      ("pS_mjpj", pS_);
    /*    
  if(U_.time().outputTime())
    {
      pS_mjpj.write();
      srnz.write();
      J1_.write();
      muJ_.write();
      pE_.write();
      
    }
    */
    

    
    // Dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff
      (
       "gammaCoeff",
       12.0*(1.0 - sqr(e_))
       *max(sqr(alpha), residualAlpha_)
       *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
       );
    
    // Drag
    volScalarField beta
      (
       refCast<const twoPhaseSystem>(phase_.fluid()).drag(phase_).K()
       );
    
    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1("J1", 3.0*beta);
    volScalarField J2
      (
       "J2",
       0.25*sqr(beta)*da*magSqr(U - Uc_)
       /(
	 max(alpha, residualAlpha_)*rho
	 *sqrtPi*(ThetaSqrt + ThetaSmallSqrt)
	 )
       );
    
    
    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);
    
    fv::options& fvOptions(fv::options::New(mesh_));
    
    // Construct the granular temperature equation (Eq. 3.20, p. 44)
    // NB. note that there are two typos in Eq. 3.20:
    //     Ps should be without grad
    //     the laplacian has the wrong sign
    fvScalarMatrix ThetaEqn
      (
       1.5*
       (
	fvm::ddt(alpha, rho, Theta_)
	+ fvm::div(alphaRhoPhi, Theta_)
	- fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
	)
       - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
       ==
       - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
       + (tau && gradU)
       + fvm::Sp(-gammaCoeff, Theta_)
       + fvm::Sp(-J1, Theta_)
       + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
       + fvOptions(alpha, rho, Theta_)
       );
    
    ThetaEqn.relax();
    fvOptions.constrain(ThetaEqn);
    ThetaEqn.solve();
    fvOptions.correct(Theta_);
    Theta_.max(0);
    Theta_.min(100);


    
    // Limit viscosity and add frictional viscosity
    //nut_.min(maxNut_);
    // Frictional pressure
    
    nuFric_ = frictionalStressModel_->nu
      (
       phase_,
       alphaMinFriction_,
       alphaMax_,
       pf/rho,
       D
       );
    
    // Limit viscosity and add frictional viscosity
    nut_.min(maxNut_);
    //    nuFric_ = min(nuFric_, maxNut_ - nut_);
    //    nut_ += nuFric_;
    
    nut_.max(minNut_);
    
    if (debug)
    {
      Info<< typeName << ':' << nl
	  << "    max(Theta) = " << max(Theta_).value() << nl
	  << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //
