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

#include "biomassModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"
#include "fvOptions.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::biomassModel::biomassModel
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

  viscosityModel_
    (
        biomassModels::viscosityModel::New
        (
         type,
         coeffDict_,
         U,
         phi,
     	 alpha,
         rho,
         alphaRhoPhi
        )
    ),
  fcicGranularPressureModel_
  (
   biomassModels::fcicGranularPressureModel::New
        (
         coeffDict_,
         U,
         alpha,
         rho
        )
        ),
    frictionalStressModel_
    (
        biomassModels::frictionalStressModel::New
        (
            coeffDict_
        )
    ),

  alphaMax_(readScalar(coeffDict_.lookup("alphaMax"))),
  alphaMinFriction_(readScalar(coeffDict_.lookup("alphaMinFriction"))),
    preAlphaExp_(readScalar(coeffDict_.lookup("preAlphaExp"))),
    expMax_(readScalar(coeffDict_.lookup("expMax"))),
    g0_
    (
        "g0",
        dimensionSet(1, -1, -2, 0, 0),
        coeffDict_.lookup("g0")
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
    )
{
 
    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::biomassModel::~biomassModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::biomassModel::read()
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
	coeffDict().lookup("alphaMinFriction") >> alphaMinFriction_;
        coeffDict().lookup("preAlphaExp") >> preAlphaExp_;
        coeffDict().lookup("expMax") >> expMax_;
        g0_.readIfPresent(coeffDict());
	fcicGranularPressureModel_->read();
        viscosityModel_->read();
	fcicGranularPressureModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::biomassModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::biomassModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::biomassModel::R() const
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
Foam::RASModels::biomassModel::pPrime() const
{
  tmp<volScalarField> tpPrime
    (
         frictionalStressModel_->frictionalPressurePrime
         (
          phase_,
          alphaMinFriction_,
          alphaMax_
          )
	  /*     +
          fcicGranularPressureModel_->granularPressureCoeff()
	  */	  
 +
	 g0_
       *min
	(
            exp(preAlphaExp_*(alpha_ - alphaMax_)),
            expMax_
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

    volScalarField tpwr
      (
	"tpwr",
	  g0_
	  *min
	  (
	   exp(preAlphaExp_*(alpha_ - alphaMax_)),
	   expMax_
	   )
       );
    
    if(U_.time().outputTime())
    {
      tpwr.write();
    }    

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::biomassModel::pPrimef() const
{
  return fvc::interpolate(pPrime());
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::biomassModel::devRhoReff() const
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
           *dev(twoSymm(fvc::grad(U_)))
          - ((rho_*lambda_)*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::biomassModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            (rho_*nut_)*dev2(T(fvc::grad(U)))
          + ((rho_*lambda_)*fvc::div(phi_))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
}


void Foam::RASModels::biomassModel::correct()
{
    volScalarField pf
        (
         "pf",
         frictionalStressModel_->frictionalPressure
         (
          phase_,
          alphaMinFriction_,
          alphaMax_
          )
         );
    volScalarField pa
      (
       "pa",
       fcicGranularPressureModel_->granularPressureCoeff()
       );

    nut_ = viscosityModel_->calcNu(pf,pa);
    lambda_ = scalar(0.0)*nut_;

}


// ************************************************************************* //
