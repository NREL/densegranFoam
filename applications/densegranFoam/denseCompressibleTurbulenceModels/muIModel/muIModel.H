/*---------------------------------------------------------------------------* \
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

Class
    Foam::muIModel

Description
    Particle-particle phase-pressure RAS model

    The derivative of the phase-pressure with respect to the phase-fraction
    is evaluated as

        g0*min(exp(preAlphaExp*(alpha - alphaMax)), expMax)

    The default model coefficients correspond to the following:
    \verbatim
        phasePressureCoeffs
        {
            preAlphaExp     500;
            expMax          1000;
            alphaMax        0.62;
            g0              1000;
        }
    \endverbatim

SourceFiles
    muIModel.C

\*---------------------------------------------------------------------------*/

#ifndef muIModel_H
#define muIModel_H


#include "RASModel.H"
#include "eddyViscosity.H"
#include "denseCompressibleTurbulenceModel.H"
#include "EddyDiffusivity.H"
#include "phaseModel.H"
#include "dragModel.H"
#include "viscosityModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class muIModel Declaration
\*---------------------------------------------------------------------------*/

class muIModel
:
    public eddyViscosity
    <
        RASModel<EddyDiffusivity<denseCompressibleTurbulenceModel>>
    >
{
    // Private data

        // Input Fields

  const phaseModel& phase_;
  

  // Sub-models
  
  
  
  dimensionedScalar Hmax_;
  dimensionedScalar L_;
  dimensionedScalar nu0_;
  dimensionedScalar muS_;
  dimensionedScalar mu2_;
  scalar nNormalCell;
  scalar alpha0_;
  scalar phiMuI_;
  scalar Cs_;
  scalar ampl_;
  scalar beta_;
  scalar theta_;
  dimensionedScalar rhoS;
  scalar relaxPa_;
  scalar nuMin_;
  dimensionedScalar Fr_;
  dimensionedScalar etaPF_;
  dimensionedScalar p_;
  dimensionedScalar alphaDeltaMin_;
  dimensionedScalar alphaMinFriction;
  dimensionedScalar Bphi;

  dimensionedScalar g0_;
  volScalarField ppres;
  uniformDimensionedVectorField g_;
  volTensorField tau;
  volScalarField eta;
  //- The granular bulk viscosity
  volScalarField lambda_;
  volScalarField href_nlgf;
  
  // Private Member Functions
  
  void correctNut()
  {}
  
  //- Disallow default bitwise copy construct
  muIModel(const muIModel&);

  //- Disallow default bitwise assignment
  void operator=(const muIModel&);
  
  
  // Sub-models
  
  
  
  
  
public:
  
  //- Runtime type information
  TypeName("muI");
  
  
  // Constructors
  
  //- Construct from components
  muIModel
  (
   const volScalarField& alpha,
   const volScalarField& rho,
   const volVectorField& U,
   const surfaceScalarField& alphaRhoPhi,
   const surfaceScalarField& phi,
   const phaseModel& transport,
   const word& propertiesName = turbulenceModel::propertiesName,
   const word& type = typeName
   );
  
  
  //- Destructor
  virtual ~muIModel();
  
  
  // Member Functions
  
  //- Re-read model coefficients if they have changed
  virtual bool read();
  
  //- Return the effective viscosity
  virtual tmp<volScalarField> nuEff() const
  {
    return this->nut();
  }
  
  //- Return the effective viscosity on patch
  virtual tmp<scalarField> nuEff(const label patchi) const
  {
    return this->nut(patchi);
  }
  
  //- Return the turbulence kinetic energy
  virtual tmp<volScalarField> k() const;
  
  //- Return the turbulence kinetic energy dissipation rate
  virtual tmp<volScalarField> epsilon() const;
  
  //- Return the Reynolds stress tensor
  virtual tmp<volSymmTensorField> R() const;
  
  //- Return the phase-pressure'
  // (derivative of phase-pressure w.r.t. phase-fraction)
  virtual tmp<volScalarField> pPrime() const;
  
  //- Return the face-phase-pressure'
  // (derivative of phase-pressure w.r.t. phase-fraction)
  virtual tmp<surfaceScalarField> pPrimef() const;
  
  //- Return the effective stress tensor
  virtual tmp<volSymmTensorField> devRhoReff() const;
  
  //- Return the source term for the momentum equation
  virtual tmp<fvVectorMatrix> divDevRhoReff(volVectorField& U) const;

  virtual tmp<volScalarField>nlgfPre() const;
  
  //- Solve the kinetic theory equations and correct the viscosity
  virtual void correct();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
