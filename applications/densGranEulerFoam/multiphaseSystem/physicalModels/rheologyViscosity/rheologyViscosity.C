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

#include "rheologyViscosity.H"
#include "addToRunTimeSelectionTable.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace physicalModels
{
    defineTypeNameAndDebug(rheology, 0);

    addToRunTimeSelectionTable
    (
        physicalModel,
        rheology,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::physicalModels::rheology::rheology
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
    J0_("J0",dimensionSet(0, 0, 0, 0, 0), dict),
    phiM_("phiM",dimensionSet(0, 0, 0, 0, 0), dict),
    S0_("S0",dimensionSet(0, 0, -1, 0, 0), dict),
    aMJ_("aMJ",dimensionSet(1, -1, -2, 0, 0), dict),
    phiRLP_("phiRLP",dimensionSet(0, 0, 0, 0, 0), dict),
    phiRCP_("phiRCP",dimensionSet(0, 0, 0, 0, 0), dict),
    alphaDeltaMin_("alphaDeltaMin",dimensionSet(0, 0, 0, 0, 0), dict),
    numax_("numax",dimensionSet(0, 2, -1, 0, 0), dict),
    numin_("numin",dimensionSet(0, 2, -1, 0, 0), dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::physicalModels::rheology::~rheology()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::physicalModels::rheology::calcNu() const
{

      // Local references
  const volScalarField& alpha_ = phase_;
  const volScalarField& da = phase_.d();
  const dimensionedScalar& rhoS_ = phase_.rho();
  const volVectorField& U = phase_.U();    
  tmp<volTensorField> tgradU(fvc::grad(U));
  const volTensorField& gradU(tgradU());
  volSymmTensorField D(symm(gradU));
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, scalar(1e-6));
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
  
  
  
  volScalarField srnz("srnz", min(mag(D)-tr(gradU)/scalar(3),S0_));
  
  tmp<volScalarField> phiPM
    (
     new volScalarField
     (
      IOobject
      (
       "phiPM",
       U.mesh().time().timeName(),
       U.mesh(),
       IOobject::NO_READ,
       IOobject::NO_WRITE,
       false
       ),
      U.mesh(),
      dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
      )
     );
  // volScalarField& phiPM_ = tphiPM.ref();
  forAll(srnz,ii)
    {
      if(S0_.value() > srnz[ii])
	{
	  phiPM.ref()[ii] = phiM_.value()+(phiRCP_.value()-phiM_.value())*(S0_.value()-srnz[ii]);
	}
      else
	{
	  phiPM.ref()[ii] = phiM_.value();
	}
    }
  
  const fvPatchList& patches = U.mesh().boundary();
  forAll(patches, patchID)
    {
      forAll (srnz.boundaryFieldRef()[patchID],facei)
	{
	  if(S0_.value() > srnz.boundaryFieldRef()[patchID][facei])
	    {
	      phiPM.ref().boundaryFieldRef()[patchID][facei] =
		phiM_.value()+(phiRCP_.value()-phiM_.value())*
		(S0_.value()-srnz.boundaryField()[patchID][facei]);
	    }
	  else
	    {
	      phiPM.ref().boundaryFieldRef()[patchID][facei] = phiM_.value();
	    }
	}
    }
  // Correct coupled BCs
  phiPM.ref().correctBoundaryConditions();    
  /*  volScalarField pS_ ("pS_", max(
				 aMJ_*max(alpha_-phiRLP_,scalar(0))/
				 max((phiRCP_-alpha_),alphaDeltaMin_)+
				 scalar(2)*nuC_*rhoC_*srnz/
				 pow(phiPM/max(alpha_-scalar(1),alphaDeltaMin_),2),smallval_p));
  */
  volScalarField pS_ ("pS_", max(U.mesh().lookupObject<volScalarField>("p"), smallval_p));
  volScalarField J1_ = scalar(2)*srnz*nuC_*rhoC_/max(pS_,smallval_p);
  volScalarField muJ_ = mu1_+(mu2_-mu1_)*J1_/max(J1_+J0_,scalar(1e-6))+
    J1_+scalar(5)/scalar(2)*phiM_*sqrt(J1_);

  //------------------------------------------------------------------------------//
  /*
  //    volScalarField& nu0_ = tnu0.ref();
    Info<<"cf.4\n";
    forAll(srnz,ii)
    {
      nu0_.ref()[ii] = muJ_[ii]*pS_[ii]/2/rhoS_.value()/max(srnz[ii],smallval_sr.value());
    }

  Info<<"cf.5\n";
  forAll(patches, patchID)
    {
      forAll (srnz.boundaryFieldRef()[patchID],facei)
        {
	  nu0_.ref().boundaryFieldRef()[patchID][facei] = muJ_.boundaryFieldRef()[patchID][facei]
	    *pS_.boundaryFieldRef()[patchID][facei]
	    /2/rhoS_.value()
	    /max(srnz.boundaryFieldRef()[patchID][facei],smallval_sr.value());
        }
    }
  

  forAll(nu0_,ii)
    {
      nu0_[ii] = muJ_[ii]*pS_[ii]/2/rhoS_.value()/max(srnz[ii],smallval_sr.value());
    }  
  */
  //  nu0_ = muJ_*pS_/2/rhoS_/max(srnz,smallval_sr);
  
  volScalarField pS_mjpj
    ("pS_mjpj", pS_);
  
  if(U.time().outputTime())
    {
      pS_mjpj.write();
      srnz.write();
    } 
  
  
  return min(max(muJ_*pS_/2/rhoS_/max(srnz,smallval_sr),numin_),numax_);
  /*tmp<Foam::volScalarField>
    (
     new volScalarField
     (
      IOobject
      (
       "nu",
       phase_.U().time().timeName(),
       phase_.U().mesh()
       ),
      phase_.U().mesh(),
      nu0_*scalar(1.0)
      )
      );*/
}
Foam::tmp<Foam::volScalarField> Foam::physicalModels::rheology::calcGranP() const
{
  const volScalarField& alpha_ = phase_;
  const volScalarField& da = phase_.d();
  const dimensionedScalar& rhoS_ = phase_.rho();
  const volVectorField& U = phase_.U();
  tmp<volTensorField> tgradU(fvc::grad(U));
  const volTensorField& gradU(tgradU());
  volSymmTensorField D(symm(gradU));
  dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, scalar(1e-6));
  dimensionedScalar smallval_p("smallval_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
  dimensionedScalar pMax("pMax", dimMass/dimLength/dimTime/dimTime, scalar(1e3));


  volScalarField srnz("srnz", min(mag(D)-tr(gradU)/scalar(3),S0_));

    tmp<volScalarField> phiPM
    (
     new volScalarField
     (
      IOobject
      (
       "phiPM",
       U.mesh().time().timeName(),
       U.mesh(),
       IOobject::NO_READ,
       IOobject::NO_WRITE,
       false
       ),
      U.mesh(),
      dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
      )
     );
  // volScalarField& phiPM_ = tphiPM.ref();
  forAll(srnz,ii)
    {
      if(S0_.value() > srnz[ii])
        {
          phiPM.ref()[ii] = phiM_.value()+(phiRCP_.value()-phiM_.value())*(S0_.value()-srnz[ii]);
        }
      else
        {
          phiPM.ref()[ii] = phiM_.value();
        }
    }

  const fvPatchList& patches = U.mesh().boundary();
  forAll(patches, patchID)
    {
      forAll (srnz.boundaryFieldRef()[patchID],facei)
        {
          if(S0_.value() > srnz.boundaryFieldRef()[patchID][facei])
            {
              phiPM.ref().boundaryFieldRef()[patchID][facei] =
                phiM_.value()+(phiRCP_.value()-phiM_.value())*
                (S0_.value()-srnz.boundaryField()[patchID][facei]);
            }
          else
            {
              phiPM.ref().boundaryFieldRef()[patchID][facei] = phiM_.value();
            }
        }
    }
  // Correct coupled BCs
  phiPM.ref().correctBoundaryConditions();

  volScalarField pS_mjpj_2
    ("pS_mjpj_2",min(max(aMJ_*max(alpha_-phiRLP_,scalar(0))/
                 max((phiRCP_-alpha_),alphaDeltaMin_)+
                 scalar(2)*nuC_*rhoC_*srnz/
                 pow(phiPM/max(alpha_-scalar(1),alphaDeltaMin_),2),smallval_p)
             , pMax));
  
  return min(max(aMJ_*max(alpha_-phiRLP_,scalar(0))/
		 max((phiRCP_-alpha_),alphaDeltaMin_)+
		 scalar(2)*nuC_*rhoC_*srnz/
		 pow(phiPM/max(alpha_-scalar(1),alphaDeltaMin_),2),smallval_p)
	     , pMax);
}

// ************************************************************************* //
