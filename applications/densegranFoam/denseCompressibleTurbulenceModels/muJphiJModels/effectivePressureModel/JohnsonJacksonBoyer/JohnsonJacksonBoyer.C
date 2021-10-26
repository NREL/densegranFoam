/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "JohnsonJacksonBoyer.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace muJphiJModels
{
namespace effectivePressureModels
{
    defineTypeNameAndDebug(JohnsonJacksonBoyer, 0);

    addToRunTimeSelectionTable
    (
        effectivePressureModel,
        JohnsonJacksonBoyer,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::muJphiJModels::effectivePressureModels::
JohnsonJacksonBoyer::JohnsonJacksonBoyer
(
    const dictionary& dict
)
:
    effectivePressureModel(dict),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    aMJ_("aMJ", dimensionSet(1, -1, -2, 0, 0), coeffDict_),
    phiRLP_("phiRLP", dimensionSet(0,0,0,0,0), coeffDict_),
    phiRCP_("phiRCP", dimensionSet(0,0,0,0,0), coeffDict_),
    alphaDeltaMin_("alphaDeltaMin", dimensionSet(0,0,0,0,0), coeffDict_)
    /*,
    phiPM_
    (
     IOobject
     (
      IOobject::groupName("phiPM", phase.name()),
      U.time().timeName(),
      U.mesh(),
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     U.mesh(),
     dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
     )*/
{/*
  forAll(srnz,ii)
    {
      if(S0_.value() > srnz[ii])
	{
	  phiPM_[ii] = phiM_.value()+(phiRCP_.value()-phiM_.value())*(S0_.value()-srnz[ii]);
	}
      else
	{
	  phiPM_[ii] = phiM_.value();
	}
    }
  
    forAll(srnz.internalField(),ii)
    {
    if(S0_.value() > srnz.internalField()[ii])
    {
    phiPM_.internalField()[ii] = phiM_+(phiRCP_-phiM_)*(S0_-srnz.internalField()[ii]);
    }
    else
    {
    phiPM_.internalField()[ii] = phiM_;
    }
    }
  
  forAll(srnz.boundaryFieldRef(), patchID)
    {
      forAll (srnz.boundaryFieldRef()[patchID],facei)
	{
	  if(S0_.value() > srnz.boundaryFieldRef()[patchID][facei])
	    {
	      phiPM_.boundaryFieldRef()[patchID][facei] =
		phiM_.value()+(phiRCP_.value()-phiM_.value())*
		(S0_.value()-srnz.boundaryField()[patchID][facei]);
	    }
	  else
	    {
	      phiPM_.boundaryFieldRef()[patchID][facei] = phiM_.value();
	    }
	}
    }
 */
}
    

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::muJphiJModels::effectivePressureModels::
JohnsonJacksonBoyer::~JohnsonJacksonBoyer()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::muJphiJModels::effectivePressureModels::
JohnsonJacksonBoyer::effectivePressure
(
 const phaseModel& phase,
 const dimensionedScalar& phiM_,
 const dimensionedScalar& nuC_,
 const dimensionedScalar& rhoC_,
 const dimensionedScalar& S0_
 ) const
{
    const volScalarField& alpha_ = phase;
    dimensionedScalar smallval_p("smallvel_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
    const volVectorField& U = phase.U();
    volScalarField srnz("srnz", max(mag(symm(fvc::grad(U))),S0_));
    tmp<volScalarField> tphiPM
      (
       new volScalarField
       (
	IOobject
	(
	 "tphiPM",
	 phase.mesh().time().timeName(),
	 phase.mesh(),
	 IOobject::NO_READ,
	 IOobject::NO_WRITE,
	 false
	 ),
	phase.mesh(),
	dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
        )
       );
    volScalarField& phiPM_ = tphiPM.ref();

    forAll(srnz,ii)
      {
	if(S0_.value() > srnz[ii])
	  {
	    phiPM_[ii] = phiM_.value()+(phiRCP_.value()-phiM_.value())*(S0_.value()-srnz[ii]);
	  }
	else
	  {
	    phiPM_[ii] = phiM_.value();
	  }
      }
    
    
    const fvPatchList& patches = phase.mesh().boundary();
    
    forAll(patches, patchID)
      {
	forAll (srnz.boundaryFieldRef()[patchID],facei)
	  {
	    if(S0_.value() > srnz.boundaryFieldRef()[patchID][facei])
	      {
		phiPM_.boundaryFieldRef()[patchID][facei] =
		  phiM_.value()+(phiRCP_.value()-phiM_.value())*
		  (S0_.value()-srnz.boundaryField()[patchID][facei]);
	      }
	    else
	      {
		phiPM_.boundaryFieldRef()[patchID][facei] = phiM_.value();
	      }
	  }
      }
    
    // Correct coupled BCs
    phiPM_.correctBoundaryConditions();    
    
    return max(
	       aMJ_*max(alpha_-phiRLP_,scalar(0))/
	       max((phiRCP_-alpha_),alphaDeltaMin_)+
	       scalar(2)*nuC_*rhoC_*srnz/
	       pow(phiPM_/max(alpha_-scalar(1),alphaDeltaMin_),2),smallval_p);
}


Foam::tmp<Foam::volScalarField>
Foam::muJphiJModels::effectivePressureModels::
JohnsonJacksonBoyer::effectivePressurePrime
(
 const phaseModel& phase,
 const dimensionedScalar& phiM_,
 const dimensionedScalar& nuC_,
 const dimensionedScalar& rhoC_,
 const dimensionedScalar& S0_
) const
{
    const volScalarField& alpha_ = phase;
    dimensionedScalar smallval_p("smallvel_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
    const volVectorField& U = phase.U();
    volScalarField srnz("srnz", min(mag(symm(fvc::grad(U))),S0_));
    
    tmp<volScalarField> tphiPM
      (
       new volScalarField
       (
        IOobject
        (
         "tphiPM",
         phase.mesh().time().timeName(),
         phase.mesh(),
         IOobject::NO_READ,
         IOobject::NO_WRITE,
         false
         ),
        phase.mesh(),
	dimensionedScalar("one", dimensionSet(0, 0, 0, 0, 0), 1.0)
        )
       );
    volScalarField& phiPM_ = tphiPM.ref();

    forAll(srnz,ii)
      {
      	if(S0_.value() > srnz[ii])
          {
            phiPM_[ii] = phiM_.value()+(phiRCP_.value()-phiM_.value())*(S0_.value()-srnz[ii]);
          }
        else
          {
            phiPM_[ii] = phiM_.value();
          }
      }

    
    const fvPatchList& patches = phase.mesh().boundary();

    forAll(patches, patchID)
      {
	forAll (srnz.boundaryFieldRef()[patchID],facei)
          {
            if(S0_.value() > srnz.boundaryFieldRef()[patchID][facei])
              {
                phiPM_.boundaryFieldRef()[patchID][facei] =
                  phiM_.value()+(phiRCP_.value()-phiM_.value())*
                  (S0_.value()-srnz.boundaryField()[patchID][facei]);
              }
            else
              {
                phiPM_.boundaryFieldRef()[patchID][facei] = phiM_.value();
              }
          }
      }

    // Correct coupled BCs
    phiPM_.correctBoundaryConditions();
    
    
    return max(
               aMJ_*(max(alpha_-phiRLP_,scalar(0)))/
	       pow(max(phiRCP_-alpha_,alphaDeltaMin_),2)-
               scalar(4)*nuC_*rhoC_*srnz*alpha_*phiPM_/
               pow(max(alpha_-phiPM_,alphaDeltaMin_),3),smallval_p);
}



bool Foam::muJphiJModels::effectivePressureModels::
JohnsonJacksonBoyer::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    aMJ_.read(coeffDict_);
    phiRLP_.read(coeffDict_);
    phiRCP_.read(coeffDict_);
    alphaDeltaMin_.read(coeffDict_);
    return true;
}


// ************************************************************************* //
