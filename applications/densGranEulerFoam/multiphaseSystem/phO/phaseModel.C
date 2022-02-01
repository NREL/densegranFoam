/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "phaseModel.H"
#include "fvOptions.H"
#include "diameterModel.H"
#include "fixedValueFvPatchFields.H"
#include "slipFvPatchFields.H"
#include "partialSlipFvPatchFields.H"
#include "surfaceInterpolate.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseModel::phaseModel
(
    const word& phaseName,
    const dictionary& phaseDict,
    const fvMesh& mesh
)
:
    volScalarField
    (
        IOobject
        (
            IOobject::groupName("alpha", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    name_(phaseName),
    phaseDict_(phaseDict),
    nu_
    (
        "nu",
        dimensionSet(0, 2, -1, 0, 0),
        phaseDict_
    ),

    rhoC_("rhoC",dimensionSet(1, -3, 0, 0, 0), phaseDict_),
    nuC_("nuC",dimensionSet(0, 2, -1, 0, 0), phaseDict_),
    mu1_("mu1",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    mu2_("mu2",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    J0_("J0",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    rhoS_("rhoS",dimensionSet(1, -3, 0, 0, 0), phaseDict_),
    phiM_("phiM",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    S0_("S0",dimensionSet(0, 0, -1, 0, 0), phaseDict_),
    
    
    kappa_
    (
        "kappa",
        dimensionSet(1, 1, -3, -1, 0),
        phaseDict_
    ),
    Cp_
    (
        "Cp",
        dimensionSet(0, 2, -2, -1, 0),
        phaseDict_
    ),
    
    aMJ_("aMJ",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    phiRLP_("phiRLP",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    phiRCP_("phiRCP",dimensionSet(0, 0, 0, 0, 0), phaseDict_),
    alphaDeltaMin_("alphaDeltaMin",dimensionSet(0, 0, 0, 0, 0), phaseDict_),


    U_
    (
        IOobject
        (
            IOobject::groupName("U", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    DDtU_
    (
        IOobject
        (
            IOobject::groupName("DDtU", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("0", dimVelocity/dimTime, Zero)
    ),
    alphaPhi_
    (
        IOobject
        (
            IOobject::groupName("alphaPhi", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
    )
{
    const word phiName = IOobject::groupName("phi", name_);

    IOobject phiHeader
    (
        phiName,
        mesh.time().timeName(),
        mesh,
        IOobject::NO_READ
    );

    if (phiHeader.typeHeaderOk<surfaceScalarField>(true))
    {
        Info<< "Reading face flux field " << phiName << endl;

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field " << phiName << endl;

        wordList phiTypes
        (
            U_.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(U_.boundaryField(), i)
        {
            if
            (
                isA<fixedValueFvPatchVectorField>(U_.boundaryField()[i])
             || isA<slipFvPatchVectorField>(U_.boundaryField()[i])
             || isA<partialSlipFvPatchVectorField>(U_.boundaryField()[i])
            )
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiPtr_.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    phiName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::flux(U_),
                phiTypes
            )
        );
    }

    dPtr_ = diameterModel::New
    (
        phaseDict_,
        *this
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseModel::~phaseModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseModel> Foam::phaseModel::clone() const
{
    NotImplemented;
    return autoPtr<phaseModel>(nullptr);
}



void Foam::phaseModel::correct(volScalarField& alpha_)
{
    // Local references
    volScalarField& da = dPtr_->d();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));
    dimensionedScalar smallval_sr("smallvel_sr", dimless/dimTime, scalar(1e-6));
    dimensionedScalar smallval_p("smallvel_p", dimMass/dimLength/dimTime/dimTime, scalar(1e-6));
    

    
    Info<<"cf.1\n";
    volScalarField srnz("srnz", min(mag(D)-tr(gradU)/scalar(3),S0_));
    Info<<"cf.2\n";

    tmp<volScalarField> tphiPM
      (
       new volScalarField
       (
	IOobject
	(
	 "tphiPM",
	 U_.mesh().time().timeName(),
	 U_.mesh(),
	 IOobject::NO_READ,
	 IOobject::NO_WRITE,
	 false
	 ),
	U_.mesh(),
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
    
    
    const fvPatchList& patches = U_.mesh().boundary();
    
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
    
    volScalarField pS_ ("pS_", max(
				   aMJ_.value()*max(alpha_-phiRLP_,scalar(0))/
	       max((phiRCP_-alpha_),alphaDeltaMin_)+
	       scalar(2)*nuC_*rhoC_*srnz/
	       pow(phiPM_/max(alpha_-scalar(1),alphaDeltaMin_),2),smallval_p));
    
    volScalarField J1_ = scalar(2)*srnz*nuC_.value()*rhoC_.value()/max(pS_,smallval_p);
    Info<<"cf.9\n";
    volScalarField muJ_ = mu1_.value()+(mu2_.value()-mu1_.value())*J1_/max(J1_+J0_,scalar(1e-6))+
      J1_+scalar(5)/scalar(2)*phiM_.value()*sqrt(J1_);
    Info<<"cf.10\n";
    
    nut_ = muJ*ps/2/rhoS_/srnz;
    Info<<"cf.11\n";
    
    volScalarField pS_mjpj
      ("pS_mjpj", pS_);
    
  if(U_.time().outputTime())
    {
      pS_mjpj.write();
      srnz.write();
    } 
}


bool Foam::phaseModel::read(const dictionary& phaseDict)
{
    phaseDict_ = phaseDict;

    //if (nuModel_->read(phaseDict_))
    {
      phaseDict_.lookup("rhoC") >> rhoC_.value();
      phaseDict_.lookup("nuC") >> nuC_.value();
      phaseDict_.lookup("mu1") >> mu1_.value();
      phaseDict_.lookup("mu2") >> mu2_.value();
      phaseDict_.lookup("J0") >> J0_.value();
      phaseDict_.lookup("rhoS") >> rhoS_.value();
      phaseDict_.lookup("phiM") >> phiM_.value();
      phaseDict_.lookup("S0") >> S0_.value();
      phaseDict_.lookup("kappa") >> kappa_.value();
      phaseDict_.lookup("Cp") >> Cp_.value();
      phaseDict_.lookup("aMJ") >> aMJ_.value();
      phaseDict_.lookup("phiRLP") >> phiRLP_.value();
      phaseDict_.lookup("phiRCP") >> phiRCP_.value();
      phaseDict_.lookup("alphaDeltaMin") >> alphaDeltaMin_.value();
      
      return true;
    }
    // else
    // {
    //     return false;
    // }

    return true;
}


void Foam::phaseModel::correctInflowOutflow(surfaceScalarField& alphaPhi) const
{
    surfaceScalarField::Boundary& alphaPhiBf = alphaPhi.boundaryFieldRef();
    const volScalarField::Boundary& alphaBf = boundaryField();
    const surfaceScalarField::Boundary& phiBf = phi().boundaryField();

    forAll(alphaPhiBf, patchi)
    {
        fvsPatchScalarField& alphaPhip = alphaPhiBf[patchi];

        if (!alphaPhip.coupled())
        {
            alphaPhip = phiBf[patchi]*alphaBf[patchi];
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::phaseModel::d() const
{
    return dPtr_().d();
}

// ************************************************************************* //
