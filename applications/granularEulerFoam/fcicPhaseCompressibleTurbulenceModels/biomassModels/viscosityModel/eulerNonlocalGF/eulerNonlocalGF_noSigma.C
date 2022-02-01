/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.
    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU_.General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    OpenFOAM is distributed in the hope that it will be useful, but WITHOU_.
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICU_.AR PU_.POSE.  See the GNU_.General Public License
    for more details.
    You should have received a copy of the GNU_.General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "eulerNonlocalGF.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"
#include "fvCFD.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
namespace viscosityModels
{
    defineTypeNameAndDebug(eulerNonlocalGF, 0);
    addToRunTimeSelectionTable(viscosityModel, eulerNonlocalGF, dictionary);
}
}
}

Foam::tmp<Foam::volScalarField>
Foam::biomassModels::viscosityModels::eulerNonlocalGF::calcNu() const
{
    dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-9);
    dimensionedScalar smallval_g("smallval_g", dimless/dimTime, 1e-9);
    dimensionedScalar p0("p_offset", dimMass/dimLength/dimTime/dimTime, pres_offset_);


    dimensionedScalar b_nlgf = delfc_/in0_;
    dimensionedScalar mu2 = delfc_+fcs_;

    dimensionedScalar gravity_acc = mag(g_);
    double xref, zref, yc, yref, dy;
    double trapsum = 0.0;
 
    const volVectorField& C = U_.mesh().C();

     volScalarField alpha_h("alpha_h", alpha_*href_nlgf);

    alpha_h *= scalar(0.0);

    int cellid, cellid0;
    
    Info<< "Calculating field g.h_NLGF\n" << endl;    
    forAll(C,i)
      {
	xref = C[i].x();
	yref = C[i].y();
	zref = C[i].z();
	dy = (href_nlgf[i]-yref)/(nNormalCell-1);
	yc = yref;
	trapsum = 0.0;

	for( int ii = 0; ii < nNormalCell; ii++)
	  {
	    if(alpha_[i]<alphaBoundary){
              break;
            }
	    
	    yc = yc + ii*dy;
	    vector pointloc(xref,yc,zref);
	    cellid = U_.mesh().findCell(pointloc);
	    if(ii == 0){
	      cellid0 = cellid;
	    }
	    else{
	      trapsum += 0.5*(alpha_[cellid]+alpha_[cellid0])*dy;
	      cellid0 = cellid;
	    }
	  }
	alpha_h[i] = trapsum;
      }
    Info<< "Ended Calculating field g.h_NLGF\n" << endl;
    
    volScalarField pV("pV", (0.5)*alpha_*rho_*magSqr(U_));
    volScalarField rho_gh("rho_gh", alpha_h*gravity_acc*rho_);

    dimensionedScalar smallp("smallp", dimMass/dimLength/dimTime/dimTime, 1e0);
      
    volScalarField pres_reg
      (
       max(rho_gh+pV+p0,smallp)
       );
    // strain-rate adjusted for zero components
    // note:VSMALL in openfoam is 1e-37 (HS: Aug 28, 2018)
    // maybe do not need to limit strain rate? JJS 10/5/18
    /*
    volSymmTensorField sigmaD
      (
       "sigmaD",
       rho_*(2.0*sigNu_*symm(fvc::grad(U_)) - (2.0/3.0)*sigNu_*tr(symm(fvc::grad(U_)))*symmTensor::I)
       );
    */
    volScalarField srnz
        (
         "srnz",
         strainRate()
        );

    volScalarField gnz
        (
         "gnz",
         max(g_nlgf, smallval_g)
        );

      
    if(U_.time().outputTime())
      {
        pres_reg.write();
	srnz.write();
	gnz.write();
	pV.write();
	//	sigmaD.write();
	rho_gh.write();
      }
    
    //mu is set in a way that it is between mus and mu2 -- limiting those
    //bounds by another fraction (10%)
    //dimensionedScalar delmu = scalar(0.1)*delfc_;
    dimensionedScalar delmu = scalar(0.05)*delfc_;
    volScalarField mu
        (
         "mu",
         //max(min(srnz/gnz, mu2-delmu), fcs_+delmu)
         // do not limit lower value, JJS 10/5/18
         min(srnz/gnz, mu2-delmu)
        );

    Info << "GNLGF eq" << endl;
    fvScalarMatrix gEqn
      (fvm::ddt(g_nlgf) + fvm::div(phi_, g_nlgf) + fvm::SuSp(-fvc::div(phi_), g_nlgf)   
       - fvm::laplacian(A_*A_*graind_*graind_*t0inv_, g_nlgf) ==
       - fvm::Sp(t0inv_*delfc_*(fcs_*g_nlgf - srnz)/(mu2*g_nlgf - srnz), g_nlgf)
       -fvm::Sp(t0inv_*b_nlgf*graind_*srnz * sqrt(rho_/pres_reg), g_nlgf)
       );
    
    gEqn.relax();
    gEqn.solve();
    
    const volScalarField etasolid
      (
       "etasolid",
       pres_reg / max(g_nlgf,smallval_g)
       );
    
    bool debug=true;
    
    //if(debug)
    //{
        const volScalarField justone
            (
             "justone",pres_reg/pres_reg
            );
        double volume=fvc::domainIntegrate(justone).value();
     Info << min(rho_gh) <<"\n";
     Info << max(rho_gh) <<"\n";
     Info << min(pV) <<"\n";
     Info << max(pV) <<"\n";
     Info << min(etasolid)<<"\n";
     Info << max(etasolid)<<"\n";
     Info << "average(etasolid_):"<<fvc::domainIntegrate(etasolid).value()/volume<<"\n";
     Info << min(pres_reg)<<"\n";
     Info << max(pres_reg)<<"\n";
     Info << "average(pres_reg):"<<fvc::domainIntegrate(pres_reg).value()/volume<<"\n";
     Info << min(g_nlgf)<<"\n";
     Info << max(g_nlgf)<<"\n";
     Info << min(mu)<<"\n";
     Info << max(mu)<<"\n";
     Info << max(href_nlgf) << "\n";
     Info << min(href_nlgf) << "\n";
     //}
    
    return
      (
       min(nu0_, etasolid/rho_)//nu0_/justone//min(nu0_, etasolid/rho_)
       );

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::viscosityModels::eulerNonlocalGF::eulerNonlocalGF
(
 const word& name,
 const dictionary& viscosityProperties,
 const volVectorField& U,
 const surfaceScalarField& phi,
 const volScalarField& alpha,
 const volScalarField& rho,
 const surfaceScalarField& alphaRhoPhi
 )
  :
  viscosityModel(name, viscosityProperties, U, phi, alpha, rho, alphaRhoPhi),
  eulerNonlocalGFCoeffs_
  (
   viscosityProperties.optionalSubDict(typeName + "Coeffs")
   ),
  graind_("graind", dimLength, eulerNonlocalGFCoeffs_),
  fcs_("fcs", dimless, eulerNonlocalGFCoeffs_),
  delfc_("delfc", dimless, eulerNonlocalGFCoeffs_),
  in0_("in0", dimless, eulerNonlocalGFCoeffs_),
  nu0_("nu0", dimViscosity, eulerNonlocalGFCoeffs_),
  rho_("rho", dimDensity, eulerNonlocalGFCoeffs_),
  A_("A_nlgf", dimless, eulerNonlocalGFCoeffs_),
  t0inv_("t0inv_nlgf", dimless/dimTime, eulerNonlocalGFCoeffs_),
  alphaBoundary(eulerNonlocalGFCoeffs_.lookupOrDefault<scalar>("alphaThreshold",0.01)),
  nNormalCell(eulerNonlocalGFCoeffs_.lookupOrDefault<scalar>("nCell",10)),
  pres_offset_(eulerNonlocalGFCoeffs_.lookupOrDefault<scalar>("pres_offset",0.0)),
  href_nlgf
  (
   IOobject
   (
    "href_nlgf",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
   U_.mesh()
   ),
  g_nlgf
  (
   IOobject
   (
    "g_nlgf",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
    ),
   U_.mesh()
   ),
  gsrc_expl
  (
   IOobject
   (
    "gsrc_expl",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ), U_.mesh(), dimensionedScalar("gsrc_expl", dimless/dimTime/dimTime, scalar(0.0))
   ),
  gdcoeff
  (
   IOobject
   (
    "gdcoeff",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ), U_.mesh(), dimensionedScalar("gdcoeff", dimLength*dimLength/dimTime, scalar(1.0))
   ),
  gsrc_impl
  (
   IOobject
   (
    "gsrc_impl",
    U_.time().timeName(),
    U_.mesh(),
    IOobject::NO_READ,
    IOobject::NO_WRITE
    ), U_.mesh(), dimensionedScalar("gsrc_impl", dimless/dimTime, scalar(0.0))
   ),
  g_(U_.mesh().lookupObject<uniformDimensionedVectorField>("g")),
  href_(U.db().lookupObject<uniformDimensionedScalarField>("hRef")),
  nu_
  (
   IOobject
   (
    name,
    U_.time().timeName(),
    U_.db(),
    IOobject::NO_READ,
    IOobject::AUTO_WRITE
    ),
   calcNu()
   )
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::biomassModels::viscosityModels::eulerNonlocalGF::~eulerNonlocalGF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
Foam::tmp<Foam::volScalarField>
Foam::biomassModels::viscosityModels::eulerNonlocalGF::strainRate() const
{
  return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}
*/


bool Foam::biomassModels::viscosityModels::eulerNonlocalGF::read
(
)
{

    eulerNonlocalGFCoeffs_ <<= dict_.optionalSubDict(typeName + "Coeffs");


    eulerNonlocalGFCoeffs_.lookup("graind") >> graind_;
    eulerNonlocalGFCoeffs_.lookup("fcs") >> fcs_;
    eulerNonlocalGFCoeffs_.lookup("delfc") >> delfc_;
    eulerNonlocalGFCoeffs_.lookup("in0") >> in0_;
    eulerNonlocalGFCoeffs_.lookup("nu0") >> nu0_;
    eulerNonlocalGFCoeffs_.lookup("A_nlgf") >> A_;
    eulerNonlocalGFCoeffs_.lookup("t0inv_nlgf") >> t0inv_;
    eulerNonlocalGFCoeffs_.lookup("alphaThreshold") >> alphaBoundary;
    eulerNonlocalGFCoeffs_.lookup("nCell") >> nNormalCell;
    eulerNonlocalGFCoeffs_.lookup("rho") >> rho_;

    return true;
}

// ************************************************************************* //
