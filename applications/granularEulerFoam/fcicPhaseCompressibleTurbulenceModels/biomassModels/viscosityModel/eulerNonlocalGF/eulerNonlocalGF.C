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
Foam::biomassModels::viscosityModels::eulerNonlocalGF::calcNu(const volScalarField& pf, const volScalarField& pa) const
{
    dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-8);
    dimensionedScalar smallval_g("smallval_g", dimless/dimTime, 1e-8);
    dimensionedScalar p0("p_offset", dimMass/dimLength/dimTime/dimTime, pres_offset_);


    dimensionedScalar b_nlgf = delfc_/in0_;
    dimensionedScalar mu2 = delfc_+fcs_;

    dimensionedScalar smallp("smallp", dimMass/dimLength/dimTime/dimTime, 1e0);

    volScalarField pres_("pres_", U_.mesh().lookupObject<volScalarField>("p_IR")); 

    
    
    volScalarField pV("pV", (0.5)*rho_*magSqr(U_));

        
    volScalarField srnz
        (
         "srnz",
         strainRate()
        );

    volScalarField tau_inv_par = relaxPa_*alpha_*srnz;

    /*
      volScalarField granP
    (
     "granP",
     rhop_*pow(delfc_*graind_*srnz_*alpha_/in0_/max(alphaMax_-alpha_,smallal),2)
     );
  volScalarField granPch
    (
     "granPch",
     rhop_*pow(B_phi_*graind_*srnz_*alpha_/max(alphaMax_-alpha_,smallal),2)
     );
    */


    fvScalarMatrix ppresEqn
      (
       fvm::ddt(ppres)
       + fvm::div(phi_, ppres)
       - fvm::Sp(fvc::div(phi_), ppres)
       ==
       tau_inv_par*(pa+pf)
       -fvm::Sp(tau_inv_par, ppres)
       );
    ppresEqn.relax();
    ppresEqn.solve();
    
    volScalarField pres_reg
      (
       "pres_reg",
       max(ppres,smallp)
       );
    
    if(U_.time().outputTime())
    {
      pres_reg.write();
    }


    volScalarField gnz
        (
         "gnz",
         max(g_nlgf, smallval_g)
        );

    dimensionedScalar delmu = scalar(0.05)*delfc_;
    volScalarField mu
        (
         "mu",
         //max(min(srnz/gnz, mu2-delmu), fcs_+delmu)
         // do not limit lower value, JJS 10/5/18
         min(srnz/gnz, mu2-delmu)
	 //srnz/gnz
        );

    Info << "GNLGF eq" << endl;
    /*    
    fvScalarMatrix gEqn
      (fvm::ddt(g_nlgf) + fvm::div(phi_, g_nlgf)// + fvm::SuSp(-fvc::div(phi_), g_nlgf) ==   
       - fvm::laplacian(A_*A_*graind_*graind_*t0inv_, g_nlgf) == 
       - fvm::Sp(t0inv_*delfc_*(fcs_*g_nlgf - srnz)/(mu2*g_nlgf - srnz), g_nlgf)
       - fvm::Sp(t0inv_*b_nlgf*graind_*srnz * sqrt(rho_/pres_reg), g_nlgf)
       );
    */


    //    fv::options& fvOptions(fv::options::New(U_.mesh()));

    fvScalarMatrix gEqn
      (fvm::ddt(g_nlgf) + fvm::div(phi_, g_nlgf) //+ fvm::SuSp(-fvc::div(phi_), g_nlgf)
       - fvm::laplacian(A_*A_*graind_*graind_*t0inv_, g_nlgf) ==
       - fvm::Sp(t0inv_*delfc_*(fcs_ - mu)/(mu2 - mu), g_nlgf)
       - fvm::Sp(t0inv_*b_nlgf*graind_*srnz * sqrt(rho_/pres_reg), g_nlgf)
       );

    
    gEqn.relax();
    //fvOptions.constrain(gEqn);
    gEqn.solve();
    //fvOptions.correct(g_nlgf);



    //    g_nlgf.min(smallval_g);
    
    const volScalarField etasolid
      (
       "etasolid",
       pres_reg / max(g_nlgf,smallval_g)
       );
    
    bool debug=false;
    
    //if(debug)
    //{
        const volScalarField justone
            (
             "justone",pres_reg/pres_reg
            );
        double volume=fvc::domainIntegrate(justone).value();
	//     Info << min(rho_gh) <<"\n";
	//Info << max(rho_gh) <<"\n";
	//Info << min(pV) <<"\n";
	//Info << max(pV) <<"\n";
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
       min(nu0_, etasolid/rho_)
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
  relaxPa_("relaxppres", dimless, eulerNonlocalGFCoeffs_),
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
  ppres
  (
   IOobject
   (
    "ppres",
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
  //pres_(U.db().lookupObject<volScalarField>("p_IR")),
  //pres_(U.db().lookupObject<volScalarField>("p_rgh")),
  g_(U_.mesh().lookupObject<uniformDimensionedVectorField>("g")),
  //  srnz(U.mesh().lookupObject<volScalarField>("srnz")),
  href_(U.db().lookupObject<uniformDimensionedScalarField>("hRef"))/*,
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
   )*/
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
    eulerNonlocalGFCoeffs_.lookup("relaxppres") >> relaxPa_;

    return true;
}

// ************************************************************************* //
