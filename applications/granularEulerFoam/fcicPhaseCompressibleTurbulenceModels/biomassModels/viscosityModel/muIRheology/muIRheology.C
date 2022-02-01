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

#include "muIRheology.H"
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
    defineTypeNameAndDebug(muIRheology, 0);
    addToRunTimeSelectionTable(viscosityModel, muIRheology, dictionary);
}
}
}

Foam::tmp<Foam::volScalarField>
Foam::biomassModels::viscosityModels::muIRheology::calcNu(const volScalarField& pf, const volScalarField& pa) const
{
    //need to play with the tolerances (HS: Aug 28, 2018)
    dimensionedScalar smallval_sr("smallval_sr", dimless/dimTime, 1e-8);

    // strain-rate adjusted for zero components
    // note:VSMALL in openfoam is 1e-37 (HS: Aug 28, 2018)
    /*    const volScalarField srnz
        (
         "srnz",
         //max(strainRate(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL))
         max(strainRate(), smallval_sr)
        );
    */
    // do we need "const" for all these temporary variables? does it matter? JJS 8/24/18
    //
    // HS (aug 28, 2018)
    // const speeds up computation because the compiler now knows there are no updates going to happen to 
    // this variable in memory.
    // const will also ensure we don't change the variable accidentally (debugging)
    //
    //HS (aug 28, 2018)
    // gravitational pressure (note initialize href in const)
    // also need to be fixed to account for varying volume fractions along gravity directions
    // ideally it should be integral(alpha rho g dy)
    /*volScalarField presgrav
      (
      "presgrav",
      rho_*gh_
      );
    // total pressure
    // HS (aug 28, 2018)
    // we need to be smart here - pres_ really has no meaning 
    // in incompressible flow, only grad(pres_) is important as it changes velocity
    // therefore we need some sort of a reference pressure - we could initialize
    // the domain with atmospheric pressure instead of 0
    volScalarField presT
    (
    "presT",
    pres_ + presgrav
    );*/
    // prevent divide-by-zero in next formula
    // improve this so that the sign is correct? 
    /*volScalarField presTnz
      (
      "presTnz",
      max(mag(presT),
      dimensionedScalar ("VSMALL", dimMass/dimLength/dimTime/dimTime, VSMALL))
      );*/
    //Info << min(mag(presTnz)) << "\n";

    //regularized pressure
    volScalarField pres_reg
        (
         "pres_reg",
         max(pres_,dimensionedScalar ("smallp", dimMass/dimLength/dimTime/dimTime, 1e1))
        );

    // inertial number "I" -- careful of divide by zero!
    const volScalarField inrtlnum
        (
         "inrtlnum",
         srnz_*graind_ /
         sqrt( pres_reg/rho_ )
        );
    // friction coefficient "mu" -- also potential divide by zero
    const volScalarField friccoef
        (
         "friccoef",
         fcs_ + delfc_/(scalar(1) + in0_/inrtlnum)
        );
    // the bulk mu-I rheology viscosity
    const volScalarField etasolid
        (
         "etasolid",
         friccoef*pres_reg / srnz_
        );

    // // temporary test of functionality of dimensions dynamic viscosity
    // dimensionedScalar var("var", dimMass/dimLength, 1.0);
    // volScalarField etasolid
    //   (
    //    "etasolid",
    //    var*max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL))
    //    );
    // // temporary test of functionality of dimensions kinematic viscosity
    // dimensionedScalar var("var", dimLength*dimLength, 1.0);
    // volScalarField nutemp
    //   (
    //    "nutemp",
    //    var*max(sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL))
    //    );


    const volScalarField justone
        (
         "justone",pres_reg/pres_reg
        );
    double volume=fvc::domainIntegrate(justone).value();

    Info << min(etasolid)<<"\n";
    Info << max(etasolid)<<"\n";

    Info << min(friccoef)<<"\n";
    Info << max(friccoef)<<"\n";

    Info << min(inrtlnum)<<"\n";
    Info << max(inrtlnum)<<"\n";

    Info << min(pres_reg)<<"\n";
    Info << max(pres_reg)<<"\n";

    Info << "average(pres_reg):"<<fvc::domainIntegrate(pres_reg).value()/volume<<"\n";
    Info << "average(friccoef_):"<<fvc::domainIntegrate(friccoef).value()/volume<<"\n";
    Info << "average(inrtlnum_):"<<fvc::domainIntegrate(inrtlnum).value()/volume<<"\n";
    Info << "average(etasolid_):"<<fvc::domainIntegrate(etasolid).value()/volume<<"\n";

    return
        (
         min(nu0_, etasolid/rho_)
         //min(nu0_, nutemp)
        );
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::viscosityModels::muIRheology::muIRheology
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
  muIRheologyCoeffs_
  (
   viscosityProperties.optionalSubDict(typeName + "Coeffs")
   ),
        graind_("graind", dimLength, muIRheologyCoeffs_),
        fcs_("fcs", dimless, muIRheologyCoeffs_),
        delfc_("delfc", dimless, muIRheologyCoeffs_),
        in0_("in0", dimless, muIRheologyCoeffs_),
        nu0_("nu0", dimViscosity, muIRheologyCoeffs_),
        rho_("rho", dimDensity, viscosityProperties),
        pres_(U.db().lookupObject<volScalarField>("p_rgh")),
        gh_(U.db().lookupObject<volScalarField>("gh")),
        srnz_(U.mesh().lookupObject<volScalarField>("srnz"))/*,
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

Foam::biomassModels::viscosityModels::muIRheology::~muIRheology()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
Foam::tmp<Foam::volScalarField>
Foam::biomassModels::viscosityModels::muIRheology::strainRate() const
{
  return sqrt(2.0)*mag(symm(fvc::grad(U_)));
}
*/


bool Foam::biomassModels::viscosityModels::muIRheology::read
(
)
{

    muIRheologyCoeffs_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    muIRheologyCoeffs_.lookup("graind") >> graind_;
    muIRheologyCoeffs_.lookup("fcs") >> fcs_;
    muIRheologyCoeffs_.lookup("delfc") >> delfc_;
    muIRheologyCoeffs_.lookup("in0") >> in0_;
    muIRheologyCoeffs_.lookup("nu0") >> nu0_;

    return true;
}

// ************************************************************************* //
