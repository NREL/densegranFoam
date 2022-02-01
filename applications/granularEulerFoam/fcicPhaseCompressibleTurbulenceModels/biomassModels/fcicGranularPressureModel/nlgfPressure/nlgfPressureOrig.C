/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
0;95;0c   \\    /   O peration     |
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

#include "nlgfPressure.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
namespace fcicGranularPressureModels
{
    defineTypeNameAndDebug(nlgf, 0);

    addToRunTimeSelectionTable
    (
        fcicGranularPressureModel,
        nlgf,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModels::nlgf::~nlgf()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



Foam::tmp<Foam::volScalarField>
Foam::biomassModels::fcicGranularPressureModels::nlgf::
granularPressureCoeffPrime
() const
{
  dimensionedScalar p0("p_offset", dimMass/dimLength/dimTime/dimTime, pres_offset_);
  dimensionedScalar gravity_acc = mag(g_);
  double xref, zref, yc, yref, dy;
  double trapsum = 0.0;
  
  const volVectorField& C = U_.mesh().C();
  
  volScalarField alpha_h("alpha_h", alpha_*href_nlgf);

  volScalarField alpha_my("alpha_my", alpha_*href_nlgf);
  alpha_h *= scalar(0.0);
  //alpha_my *= scalar(0.0);

  Info << max(rho_) << endl;
  Info << min(rho_) << endl;

  //  Info << C.size() << endl;
  
  int cellid;
  int cellid0;
  
  Info<< "Calculating field g.h_NLGF_pressure\n" << endl;
    forAll(C,i)
      {
	//	Info << "Entered forAll" << endl;
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
	    //  Info << "Entered after Break" << endl;
	    //	    alpha_my[i] = alpha_[ii];
            yc = yc + ii*dy;
            vector pointloc(xref,yc,zref);
            cellid = U_.mesh().findCell(pointloc);
	    //	    Info << cellid << endl;
            if(ii == 0){
              cellid0 = cellid;
	      
              //Info << "ii=0" << endl;
            }
            else{
              if(cellid < 0)
                {
		  //  Info << "ii<0" << endl;
                  continue;
                }
	      // Info << "Entered trapz calc" << endl;
              trapsum += 0.5*(alpha_[cellid]+alpha_[cellid0])*dy;
              cellid0 = cellid;
            }
          }
      alpha_h[i] = trapsum;
    }

  Info<< "Ended field g.h_NLGF_pressure\n" << endl;

  volScalarField pVP("pVP", (0.5)*alpha_*rho_*magSqr(U_));
  volScalarField rho_ghP("rho_ghP", alpha_h*gravity_acc*rho_);
  
  dimensionedScalar smallp("smallp", dimMass/dimLength/dimTime/dimTime, 1e0);

  //  volScalarField& pres_regN = const_cast<volScalarField&>(pres_reg);
  /*  forAll(C,i)
    {
      pres_reg[i] = max(rho_gh[i]+pV[i]+p0,smallp)
    }
  */
  
  volScalarField pres_regP
    (
     "pres_regP",
     max(rho_ghP+p0,smallp)//dimensionedScalar ("smallp", dimMass/dimLength/dimTime/dimTime, 1e1))
     );

  //pres_regN = (max(rho_gh+pV+p0,smallp));//dimensionedScalar ("smallp", dimMass/dimLength/dimTime/dimTime, 1e1));

      if(U_.time().outputTime())
      {
        pres_regP.write();
    	//      srnz.write();
        //U_.write();
        pVP.write();
        alpha_my.write();
        rho_ghP.write();
        //Cmesh.write();
        //alphaCheck.write();
        
      }
     Info << min(max(pres_regP,smallp))<<"\n";
     Info << max(max(pres_regP,smallp))<<"\n";
  return max(pres_regP,smallp);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::biomassModels::fcicGranularPressureModels::nlgf::nlgf
(
 const dictionary& dict,
 const volVectorField& U,
 const volScalarField& alpha,
 const volScalarField& rho
 )
  :
    fcicGranularPressureModel(dict, U, alpha, rho),
    coeffDict_(dict.optionalSubDict(typeName + "Coeffs")),
    pres_(U.db().lookupObject<volScalarField>("p")),
    g_(U_.mesh().lookupObject<uniformDimensionedVectorField>("g")),
    pres_offset_(coeffDict_.lookupOrDefault<scalar>("pres_offset",0.0)),
    href_nlgf(U.mesh().lookupObject<volScalarField>("href_nlgf")),
    alphaBoundary(coeffDict_.lookupOrDefault<scalar>("alphaThreshold",0.01)),
    nNormalCell(coeffDict_.lookupOrDefault<scalar>("nCell",10)),
    //    pres_reg(U.mesh().lookupObject<volScalarField>("pres_reg")),
    p_rgh(U.mesh().lookupObject<volScalarField>("p")),
    href_(U.db().lookupObject<uniformDimensionedScalarField>("hRef"))
{}

bool Foam::biomassModels::fcicGranularPressureModels::nlgf::read()
{
    coeffDict_ <<= dict_.optionalSubDict(typeName + "Coeffs");

    
    coeffDict_.lookup("alphaThreshold") >> alphaBoundary;
    coeffDict_.lookup("nCell") >> nNormalCell;
    return true;
}
// ************************************************************************* //
