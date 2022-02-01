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

Class
    Foam::viscosityModel

Description

SourceFiles
    viscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef viscosityModel_H
#define viscosityModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{

/*---------------------------------------------------------------------------*\
                           Class viscosityModel Declaration
\*---------------------------------------------------------------------------*/

class viscosityModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        viscosityModel(const viscosityModel&);

        //- Disallow default bitwise assignment
        void operator=(const viscosityModel&);


  

protected:

    // Protected data

  word name_;
  dictionary dict_;

  const volVectorField& U_;
  const surfaceScalarField& phi_;
  const volScalarField& alpha_;
  const volScalarField& rho_;
  const surfaceScalarField& alphaRhoPhi_;
  
public:

    //- Runtime type information
    TypeName("viscosityModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
        (
            autoPtr,
            viscosityModel,
            dictionary,
            (
	     const word& name,
	     const dictionary& dict,
	     const volVectorField& U,
	     const surfaceScalarField& phi,
	     const volScalarField& alpha,
	     const volScalarField& rho,
	     const surfaceScalarField& alphaRhoPhi
	     ),
	    (name, dict, U, phi, alpha, rho, alphaRhoPhi)
	 );


    // Constructors

        //- Construct from components
  viscosityModel
  (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi,
	    const volScalarField& alpha,
            const volScalarField& rho,
            const surfaceScalarField& alphaRhoPhi
  );    //SNA_edit


    // Selectors

        static autoPtr<viscosityModel> New
        (
            const word& name,
            const dictionary& dict,
            const volVectorField& U,
            const surfaceScalarField& phi,
	    const volScalarField& alpha,
             const volScalarField& rho,
             const surfaceScalarField& alphaRhoPhi
	);


    //- Destructor
  
  virtual ~viscosityModel();


    // Member Functions

        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return the strain rate
       virtual tmp<volScalarField> strainRate() const;


virtual tmp<volScalarField> calcNu(const volScalarField& pf, const volScalarField& pa) const = 0;
  
    
        //- Read transportProperties dictionary
 virtual bool read()
  {
    return true;
  }
};

  //  tmp<volScalarField> strainRate();    //SNA_edit


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace biomassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
