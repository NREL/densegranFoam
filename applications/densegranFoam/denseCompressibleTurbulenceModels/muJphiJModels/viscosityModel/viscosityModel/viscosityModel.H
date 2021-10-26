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
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace muJphiJModels
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

  const dictionary& dict_;
  const volVectorField& U_;
  const surfaceScalarField& phi_;
  const volScalarField& alpha_;

  

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
	 const dictionary& dict,
	 const volVectorField& U,
	 const surfaceScalarField& phi,
	 const volScalarField& alpha
        ),
        (dict,U,phi,alpha)
    );


    // Constructors

        //- Construct from components
        viscosityModel(
		       const dictionary& dict,
		       const volVectorField& U,
		       const surfaceScalarField& phi,
		       const volScalarField& alpha
		       );


    // Selectors

        static autoPtr<viscosityModel> New
        (
	 const dictionary& dict,
	 const volVectorField& U,
	 const surfaceScalarField& phi,
	 const volScalarField& alpha
        );


    //- Destructor
    virtual ~viscosityModel();


    // Member Functions
        virtual tmp<volScalarField> nu
        (
	 const volScalarField& muJ,
	 const dimensionedScalar& rho1,
	 const volScalarField& srnz,
	 const volScalarField& ps,
	 const volScalarField& da,
	 const volScalarField& g_nlgf
	 ) const = 0;

        virtual bool read()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace muJphiJModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
