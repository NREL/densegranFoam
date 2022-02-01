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
    Foam::biomassModels::fcicGranularPressureModel

SourceFiles
    fcicGranularPressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef fcicGranularPressureModel_H
#define fcicGranularPressureModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{

/*---------------------------------------------------------------------------*\
                           Class fcicGranularPressureModel Declaration
\*---------------------------------------------------------------------------*/

class fcicGranularPressureModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        fcicGranularPressureModel(const fcicGranularPressureModel&);

        //- Disallow default bitwise assignment
        void operator=(const fcicGranularPressureModel&);


protected:

    // Protected data

  const dictionary& dict_;
  const volVectorField& U_;
  const volScalarField& alpha_;
  const volScalarField& rho_;


public:

    //- Runtime type information
    TypeName("fcicGranularPressureModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        fcicGranularPressureModel,
        dictionary,
        (
	 const dictionary& dict,
	 const volVectorField& U,
	 const volScalarField& alpha,
	 const volScalarField& rho
        ),
        (dict, U, alpha, rho)
    );


    // Constructors

        //- Construct from components
        fcicGranularPressureModel
	(
	  const dictionary& dict,
	  const volVectorField& U,
	  const volScalarField& alpha,
	  const volScalarField& rho
	 );
  

    // Selectors

        static autoPtr<fcicGranularPressureModel> New
        (
             const dictionary& dict,
	     const volVectorField& U,
	     const volScalarField& alpha,
	     const volScalarField& rho
        );
  
        //- Return the strain rate
       virtual tmp<volScalarField> strainRate() const;

    //- Destructor
    virtual ~fcicGranularPressureModel();


    // Member Functions


        //- Derivative of the granular pressure coefficient
        virtual tmp<volScalarField> granularPressureCoeff
        () const = 0;

        virtual bool read()
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace biomassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
