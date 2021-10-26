/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::muJphiJModels::effectivePressureModel

SourceFiles
    effectivePressureModel.C

\*---------------------------------------------------------------------------*/

#ifndef effectivePressureModel_H
#define effectivePressureModel_H

#include "dictionary.H"
#include "volFields.H"
#include "dimensionedTypes.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace muJphiJModels
{

/*---------------------------------------------------------------------------*\
                           Class effectivePressureModel Declaration
\*---------------------------------------------------------------------------*/

class effectivePressureModel
{
    // Private member functions

        //- Disallow default bitwise copy construct
        effectivePressureModel(const effectivePressureModel&);

        //- Disallow default bitwise assignment
        void operator=(const effectivePressureModel&);


protected:

    // Protected data

        //- Reference to higher-level dictionary for re-read
        const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("effectivePressureModel");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        effectivePressureModel,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );


    // Constructors

        //- Construct from components
        effectivePressureModel(const dictionary& dict);


    // Selectors

        static autoPtr<effectivePressureModel> New
        (
            const dictionary& dict
        );


    //- Destructor
    virtual ~effectivePressureModel();


    // Member Functions

        virtual tmp<volScalarField> effectivePressure
        (
	 const phaseModel& phase,
	 const dimensionedScalar& phiM_,
	 const dimensionedScalar& nuC_,
	 const dimensionedScalar& rhoC_,
	 const dimensionedScalar& S0_
	 ) const = 0;

        virtual tmp<volScalarField> effectivePressurePrime
        (
	 const phaseModel& phase,
         const dimensionedScalar& phiM_,
	 const dimensionedScalar& nuC_,
	 const dimensionedScalar& rhoC_,
	 const dimensionedScalar& S0_
	 ) const = 0;

        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace muJphiJModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
