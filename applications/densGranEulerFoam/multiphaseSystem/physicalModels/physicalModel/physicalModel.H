/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::physicalModel

Description
    Abstract base-class for dispersed-phase particle physical models.

SourceFiles
    physicalModel.C
    newViscosityModel.C

\*---------------------------------------------------------------------------*/

#ifndef physicalModel_H
#define physicalModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class physicalModel Declaration
\*---------------------------------------------------------------------------*/

class physicalModel
{
protected:

    // Protected data

        const dictionary& dict_;
        const phaseModel& phase_;


public:

    //- Runtime type information
    TypeName("physicalModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            physicalModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseModel& phase
            ),
            (dict, phase)
        );


    // Constructors

        physicalModel
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~physicalModel();


    // Selectors

        static autoPtr<physicalModel> New
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    // Member Functions

        //- Return the phase mean physical field
        virtual tmp<volScalarField> calcNu() const = 0;

        virtual tmp<volScalarField> calcGranP() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
