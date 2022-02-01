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

\*---------------------------------------------------------------------------*/

#include "fcicGranularPressureModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::biomassModels::fcicGranularPressureModel>
Foam::biomassModels::fcicGranularPressureModel::New
(
 const dictionary& dict,
 const volVectorField& U,
 const volScalarField& alpha,
 const volScalarField& rho
)
{
    word fcicGranularPressureModelType(dict.lookup("fcicGranularPressureModel"));

    Info<< "Selecting fcicGranularPressureModel "
        << fcicGranularPressureModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fcicGranularPressureModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalError
            << "fcicGranularPressureModel::New(const dictionary&) : " << endl
            << "    unknown fcicGranularPressureModelType type "
            << fcicGranularPressureModelType
            << ", constructor not in hash table" << endl << endl
            << "    Valid fcicGranularPressureModelType types are :" << endl;
        Info<< dictionaryConstructorTablePtr_->sortedToc()
            << abort(FatalError);
    }

    return autoPtr<fcicGranularPressureModel>(cstrIter()(dict, U, alpha, rho));
}


// ************************************************************************* //
