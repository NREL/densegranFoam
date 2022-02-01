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

\*---------------------------------------------------------------------------*/
        // Calculating exceeding volume fractions
        volScalarField alpha1Ex(max(alpha1 - alphaMax, scalar(0)));

        // Finding neighbouring cells of the whole domain
labelListList neighbour = mesh().cellCells();
scalarField cellVolumes(mesh().cellVolumes());

        forAll (alpha1Ex, celli)
        {
            // Finding the labels of the neighbouring cells
            labelList neighbourCell = neighbour[celli];

            // Initializing neighbouring cells contribution
            scalar neighboursEx = 0.0;

            forAll (neighbourCell, cellj)
            {
                labelList neighboursNeighbour = neighbour[neighbourCell[cellj]];
                scalar neighboursNeighbourCellVolumes = 0.0;

                forAll (neighboursNeighbour, cellk)
                {
                    neighboursNeighbourCellVolumes +=
                        cellVolumes[neighboursNeighbour[cellk]];
                }

                neighboursEx +=
                    alpha1Ex[neighbourCell[cellj]]*cellVolumes[celli]
                   /neighboursNeighbourCellVolumes;
            }

            alpha1[celli] += neighboursEx - alpha1Ex[celli];
        }
