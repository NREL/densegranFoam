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
    Foam::nlgf

Description

SourceFiles
    nlgfPressure.C

\*---------------------------------------------------------------------------*/

#ifndef nlgf_H
#define nlgf_H

#include "fcicGranularPressureModel.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace biomassModels
{
namespace fcicGranularPressureModels
{

/*---------------------------------------------------------------------------*\
                           Class nlgf Declaration
\*---------------------------------------------------------------------------*/

class nlgf
:
    public fcicGranularPressureModel
{

  dictionary coeffDict_;
  volScalarField pres_;
  uniformDimensionedVectorField g_;
  scalar pres_offset_;
  volScalarField href_nlgf;
  volScalarField p_rgh;
  //  volScalarField pres_reg;
  const uniformDimensionedScalarField& href_;
  scalar alphaBoundary;
  scalar nNormalCell;
  scalar theta_;
public:

    //- Runtime type information
    TypeName("nlgf");


    // Constructors

        //- Construct from components
  nlgf(
       const dictionary& dict,
       const volVectorField& U,
       const volScalarField& alpha,
       const volScalarField& rho
       );


    //- Destructor
    virtual ~nlgf();


    // Member Functions


        tmp<volScalarField> granularPressureCoeff
        () const;

  virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace granularPressureModels
} // End namespace biomassModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
