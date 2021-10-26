/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::muJphiJModels::effectivePressureModels::JohnsonJacksonBoyer

Description

SourceFiles
    JohnsonJacksonBoyer.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonJacksonBoyer_H
#define JohnsonJacksonBoyer_H

#include "effectivePressureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace muJphiJModels
{
namespace effectivePressureModels
{

/*---------------------------------------------------------------------------*\
               Class JohnsonJacksonBoyer Declaration
\*---------------------------------------------------------------------------*/

class JohnsonJacksonBoyer
:
    public effectivePressureModel
{
    // Private data

        dictionary coeffDict_;

        //- Material constant for frictional normal stress
  dimensionedScalar aMJ_;
  dimensionedScalar phiRLP_;
  dimensionedScalar phiRCP_;
  dimensionedScalar alphaDeltaMin_;
  
  //  volScalarField phiPM_;
  

public:

    //- Runtime type information
    TypeName("JohnsonJacksonBoyer");


    // Constructors

        //- Construct from components
        JohnsonJacksonBoyer(const dictionary& dict);


    //- Destructor
    virtual ~JohnsonJacksonBoyer();


    // Member functions

        virtual tmp<volScalarField> effectivePressure
        (
	 const phaseModel& phase,
	 const dimensionedScalar& phiM_,
	 const dimensionedScalar& nuC_,
	 const dimensionedScalar& rhoC_,
	 const dimensionedScalar& S0_
        ) const;

        virtual tmp<volScalarField> effectivePressurePrime
        (
	 const phaseModel& phase,
	 const dimensionedScalar& phiM_,
	 const dimensionedScalar& nuC_,
	 const dimensionedScalar& rhoC_,
	 const dimensionedScalar& S0_
        ) const;

        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace effectivePressureModels
} // End namespace muJphiJModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
