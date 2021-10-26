/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::solidBodyMotionFunctions::rotatingMotion

Description
    SolidBodyMotionFvMesh 6DoF motion function.

    The rotation is defined by an origin and axis of rotation and an angular
    speed.

SourceFiles
    rotatingMotion.C

\*---------------------------------------------------------------------------*/

#ifndef linearRampRotatingMotion_H
#define linearRampRotatingMotion_H

#include "solidBodyMotionFunction.H"
#include "fvMesh.H"
#include "volFields.H"
#include "kinematicMomentumTransportModel.H"
#include "primitiveFields.H"
#include "point.H"
#include "Function1.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{

/*---------------------------------------------------------------------------*\
                          Class rotatingMotion Declaration
\*---------------------------------------------------------------------------*/

class linearRampRotatingMotion
:
    public solidBodyMotionFunction
{

private:
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Model dictionary
        dictionary coeffs_;

        //- Boundary patches for force calculation
        const wordList patchSet_;

        //- Name of pressure field
        const word pName_;

        //- Reference density
        const scalar rhoRef_;

        //- Origin of the axis
        const vector origin_;

        //- Axis vector
        const vector axis_;

        //- Moment of inertia
        const scalar momentOfIntertia_;

        //- Opposing torque (Nm)
        const scalar opposingTorque_;

        //- Angular velocity (rad/sec)
        scalar omega_;

        //- Rotation angle (rad)
        scalar angle_;

        //- Angular acceleration (rad/sec/sec)
        scalar alpha_;

public:

    //- Runtime type information
    TypeName("linearRampRotatingMotion");


    // Constructors

        //- Construct from components
        linearRampRotatingMotion
        (
            const dictionary& SBMFCoeffs,
            const Time& runTime
        );

        //- Disallow default bitwise copy construction
        linearRampRotatingMotion(const linearRampRotatingMotion&);


        //- Construct and return a clone
        virtual autoPtr<solidBodyMotionFunction> clone() const
        {
            return autoPtr<solidBodyMotionFunction>
            (
                new linearRampRotatingMotion
                (
                    SBMFCoeffs_,
                    time_
                )
            );
        }


    //- Destructor
    virtual ~linearRampRotatingMotion();


    // Member Functions

        //- Return the solid-body motion transformation septernion
        virtual septernion transformation();

        //- Update properties from given dictionary
        virtual bool read(const dictionary& SBMFCoeffs);

        //- Update the angle
        void updateAngle();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearRampRotatingMotion&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidBodyMotionFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
