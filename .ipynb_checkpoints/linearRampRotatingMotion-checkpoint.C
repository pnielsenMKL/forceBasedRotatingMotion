/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "linearRampRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "surfaceMesh.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(linearRampRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        linearRampRotatingMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearRampRotatingMotion::linearRampRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    mesh_(time_.db().parent().lookupObject<fvMesh>("region0")),
    coeffs_(SBMFCoeffs.subDict("linearRampRotatingMotionCoeffs")),
    origin_(coeffs_.lookup("origin")),
    axis_(coeffs_.lookup("axis")),
    omegaIn_(coeffs_.lookupOrDefault<scalar>("omegaInitial", 0)),
    omegaFin_(coeffs_.lookupOrDefault<scalar>("omegaFinal", 0)),
    thetaIn_(coeffs_.lookupOrDefault<scalar>("thetaInitial", 0)),
    thetaFin_(coeffs_.lookupOrDefault<scalar>("thetaFinal", 0)),
//    angle_(0),
    angle_(coeffs_.lookupOrDefault<scalar>("angle", 0))//added 
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::linearRampRotatingMotion::~linearRampRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::solidBodyMotionFunctions::linearRampRotatingMotion::updateAngle()
{
    //scalar t = time_.value();

    // Get timestep data
    scalar dt = time_.deltaTValue();

    // Update the motion
    angle_ += dt*omega_;
    
    if (angle_ < thetaFin_)
    {
        omega_ = omegaIn_ - (omegaIn_ - omegaFin_)*(angle_/thetaFin_);        
    }
    else
    {
        omega_ = omegaFin_;
    }

    Info << "omega " << omega_ << endl;
    Info << "angle " << angle_ << endl;
}

Foam::septernion
Foam::solidBodyMotionFunctions::linearRampRotatingMotion::transformation()
{
    this->updateAngle();
    // Rotation around axis
    //scalar angle = omega_->integrate(0, t);

    quaternion R(axis_, angle_);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    //DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::linearRampRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
