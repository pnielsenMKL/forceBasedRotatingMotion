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

#include "forceBasedRotatingMotion.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "surfaceMesh.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidBodyMotionFunctions
{
    defineTypeNameAndDebug(forceBasedRotatingMotion, 0);
    addToRunTimeSelectionTable
    (
        solidBodyMotionFunction,
        forceBasedRotatingMotion,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::forceBasedRotatingMotion
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
:
    solidBodyMotionFunction(SBMFCoeffs, runTime),
    mesh_(time_.db().parent().lookupObject<fvMesh>("region0")),
    coeffs_(SBMFCoeffs.subDict("forceBasedRotatingMotionCoeffs")),
    patchSet_(coeffs_.lookup("patches")),
    pName_(coeffs_.lookupOrDefault<word>("pName", "p")),
    rhoRef_(readScalar(coeffs_.lookup("rhoRef"))),
    origin_(coeffs_.lookup("origin")),
    axis_(coeffs_.lookup("axis")),
    momentOfIntertia_(readScalar(coeffs_.lookup("momentOfIntertia"))),
    opposingTorque_(readScalar(coeffs_.lookup("opposingTorque"))),
    omega_(coeffs_.lookupOrDefault<scalar>("omega", 0)),
    angle_(0),
    alpha_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::~forceBasedRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::updateAngle()
{
    // Typedefs
    typedef incompressible::momentumTransportModel icoTurbModel;

    // Lookup turbulence model
    // - Note: this is currently limited to incompressible turbulence models
    const incompressible::momentumTransportModel& turb =
      mesh_.lookupObject<icoTurbModel>(momentumTransportModel::typeName);

    // Initialize forces/moments
    vector pressureForce = vector(0,0,0);
    vector pressureMoment = vector(0,0,0);
    vector viscousForce = vector(0,0,0);
    vector viscousMoment = vector(0,0,0);
    vector netForce = vector(0,0,0);
    vector netMoment = vector(0,0,0);

    // Get fields for calculating forces/moments
    const volScalarField& p = mesh_.lookupObject<volScalarField>(pName_);
    const volSymmTensorField& devSigma = turb.devSigma();

    forAll(patchSet_, i)
    {
        // Get patch index
        label patchI = mesh_.boundaryMesh().findPatchID(patchSet_[i]);

        // Get displacement vectors for faces from origin
        vectorField Md
        (
            mesh_.C().boundaryField()[patchI] - origin_
        );

        // Get normal forces on faces
        vectorField Fp
        (
            rhoRef_*p.boundaryField()[patchI]*mesh_.boundary()[patchI].Sf()
        );

        // Get moment of pressure force on each face
        vectorField Mp
        (
            Md ^ Fp
        );

        // Get viscous forces on faces
        vectorField Fv
        (
            rhoRef_*devSigma.boundaryField()[patchI]&mesh_.boundary()[patchI].Sf()
        );

        // Get moment of viscous force on each face
        vectorField Mv
        (
            Md ^ Fv
        );

        // Sum the pressure force/moment
        pressureForce += sum(Fp);
        pressureMoment += sum(Mp);
        viscousForce += sum(Fv);
        viscousMoment += sum(Mv);
    }

    // Combine for parallel case
    Pstream::combineGather(pressureForce, plusEqOp<vector>());
    Pstream::combineGather(viscousForce, plusEqOp<vector>());
    Pstream::combineGather(pressureMoment, plusEqOp<vector>());
    Pstream::combineGather(viscousMoment, plusEqOp<vector>());
    Pstream::combineScatter(pressureForce);
    Pstream::combineScatter(viscousForce);
    Pstream::combineScatter(pressureMoment);
    Pstream::combineScatter(viscousMoment);

    netForce = pressureForce + viscousForce;
    netMoment = pressureMoment + viscousMoment;

    Info << "Pressure force " << pressureForce << endl;
    Info << "Pressure moment " << pressureMoment << endl;
    Info << "Viscous force " << viscousForce << endl;
    Info << "Viscous moment " << viscousMoment << endl;
    Info << "Net force " << netForce << endl;
    Info << "Net moment " << netMoment << endl;

    //scalar t = time_.value();

    // Get timestep data
    scalar dt = time_.deltaTValue();
    scalar dto = time_.deltaT0Value();

    // Update the motion
    omega_ += 0.5*dto*alpha_;
    angle_ += dt*omega_;
    scalar appliedMoment = netMoment & axis_;
    scalar opposingTorque = sign(appliedMoment)*opposingTorque_;
    alpha_ = (appliedMoment - opposingTorque)/momentOfIntertia_;
    omega_ += 0.5*dt*alpha_;

    Info << "omega " << omega_ << endl;
}

Foam::septernion
Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::transformation()
{
    this->updateAngle();
    // Rotation around axis
    //scalar angle = omega_->integrate(0, t);

    quaternion R(axis_, angle_);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    //DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    return true;
}


// ************************************************************************* //
