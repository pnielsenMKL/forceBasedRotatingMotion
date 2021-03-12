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
    coeffs_(SBMFCoeffs.subDict("forceBasedRotatingMotionCoeffs")),
    patchSet_(coeffs_.lookup("patches")),
    pName_(coeffs_.lookup("pName")),
    rhoRef_(readScalar(coeffs_.lookup("rhoRef"))),
    origin_(coeffs_.lookup("origin")),
    axis_(coeffs_.lookup("axis")),
    omega_(Function1<scalar>::New("omega", coeffs_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::~forceBasedRotatingMotion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::septernion
Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::transformation() const
{
    // Testing forces
    vector force = vector(0,0,0);
    const fvMesh& mesh = time_.db().parent().lookupObject<fvMesh>("region0"); // Needs to be generalized and moved to a member variable
    const volScalarField& p = mesh.lookupObject<volScalarField>(pName_);
    const volVectorField& U = mesh.lookupObject<volVectorField>("U");

    forAll(patchSet_, i)
    {
        // Get patch index
        label patchI = mesh.boundaryMesh().findPatchID(patchSet_[i]);

        // Get displacement vectors for faces from origin
        vectorField Md
        (
            mesh.C().boundaryField()[patchI] - origin_
        );

        // Get normal forces on faces
        scalar pRef = 0; // Needs to be read from dict
        vectorField fN
        (
            rhoRef_*mesh.boundary()[patchI].Sf()*(p.boundaryField()[patchI] - pRef)
        );

	// Sum the pressure force
        force = sum(fN);
	Info << "Pressure force " << force << endl;
    }

    scalar t = time_.value();

    // Rotation around axis
    scalar angle = omega_->integrate(0, t);

    quaternion R(axis_, angle);
    septernion TR(septernion(-origin_)*R*septernion(origin_));

    DebugInFunction << "Time = " << t << " transformation: " << TR << endl;

    return TR;
}


bool Foam::solidBodyMotionFunctions::forceBasedRotatingMotion::read
(
    const dictionary& SBMFCoeffs
)
{
    solidBodyMotionFunction::read(SBMFCoeffs);

    omega_.reset
    (
        Function1<scalar>::New("omega", SBMFCoeffs_).ptr()
    );

    return true;
}


// ************************************************************************* //
