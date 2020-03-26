//  READING PHYSICAL PROPERTIES

    Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar epsilon
    (
        physicalProperties.lookup("epsilon0")
    );

    dimensionedScalar e
    (
        physicalProperties.lookup("e")
    );

    dimensionedScalar nRef
    (
        physicalProperties.lookup("nRef")
    );

    dimensionedScalar muPlus
    (
        physicalProperties.lookup("muPlus")
    );

    dimensionedScalar muMinus
    (
        physicalProperties.lookup("muMinus")
    );

    dimensionedScalar DPlus
    (
        physicalProperties.lookup("DPlus")
    );

    dimensionedScalar DMinus
    (
        physicalProperties.lookup("DMinus")
    );

    dimensionedScalar ZPlus
    (
        physicalProperties.lookup("ZPlus")
    );

    dimensionedScalar ZMinus
    (
        physicalProperties.lookup("ZMinus")
    );

    dimensionedScalar nMinimal
    (
        physicalProperties.lookup("nMinimal")
    );


//  READING FIELDS

    Info<< "Reading field ePhi\n" << endl;
    volScalarField ePhi
    (
        IOobject
        (
            "ePhi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field PhiGrad\n" << endl;
    volVectorField PhiGrad
    (
        IOobject
        (
            "PhiGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(ePhi)
    );

    Info<< "Creating field snPhiGrad\n" << endl;
    surfaceScalarField snPhiGrad
    (
        IOobject
        (
            "snPhiGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::snGrad(ePhi) * mesh.magSf()
    );


    Info<< "Reading field nPlus\n" << endl;
    volScalarField nPlus
    (
        IOobject
        (
            "nPlus",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field nMinus\n" << endl;

    volScalarField nMinus
    (
        IOobject
        (
            "nMinus",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Creating field netCharge\n" << endl;

    volScalarField netCharge
    (
        IOobject
        (
            "netCharge",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        nRef*e*(nPlus*ZPlus-nMinus*ZMinus)
    );

    Info<< "Creating fields muSnPhiGrad\n" << endl;

    surfaceScalarField muPlusSnPhiGrad
    (
        IOobject
        (
            "muPlusSnPhiGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        snPhiGrad*muPlus
    );

    surfaceScalarField muMinusSnPhiGrad
    (
        IOobject
        (
            "muMinusSnPhiGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        snPhiGrad*muMinus
    );

    volVectorField nPlusGrad
    (
        IOobject
        (
            "nPlusGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(nPlus)
    );

    volVectorField nMinusGrad
    (
        IOobject
        (
            "nMinusGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::grad(nMinus)
    );

    surfaceScalarField nPlusSnGrad
    (
        IOobject
        (
            "nPlusSnGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(nPlusGrad) & mesh.Sf()
    );

    surfaceScalarField nMinusSnGrad
    (
        IOobject
        (
            "nMinusSnGrad",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(nMinusGrad) & mesh.Sf()
    );

    volVectorField nPlusCurrent
    (
        IOobject
        (
            "nPlusCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (- nPlusGrad*DPlus - nPlus*muPlus*PhiGrad)
    );

    volVectorField nMinusCurrent
    (
        IOobject
        (
            "nMinusCurrent",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (- nMinusGrad*DMinus + nMinus*muMinus*PhiGrad)
    );

    Info << "Create fields nSnVelocity\n" << endl;

    surfaceScalarField nPlusSnVelocity
    (
        IOobject
        (
            "nPlusSnVelocity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (- nPlusSnGrad*DPlus/max(fvc::interpolate(nPlus),nMinimal) - muPlusSnPhiGrad)
    );

    surfaceScalarField nMinusSnVelocity
    (
        IOobject
        (
            "nMinusSnVelocity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (- nMinusSnGrad*DMinus/max(fvc::interpolate(nMinus),nMinimal) + muMinusSnPhiGrad)
    );

    Info << "Create fields nCourantNo\n" << endl;

    volScalarField nPlusCourantNo
    (
        IOobject
        (
            "nPlusCourantNo",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    volScalarField nMinusCourantNo
    (
        IOobject
        (
            "nMinusCourantNo",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
 	mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    Info << "Create field nMinimals\n" << endl;

    surfaceScalarField nMinimals
    (
        IOobject
        (
            "nMinimals",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
 	mesh,
        dimensionedScalar(nMinimal)
    );

    const label nOuterCorrIons = mesh.solutionDict().subDict("PISO").lookupOrDefault("nOuterCorrIons", 1);
    const label nCorrIons = mesh.solutionDict().subDict("PISO").lookupOrDefault("nCorrIons", 1);

