    Info<< "Reading field bubble1\n" << endl;
    volScalarField bubble1
    (
        IOobject
        (
            "bubble1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble2\n" << endl;
    volScalarField bubble2
    (
        IOobject
        (
            "bubble2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble3\n" << endl;
    volScalarField bubble3
    (
        IOobject
        (
            "bubble3",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble4\n" << endl;
    volScalarField bubble4
    (
        IOobject
        (
            "bubble4",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble5\n" << endl;
    volScalarField bubble5
    (
        IOobject
        (
            "bubble5",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble6\n" << endl;
    volScalarField bubble6
    (
        IOobject
        (
            "bubble6",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble7\n" << endl;
    volScalarField bubble7
    (
        IOobject
        (
            "bubble7",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field bubble8\n" << endl;
    volScalarField bubble8
    (
        IOobject
        (
            "bubble8",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

surfaceScalarField CSF
(
    IOobject
    (
        "CSF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mixture.surfaceTensionForce()
);

volVectorField f_sigma
(
    IOobject
    (
        "f_sigma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::reconstruct(CSF * mesh.magSf())
);

mesh.setFluxRequired(bubble1.name());
mesh.setFluxRequired(bubble2.name());
mesh.setFluxRequired(bubble3.name());
mesh.setFluxRequired(bubble4.name());
mesh.setFluxRequired(bubble5.name());
mesh.setFluxRequired(bubble6.name());
mesh.setFluxRequired(bubble7.name());
mesh.setFluxRequired(bubble8.name());
