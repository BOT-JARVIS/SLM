#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


dimensionedScalar Tl2("Tl2", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);
dimensionedScalar Tl3("Tl3", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);

dimensionedScalar Ts2("Ts2", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);
dimensionedScalar Ts3("Ts3", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);

dimensionedScalar cp1("cp1", dimensionSet(0, 2, -2, -1, 0, 0, 0),transportProperties);
dimensionedScalar cp2s("cp2s", dimensionSet(0, 2, -2, -1, 0, 0, 0),transportProperties);
dimensionedScalar cp2l("cp2l", dimensionSet(0, 2, -2, -1, 0, 0, 0),transportProperties);
dimensionedScalar cp3s("cp3s", dimensionSet(0, 2, -2, -1, 0, 0, 0),transportProperties);
dimensionedScalar cp3l("cp3l", dimensionSet(0, 2, -2, -1, 0, 0, 0),transportProperties);

dimensionedScalar k1("k1", dimensionSet(1, 1, -3, -1, 0, 0, 0),transportProperties);
dimensionedScalar k2s("k2s", dimensionSet(1, 1, -3, -1, 0, 0, 0),transportProperties);
dimensionedScalar k2l("k2l", dimensionSet(1, 1, -3, -1, 0, 0, 0),transportProperties);
dimensionedScalar k3s("k3s", dimensionSet(1, 1, -3, -1, 0, 0, 0),transportProperties);
dimensionedScalar k3l("k3l", dimensionSet(1, 1, -3, -1, 0, 0, 0),transportProperties);

dimensionedScalar beta1("beta1", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);
dimensionedScalar beta2("beta2", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);
dimensionedScalar beta3("beta3", dimensionSet(0, 0, 0, -1, 0, 0, 0), transportProperties);

dimensionedScalar L2("L2", dimensionSet(0, 2, -2, 0, 0, 0, 0), transportProperties);
dimensionedScalar L3("L3", dimensionSet(0, 2, -2, 0, 0, 0, 0), transportProperties);

dimensionedScalar dsigmadT2("dsigmadT2", dimensionSet(1, 0, -2, -1, 0, 0, 0),transportProperties);
dimensionedScalar dsigmadT3("dsigmadT3", dimensionSet(1, 0, -2, -1, 0, 0, 0),transportProperties);

dimensionedScalar L_v2("L_v2", dimensionSet(0, 2, -2, 0, 0, 0, 0),transportProperties);
dimensionedScalar M2("M2", dimensionSet(1, 0, 0, 0, -1, 0, 0),transportProperties);
dimensionedScalar T_v2("T_v2", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);
dimensionedScalar L_v3("L_v3", dimensionSet(0, 2, -2, 0, 0, 0, 0),transportProperties);
dimensionedScalar M3("M3", dimensionSet(1, 0, 0, 0, -1, 0, 0),transportProperties);
dimensionedScalar T_v3("T_v3", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);
dimensionedScalar R("R", dimensionSet(1, 2, -2, -1, -1, 0, 0),transportProperties);
dimensionedScalar P_0("P_0", dimensionSet(1, -1, -2, 0, 0, 0, 0),transportProperties);


dimensionedScalar C("C", dimensionSet(1, -3, -1, 0, 0, 0, 0),transportProperties);
dimensionedScalar b("b", dimensionSet(0, 0, 0, 0, 0, 0, 0),transportProperties);

dimensionedScalar E2("E2", dimensionSet(0, 0, 0, 0, 0, 0, 0),transportProperties);
dimensionedScalar E3("E3", dimensionSet(0, 0, 0, 0, 0, 0, 0),transportProperties);
dimensionedScalar sig("sig", dimensionSet(1, 0, -3, -4, 0, 0, 0),transportProperties);
dimensionedScalar Tref("Tref", dimensionSet(0, 0, 0, 1, 0, 0, 0),transportProperties);
dimensionedScalar h("h", dimensionSet(1, 0, -3, -1, 0, 0, 0),transportProperties);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "createPhi.H"

immiscibleIncompressibleThreePhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());
volScalarField& alpha3(mixture.alpha3());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();
const dimensionedScalar& rho3 = mixture.rho3();

dimensionedScalar D23("D23", dimViscosity, mixture);

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2 + alpha3*rho3
);
rho.oldTime();


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho1*phi
);


// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha2.name());

#include "createMRF.H"
#include "createFvOptions.H"


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Tl
(
    IOobject
    (
        "Tl",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Ts
(
    IOobject
    (
        "Ts",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField lf
(
    IOobject
    (
        "lf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField cp2 = cp2s*(1-lf) + cp2l*lf;
volScalarField cp3 = cp3s*(1-lf) + cp3l*lf;

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*rho1*cp1 +  alpha2*rho2*cp2 + alpha3*rho3*cp3
);
rhoCp.oldTime();

surfaceScalarField rhoCpphi
(
    IOobject
    (
        "rhoCpphi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rhoCp)*phi
);


volScalarField k2 = k2s*(1-lf) + k2l*lf;
volScalarField k3 = k3s*(1-lf) + k3l*lf;
volScalarField k = alpha1*k1 +  alpha2*k2 + alpha3*k3;
surfaceScalarField kappaf =  fvc::interpolate(k);
volScalarField beta = alpha1*beta1 +  alpha2*beta2 +  alpha3*beta3;

volScalarField A =   -C*(1-alpha1)*sqr(scalar(1)-lf)/ (pow(lf,scalar(3))+b);
volScalarField factor1 = (scalar(2.0)*rho)/(rho1+rho2+rho3);
volScalarField factor2 = (scalar(2.0)*rhoCp)/(rho1*cp1 + rho2*cp2);
volVectorField gradAlpha = fvc::grad(alpha1);
volVectorField nHat = gradAlpha/(mag(gradAlpha) + mixture.deltaN());
volScalarField kappa = mixture.K();
volScalarField B = lf;

volScalarField vaploss2 = 0.82 * ((L_v2 * M2) / (sqrt(2 * 3.14 * M2 * R * T))) * P_0 * exp(( L_v2 * M2 * (T - T_v2))/( R * T * T_v2));
volScalarField vaploss3 = 0.82 * ((L_v3 * M3) / (sqrt(2 * 3.14 * M3 * R * T))) * P_0 * exp(( L_v3 * M3 * (T - T_v3))/( R * T * T_v3));
volScalarField vaploss = alpha2*vaploss2 + alpha3*vaploss3;
volScalarField radloss2= E2*sig*(pow(T,scalar(4.0)) - pow(Tref,scalar(4.0)));
volScalarField radloss3= E3*sig*(pow(T,scalar(4.0)) - pow(Tref,scalar(4.0)));
volScalarField radloss=  alpha2*radloss2 + alpha3*radloss3;
volScalarField convloss = h * (T - Tref);



///////////////////////////////////////////////////// Radiation model  //////////////////////////////////////////////////////////////////////////////////

volScalarField alphaT = alpha1;
volScalarField s1 = alpha1;
// Creating radiation model
autoPtr<radiation::radiationModel> radiation
(
        radiation::radiationModel::New(alphaT)
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////










