#pragma once

#include "elUtil.h"

#define X mass
#define X_all mass_all
#define rest_Psi rest_density
#define rest_Psi_all rest_density_all

//#define X rest_volume
//#define X_all rest_volume_all
//#define rest_Psi rest_compressionRate
//#define rest_Psi_all rest_compressionRate_all


// Grid config mode

constexpr val_i basicGridConfig[3] = { 0, 3, 0 }; // coressponding to 3 blocks, each of which has 1, 3, 3 elements respectviely
constexpr val_i objectGridConfig[3] = { 3, 3, 3 };
constexpr val_i fluidGridConfig[4] = { 3, 3, 45 };
constexpr val_i boundGridConfig[3] = { 3, 3, 24 };
constexpr val_i phaseGridConfig[3] = { 0, 0, 24 };
constexpr val_i neighbAuxConfig[3] = { 6, 0, 0 };

// Simulation domain value

class Elvari {
public:
    long int frameNumber = 0;
    Array<val_f, 3, 1> gravity;
    Array<val_f, 3, 1> spaceMin;
    Array<val_f, 3, 1> spaceMax;
    val_f currentTime;
    val_i phaseNumber;

    Elvari() {
        gravity = Array < val_f, 3, 1>(0, -9.8, 0);
        spaceMin = Array < val_f, 3, 1>(-8, -8, -8);
        spaceMax = Array < val_f, 3, 1>(8, 12, 8);
        phaseNumber = 1;
        currentTime = 0;
    }
};

constexpr val_i MaxPartNum_fluid = 2e6;        // constraint of max particle num (fluid)
constexpr val_i MaxPartNum_bound = 2e5;        // constraint of max particle num (boundary)
constexpr val_i MaxPartNum = MaxPartNum_fluid + MaxPartNum_bound;

constexpr val_f diamPart = 0.05;              // particle cuboid length
constexpr val_f diamPart2 = diamPart * diamPart;
constexpr val_f diamPart3 = diamPart2 * diamPart;
constexpr val_f diamPart4 = diamPart3 * diamPart;
constexpr val_f smoothRadius = diamPart * 2; // smooth radius
constexpr val_f smoothRadius2 = smoothRadius * smoothRadius;
constexpr val_f smoothRadius3 = smoothRadius2 * smoothRadius;
constexpr val_f smoothRadius4 = smoothRadius3 * smoothRadius;
constexpr val_f sig3 = M_1_PI / diamPart3;
constexpr val_f sig3_grad = M_1_PI / diamPart4;
constexpr val_f restVolume = diamPart3;


constexpr val_i speedSound = 100;            // speed of sound (simulation time step related)
constexpr val_f timeStep = val_f(diamPart / speedSound); // fixed time step
constexpr val_f timeStep2 = val_f(diamPart * diamPart / speedSound / speedSound);
constexpr val_f timeStep_1 = 1 / timeStep;
constexpr val_f timeStep2_1 = 1 / timeStep2;

constexpr val_f viscosity_nu = 0.1;
constexpr val_i maxIncompressibleSolverIter = 50;
constexpr val_i maxDivergenceFreeSolverIter = 30;
constexpr val_i minIncompressibleSolverIter = 2;
constexpr val_i minDivergenceFreeSolverIter = 1;
constexpr val_f incompressibleThreshold = 1e-4;
constexpr val_f divergenceFreeThreshold = 1e-3; 
constexpr val_f incompressibleThreshold_ii = 1e-4;
constexpr val_f relaxingFactor = 0.5;
constexpr val_f gamma_pb = 0.7;

constexpr val_i gridMembers = 70;           // element contains in one neighbour search grid node
constexpr val_i nonZeros = gridMembers - 1;
constexpr val_i colLock = gridMembers - 2;
constexpr val_i scanRange[2] = { -1,2 };
constexpr val_f neighbSearchGridSize = smoothRadius; // size of the neighbour search grid

// Multiphase
#define Multiphase // undefine this to unable volume transaction and diffusion
constexpr val_i phaseNumber = 3; // <=5
constexpr val_f transactionCoefficient = 0.1;
constexpr val_f diffusionCoefficinet = 0;


// SPH kernel computation

inline val_f W(val_f r) {
    return cubicSplineKernel3d(r, diamPart, sig3);
}
inline val_f grad_W(val_f r) {
    return gradCubicSplineKernel3d(r, diamPart, sig3_grad);
}

inline Array3_f laplacian_W(Array3_f& Aij, Array3_f& xij, val_f dis, Array3_f& grad_W, val_f Vj) {
    return 10 * Vj * (Aij.matrix().transpose() * (xij.matrix()))(0, 0) / (pow(dis, 2) + 10e-6) * grad_W;
}

inline val_i bid(val_i uid) {
    return uid - MaxPartNum_fluid;
}

inline bool isFluid(val_i uid) {
    if (uid < MaxPartNum_fluid) {
        return true;
    }
    else {
        return false;
    }
}