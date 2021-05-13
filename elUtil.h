#pragma once

#include <iostream>
#include <iomanip>
#include <math.h>
#include <omp.h>
#include <time.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"

using namespace std;
using namespace Eigen;

#define M_PI        3.14159265358979323846264338327950288   // pi  
#define M_1_PI      0.318309886183790671537767526745028724  // 1/pi 

static omp_lock_t lock;

// Universal parameter type 
// val_i
typedef int val_i;
typedef double val_f;

typedef Array<val_i, 3, 1> Array3_i;
typedef Array<val_f, 3, 1> Array3_f;
typedef Block<Array<val_i, Dynamic, Dynamic>> IntegerBlockRef;
typedef Block<Array<val_f, Dynamic, Dynamic>> FloatBlockRef;

// Marker definition

constexpr val_i POS_X = 0;
constexpr val_i POS_Y = 1;
constexpr val_i POS_Z = 2;

constexpr val_i FLUID = 0;
constexpr val_i BOUND = 1;

constexpr val_i ROW = 0;
constexpr val_i COL = 1;

constexpr val_i ELE = 0;
constexpr val_i NUM = 1;

// Timing

#define CALL_TIME(a, name) \
    double et1, et2;\
    et1 = omp_get_wtime();  \
    a;                     \
    et2 = omp_get_wtime();  \
    cout << "time execute " << name << " is: " << (et2 - et1) * 1000 << "ms" << endl;

#define PAUSE cout << "input a character & press enter to end"; char atemp; cin >> atemp;

#define cout_size(a, name) cout << "sizeof(" << name << "): " << sizeof(a) << " (float is 4)" << endl;

#define toArray3f(a) Map<Array<val_f,3,1>>(a)

// Matrix level calculation

inline val_i iterProduct(Array3_i& arr)
{
    val_i tmp = 1;
    for (val_i i = 0; i < arr.rows(); i++) {
        for (val_i j = 0; j < arr.cols(); j++) {
            tmp *= arr(i, j);
        }
    }
    return tmp;
}

// Matrix level operation

template <class T, val_i I, val_i J>
inline void expandMat(Array<T, I, J>& arr, val_i rc, val_i num, val_f value = 0)
{
    if (rc == ROW) {
        arr.conservativeResize(arr.rows() + num, arr.cols());
        arr.block(arr.rows() - num, 0, num, arr.cols()).setConstant(value);
    }
    else if (rc == COL) {
        arr.conservativeResize(arr.rows(), arr.cols() + num);
        arr.block(0, arr.cols() - num, arr.rows(), num).setConstant(value);
    }
}

template <class T, val_i I, val_i J>
inline void contractMat(Array<T, I, J>& arr, val_i rc, val_i num)
{
    if (rc == ROW) {
        arr.conservativeResize(arr.rows() - num, arr.cols());
    }
    else if (rc == COL) {
        arr.conservativeResize(arr.rows(), arr.cols() - num);
    }
}

template <class T, val_i I, val_i J>
inline void appendElement(Array<T, I, J>& arr, val_f value = 0)
{
    expandMat(arr, ROW, 1, value);
}

template <class T, val_i I, val_i J>
inline void pullElement(Array<T, I, J>& arr, val_i num = 1)
{
    contractMat(arr, ROW, num);
}

// SPH calculation

inline val_f cubicSplineKernel3d(val_f dis, val_f diam, val_f sig)
// NOT intend to be used directly, just as a sample. 
// New one is in globalValue.h for better performance
{
    val_f q = dis / diam;
    val_f tmp = 0;
    if (q < 1) { tmp = 1 - 1.5 * pow(q, 2) + 0.75 * pow(q, 3); }
    else if (q > 1 && q < 2) { tmp = 0.25 * pow(2 - q, 3); }
    else if (q == 1) { tmp = 0.25; }
    else { return 0; }
    tmp *= sig;
    return tmp;
}

inline val_f gradCubicSplineKernel3d(val_f dis, val_f diam, val_f sig_drad)
// NOT intend to be used directly, just as a sample. 
// New one is in globalValue.h for better performance
{
    val_f q = dis / diam;
    val_f tmp;
    if (q < 1) { tmp = 2.25 * pow(q, 2) - 3 * q; }
    else if (q > 1 && q < 2) { tmp = -0.75 * pow(2 - q, 2); }
    else if (q == 1) { tmp = -0.75; }
    else { return 0; }
    tmp *= sig_drad;
    return tmp;
}

// NOT intend to be used directly, just as a sample. 
// New one is in globalValue.h for better performance
inline Array3_f artificialLaplacian3d(FloatBlockRef xij, val_f xijDis, val_f diam, val_f sig_drad, val_f Vj, FloatBlockRef Aij)
{
    return 10 * Vj * (Aij.matrix().transpose() * (xij.matrix())) / (pow(xijDis,2)+10e-6) * gradCubicSplineKernel3d(xijDis, diam, sig_drad);
}