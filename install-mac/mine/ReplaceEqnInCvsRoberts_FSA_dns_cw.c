/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2021, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Example problem:
 *
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODES for Forward Sensitivity 
 * Analysis. The problem is from chemical kinetics, and consists
 * of the following three rate equations:
 *    dy1/dt = -p1*y1 + p2*y2*y3
 *    dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
 *    dy3/dt =  p3*(y2)^2
 * on the interval from t = 0.0 to t = 4.e10, with initial
 * conditions y1 = 1.0, y2 = y3 = 0. The reaction rates are: p1=0.04,
 * p2=1e4, and p3=3e7. The problem is stiff.
 * This program solves the problem with the BDF method, Newton
 * iteration with the DENSE linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance.
 * Output is printed in decades from t = .4 to t = 4.e10.
 * Run statistics (optional outputs) are printed at the end.
 *
 * Optionally, CVODES can compute sensitivities with respect to the
 * problem parameters p1, p2, and p3.
 * The sensitivity right hand side is given analytically through the
 * user routine fS (of type SensRhs1Fn).
 * Any of three sensitivity methods (SIMULTANEOUS, STAGGERED, and
 * STAGGERED1) can be used and sensitivities may be included in the
 * error test or not (error control set on SUNTRUE or SUNFALSE,
 * respectively).
 *
 * Execution:
 *
 * If no sensitivities are desired:
 *    % cvsRoberts_FSA_dns -nosensi
 * If sensitivities are to be computed:
 *    % cvsRoberts_FSA_dns -sensi sensi_meth err_con
 * where sensi_meth is one of {sim, stg, stg1} and err_con is one of
 * {t, f}.
 * -----------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cvodes/cvodes.h>             /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>    /* definition of ABS */
#include "helpers.h"

/* Accessor macros */

#define Ith(v, i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A, i, j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Problem Constants */
#define NEQ 4 // num equations.
#define NP 10 // num parameters in the model
#define NS 1  // Number of parameters we want sensitivities for

#define ATOL  RCONST(1.e-5) /* scalar absolute tolerance */
#define T0    RCONST(0.0)   /* initial time              */
#define T1    RCONST(0.5)   /* first output time         */
#define DTOUT RCONST(0.5)   /* output time increment     */
#define NOUT  10            /* number of output times    */
#define TMULT 10            /* output time factor */

#define RTOL  RCONST(1e-4)  /* scalar relative tolerance */
#define ATOL1 RCONST(1e-8)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1e-14)
#define ATOL3 RCONST(1e-6)
#define ATOL4 RCONST(1e-6)

#define keff1    0.0017
#define keff2    1
#define keff3    0.03
#define n        3
#define u1        0.0001
#define u3        0.0001
#define u2        0.001
#define u4        0.001
#define alpha1    9e-5
#define alpha2    0.001


#define scUPA0    1.16213e-8
#define PLG0    0.0262792
#define PLS0    19.7847
#define tcUPA0    19.9809


/* Type : UserData */

typedef struct {
    realtype p[NEQ];           /* problem parameters */
} *UserData;

/* Prototypes of functions by CVODES */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

/* Prototypes of private functions */

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int main(int argc, char *argv[]) {
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    UserData data;
    realtype t, tout;
    N_Vector y;
    int iout, retval;

    realtype pbar[NS];
    int is;
    N_Vector *yS;
    booleantype sensi, err_con;
    int sensi_meth;

    cvode_mem = NULL;
    data = NULL;
    y = NULL;
    yS = NULL;
    A = NULL;
    LS = NULL;

    sensi = 1;

    /* Process arguments */
//  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);
    sensi_meth = CV_SIMULTANEOUS;
//    sensi_meth = CV_STAGGERED;
    err_con = 1;
//    err_con = 0;

    /* User data structure */
    data = (UserData) malloc(sizeof *data);
    if (check_retval((void *) data, "malloc", 2)) return (1);
    data->p[0] = RCONST(0.04);
    data->p[1] = RCONST(1.0e4);
    data->p[2] = RCONST(3.0e7);

    /* Initial conditions */
    y = N_VNew_Serial(NEQ);
    if (check_retval((void *) y, "N_VNew_Serial", 0)) return (1);

    double* yData = y->ops->nvgetarraypointer(y);

    yData[0] = scUPA0;
    yData[1] = PLG0;
    yData[2] = PLS0;
    yData[3] = tcUPA0;

    /* Create CVODES object */
//    cvode_mem = CVodeCreate(CV_ADAMS);
    cvode_mem = CVodeCreate(CV_BDF);
    if (check_retval((void *) cvode_mem, "CVodeCreate", 0)) return (1);

    /* Allocate space for CVODES */
    retval = CVodeInit(cvode_mem, f, T0, y);
    if (check_retval(&retval, "CVodeInit", 1)) return (1);

    /* Use private function to compute error weights */
//  retval = CVodeWFtolerances(cvode_mem, ewt);
//  if (check_retval(&retval, "CVodeSetEwtFn", 1)) return(1);
    retval = CVodeSStolerances(cvode_mem, 1e-6, 1e-12);
    if (check_retval(&retval, "CVodeSStolerances", 1)) return (1);

    /* Attach user data */
    retval = CVodeSetUserData(cvode_mem, data);
    if (check_retval(&retval, "CVodeSetUserData", 1)) return (1);

    /* Create dense SUNMatrix */
    A = SUNDenseMatrix(NEQ, NEQ);
    if (check_retval((void *) A, "SUNDenseMatrix", 0)) return (1);

    /* Create dense SUNLinearSolver */
    LS = SUNLinSol_Dense(y, A);
    if (check_retval((void *) LS, "SUNLinSol_Dense", 0)) return (1);

    /* Attach the matrix and linear solver */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if (check_retval(&retval, "CVodeSetLinearSolver", 1)) return (1);

    /* Set the user-supplied Jacobian routine Jac */
//  retval = CVodeSetJacFn(cvode_mem, Jac);
    retval = CVodeSetJacFn(cvode_mem, NULL);
    if (check_retval(&retval, "CVodeSetJacFn", 1)) return (1);

    retval = CVodeSetMaxNumSteps(cvode_mem, (long)1e8);
    if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return (1);

    printf("\n3-species chemical kinetics problem\n");

    /* Sensitivity-related settings */
    if (sensi) {

        /* Set parameter scaling factor */
        pbar[0] = data->p[0];
//        pbar[1] = data->p[1];
//        pbar[2] = data->p[2];

        /* Set sensitivity initial conditions */
        yS = N_VCloneVectorArray(NS, y);
        if (check_retval((void *) yS, "N_VCloneVectorArray", 0)) return (1);
//        for (is = 0; is < NS; is++) N_VConst(ZERO, yS[is]);
        setSensInitTo(yS, NS, NEQ, 0.0);


        /* Call CVodeSensInit1 to activate forward sensitivity computations
           and allocate internal memory for COVEDS related to sensitivity
           calculations. Computes the right-hand sides of the sensitivity
           ODE, one at a time */

        retval = CVodeSensInit(cvode_mem, NS, sensi_meth, NULL, yS);
        if (check_retval(&retval, "CVodeSensInit", 1)) return (1);

        /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
           variables based on the rolerances supplied for states variables and
           the scaling factor pbar */
//        double reltol = 1e-4;
//        double abstol = 1e-6;
//        retval = CVodeSensSStolerances(cvode_mem, reltol, &abstol);
        retval = CVodeSensEEtolerances(cvode_mem);
        if (check_retval(&retval, "CVodeSensEEtolerances", 1)) return (1);

        /* Set sensitivity analysis optional inputs */
        /* Call CVodeSetSensErrCon to specify the error control strategy for
           sensitivity variables */
        retval = CVodeSetSensErrCon(cvode_mem, err_con);
        if (check_retval(&retval, "CVodeSetSensErrCon", 1)) return (1);

        /* Call CVodeSetSensParams to specify problem parameter information for
           sensitivity calculations */
        retval = CVodeSetSensParams(cvode_mem, data->p, pbar, NULL);
        if (check_retval(&retval, "CVodeSetSensParams", 1)) return (1);

        printf("Sensitivity: YES ");
        if (sensi_meth == CV_SIMULTANEOUS)
            printf("( SIMULTANEOUS +");
        else if (sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
        else printf("( STAGGERED1 +");
        if (err_con) printf(" FULL ERROR CONTROL )");
        else printf(" PARTIAL ERROR CONTROL )");

    } else {

        printf("Sensitivity: NO ");

    }

    /* In loop over output points, call CVode, print results, test for error */

    printf("\n\n");

    for (iout = 1, tout = T1; iout <= NOUT; iout++, tout *= TMULT) {

        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (check_retval(&retval, "CVode", 1)) break;

        printf("t=%f\n", t);

        PrintVectorArray(&y, 1, 4);

        /* Call CVodeGetSens to get the sensitivity solution vector after a
           successful return from CVode */
        if (sensi) {
            retval = CVodeGetSens(cvode_mem, &t, yS);
            if (check_retval(&retval, "CVodeGetSens", 1)) break;
            PrintVectorArray(yS, 1, 4);
        }

    }

    /* Print final statistics */
    PrintFinalStats(cvode_mem, sensi);

    /* Free memory */

    N_VDestroy(y);                    /* Free y vector */
    if (sensi) {
        N_VDestroyVectorArray(yS, NS);  /* Free yS vector */
    }
    free(data);                              /* Free user data */
    CVodeFree(&cvode_mem);                   /* Free CVODES memory */
    SUNLinSolFree(LS);                       /* Free the linear solver memory */
    SUNMatDestroy(A);                        /* Free the matrix memory */

    return (0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY CVODES
 *--------------------------------------------------------------------
 */

/*
 * f routine. Compute f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {

    realtype scUPA, PLG, PLS, tcUPA;

    double *yData = y->ops->nvgetarraypointer(y);

    scUPA = yData[0];
    PLG = yData[1];
    PLS = yData[2];
    tcUPA = yData[3];

    double v1scUPAIn = alpha1;
    double v2scUPAOut = u1 * scUPA;
    double v3PLGIn = alpha2;
    double v4PLGOut = u2 * PLG;
    double v5PLGtoPLSByscUPA = keff1 * scUPA * PLG;
    double v6scUPATotcUPA = keff2 * scUPA * (pow(PLS, n));
    double v7PLGtoPLSBytcUPA = keff3 * tcUPA * PLG;
    double v8tcUPAOut = u3 * tcUPA;
    double v9PLSOut = u4 * PLS;

    double *ydotData = ydot->ops->nvgetarraypointer(ydot);

    ydotData[0] = +v1scUPAIn - v2scUPAOut - v6scUPATotcUPA;
    ydotData[1] = +v3PLGIn - v4PLGOut - v5PLGtoPLSByscUPA - v7PLGtoPLSBytcUPA;
    ydotData[2] = -v9PLSOut + v5PLGtoPLSByscUPA + v7PLGtoPLSBytcUPA;
    ydotData[3] = -v8tcUPAOut + v6scUPATotcUPA;

    return (0);
}


/* 
 * Jacobian routine. Compute J(t,y). 
 */

//static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
//               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
//{
//  realtype y2, y3;
//  UserData data;
//  realtype p1, p2, p3;
//
//  y2 = Ith(y,2); y3 = Ith(y,3);
//  data = (UserData) user_data;
//  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
//
//  IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
//  IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
//                      IJth(J,3,2) = 2*p3*y2;
//
//  return(0);
//}
//
/* 
 * fS routine. Compute sensitivity r.h.s. 
 */

//static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
//              int iS, N_Vector yS, N_Vector ySdot,
//              void *user_data, N_Vector tmp1, N_Vector tmp2)
//{
//  UserData data;
//  realtype p1, p2, p3;
//  realtype y1, y2, y3;
//  realtype s1, s2, s3;
//  realtype sd1, sd2, sd3;
//
//  data = (UserData) user_data;
//  p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
//
//  y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
//  s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);
//
//
//  // this is their equations
////  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
////  sd3 = 2*p3*y2*s2;
////  sd2 = -sd1-sd3;
//
//  // these are equivalent - without shortcuts.
//  sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
//  sd2 = p1*s1 - (p2*y3 + 2*p3*y2)*s2 - p2*y2*s3;
//  sd3 = 2*p3*y2*s2;
//
//  /**
//   * Solution:
//=======================================================================
//     T     Q       H      NST           y1           y2           y3
//=======================================================================
//4.000e-01  3  4.881e-02   115
//                  Solution         9.8517e-01   3.3864e-05   1.4794e-02
//                  Sensitivity 1   -3.5595e-01   3.9025e-04   3.5556e-01
//                  Sensitivity 2    9.5431e-08  -2.1309e-10  -9.5218e-08
//                  Sensitivity 3   -1.5833e-11  -5.2900e-13   1.6362e-11
//   */
//
//  switch (iS) {
//  case 0:
//    sd1 += -y1;
//    sd2 +=  y1;
//    break;
//  case 1:
//    sd1 +=  y2*y3;
//    sd2 += -y2*y3;
//    break;
//  case 2:
//    sd2 += -y2*y2;
//    sd3 +=  y2*y2;
//    break;
//  }
//
//  Ith(ySdot,1) = sd1;
//  Ith(ySdot,2) = sd2;
//  Ith(ySdot,3) = sd3;
//
//  return(0);
//}

/*
 * EwtSet function. Computes the error weights at the current solution.
 */

static int ewt(N_Vector y, N_Vector w, void *user_data) {
    int i;
    realtype yy, ww, rtol, atol[NEQ];

    rtol = RTOL;
    atol[0] = ATOL1;
    atol[1] = ATOL2;
    atol[2] = ATOL3;
    atol[3] = ATOL4;

    for (i = 1; i <= 3; i++) {
        yy = Ith(y, i);
        ww = rtol * SUNRabs(yy) + atol[i - 1];
        if (ww <= 0.0) return (-1);
        Ith(w, i) = 1.0 / ww;
    }

    return (0);
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Process and verify arguments to cvsfwddenx.
 */


/**
 *
 *

 This is outoput of the program when the rhs of the sensitivity equations are provided
 (i.e. not what we're aiming for in this program)
3-species chemical kinetics problem
Sensitivity: YES ( SIMULTANEOUS + FULL ERROR CONTROL )

=======================================================================
     T     Q       H      NST           y1           y2           y3
=======================================================================
4.000e-01  3  4.881e-02   115
                  Solution         9.8517e-01   3.3864e-05   1.4794e-02
                  Sensitivity 1   -3.5595e-01   3.9025e-04   3.5556e-01
                  Sensitivity 2    9.5431e-08  -2.1309e-10  -9.5218e-08
                  Sensitivity 3   -1.5833e-11  -5.2900e-13   1.6362e-11
-----------------------------------------------------------------------
4.000e+00  5  2.363e-01   138
                  Solution         9.0552e-01   2.2405e-05   9.4459e-02
                  Sensitivity 1   -1.8761e+00   1.7922e-04   1.8759e+00
                  Sensitivity 2    2.9614e-06  -5.8305e-10  -2.9608e-06
                  Sensitivity 3   -4.9334e-10  -2.7626e-13   4.9362e-10
-----------------------------------------------------------------------
4.000e+01  3  1.485e+00   219
                  Solution         7.1583e-01   9.1856e-06   2.8416e-01
                  Sensitivity 1   -4.2475e+00   4.5913e-05   4.2475e+00
                  Sensitivity 2    1.3731e-05  -2.3573e-10  -1.3730e-05
                  Sensitivity 3   -2.2883e-09  -1.1380e-13   2.2884e-09
-----------------------------------------------------------------------
4.000e+02  3  8.882e+00   331
                  Solution         4.5052e-01   3.2229e-06   5.4947e-01
                  Sensitivity 1   -5.9584e+00   3.5431e-06   5.9584e+00
                  Sensitivity 2    2.2738e-05  -2.2605e-11  -2.2738e-05
                  Sensitivity 3   -3.7896e-09  -4.9948e-14   3.7897e-09
-----------------------------------------------------------------------
4.000e+03  3  1.219e+02   490
                  Solution         1.8318e-01   8.9413e-07   8.1682e-01
                  Sensitivity 1   -4.7500e+00  -5.9948e-06   4.7500e+00
                  Sensitivity 2    1.8809e-05   2.3132e-11  -1.8809e-05
                  Sensitivity 3   -3.1348e-09  -1.8757e-14   3.1348e-09
-----------------------------------------------------------------------
4.000e+04  3  2.776e+03   598
                  Solution         3.8980e-02   1.6216e-07   9.6102e-01
                  Sensitivity 1   -1.5749e+00  -2.7620e-06   1.5749e+00
                  Sensitivity 2    6.2872e-06   1.1003e-11  -6.2872e-06
                  Sensitivity 3   -1.0479e-09  -4.5364e-15   1.0479e-09
-----------------------------------------------------------------------
4.000e+05  3  1.398e+04   665
                  Solution         4.9394e-03   1.9854e-08   9.9506e-01
                  Sensitivity 1   -2.3639e-01  -4.5856e-07   2.3639e-01
                  Sensitivity 2    9.4527e-07   1.8331e-12  -9.4527e-07
                  Sensitivity 3   -1.5753e-10  -6.3637e-16   1.5753e-10
-----------------------------------------------------------------------
4.000e+06  4  1.903e+05   716
                  Solution         5.1680e-04   2.0682e-09   9.9948e-01
                  Sensitivity 1   -2.5664e-02  -5.1053e-08   2.5664e-02
                  Sensitivity 2    1.0265e-07   2.0419e-13  -1.0265e-07
                  Sensitivity 3   -1.7109e-11  -6.8507e-17   1.7109e-11
-----------------------------------------------------------------------
4.000e+07  4  2.829e+06   771
                  Solution         5.2026e-05   2.0812e-10   9.9995e-01
                  Sensitivity 1   -2.5988e-03  -5.1935e-09   2.5988e-03
                  Sensitivity 2    1.0395e-08   2.0774e-14  -1.0395e-08
                  Sensitivity 3   -1.7326e-12  -6.9311e-18   1.7326e-12
-----------------------------------------------------------------------
4.000e+08  3  2.036e+07   816
                  Solution         5.1952e-06   2.0781e-11   9.9999e-01
                  Sensitivity 1   -2.6010e-04  -5.2113e-10   2.6010e-04
                  Sensitivity 2    1.0404e-09   2.0845e-15  -1.0404e-09
                  Sensitivity 3   -1.7315e-13  -6.9262e-19   1.7315e-13
-----------------------------------------------------------------------
4.000e+09  3  4.107e+08   849
                  Solution         5.1788e-07   2.0715e-12   1.0000e+00
                  Sensitivity 1   -2.5993e-05  -5.2215e-11   2.5993e-05
                  Sensitivity 2    1.0397e-10   2.0886e-16  -1.0397e-10
                  Sensitivity 3   -1.7262e-14  -6.9049e-20   1.7262e-14
-----------------------------------------------------------------------
4.000e+10  2  5.788e+09   878
                  Solution         5.0813e-08   2.0325e-13   1.0000e+00
                  Sensitivity 1   -2.4876e-06  -4.8465e-12   2.4876e-06
                  Sensitivity 2    9.9502e-12   1.9386e-17  -9.9502e-12
                  Sensitivity 3   -1.6938e-15  -6.7751e-21   1.6938e-15
-----------------------------------------------------------------------

Final Statistics

nst     =   878

nfe     =  1231
netf    =    30    nsetups  =   148
nni     =  1228    ncfn     =     4

nfSe    =  3693    nfeS     =     0
netfs   =     0    nsetupsS =     0
nniS    =     0    ncfnS    =     0

nje    =    24    nfeLS     =     0

Process finished with exit code 0

 */