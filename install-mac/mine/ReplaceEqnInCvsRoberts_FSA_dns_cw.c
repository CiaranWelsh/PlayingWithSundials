
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
#define NS NP  // Number of parameters we want sensitivities for

#define ATOL  RCONST(1.e-5) /* scalar absolute tolerance */
#define T0    RCONST(0.0)   /* initial time              */
#define T1    RCONST(0.5)   /* first output time         */
#define DTOUT RCONST(0.5)   /* output time increment     */
#define NOUT  10            /* number of output times    */
#define TMULT 10            /* output time factor */
#define STEP  1

#define RTOL  RCONST(1e-4)  /* scalar relative tolerance */

// used in error weight fn, which isn't being used right now
#define ATOL1 RCONST(1e-8)  /* vector absolute tolerance components */
#define ATOL2 RCONST(1e-7)
#define ATOL3 RCONST(1e-6)
#define ATOL4 RCONST(1e-6)

//#define keff1    0.0017    // param 0
//#define keff1    0.006829    // param 0
#define keff1    0.001717    // param 0


#define keff2    1         // param 1
#define keff3    0.03      // param 2
#define n        3         // param 3
#define u1        0.0001   // param 4
#define u3        0.0001   // param 5
#define u2        0.001    // param 6
#define u4        0.001    // param 7
#define alpha1    9e-5     // param 8
#define alpha2    0.001    // param 9


#define scUPA0    1.16213e-8
#define PLG0    0.0262792
#define PLS0    19.7847
#define tcUPA0    19.9809

#define SUNDIALS_EXTENDED_PRECISION


/* Type : UserData */

typedef struct {
    realtype p[NS];           /* problem parameters */
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
    int plist[NS];
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

    sensi = 0;

    /* Process arguments */
//  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);
//    sensi_meth = CV_SIMULTANEOUS;
    sensi_meth = CV_STAGGERED;
//    err_con = 1;
    err_con = 0;

    /* User data structure */
    data = (UserData) malloc(sizeof *data);
    if (check_retval((void *) data, "malloc", 2)) return (1);
    // fill the parameters in the user data struct
    data->p[0] = keff1;
    data->p[1] = keff2;
    data->p[2] = keff3;
    data->p[3] = n;
    data->p[4] = u1;
    data->p[5] = u3;
    data->p[6] = u2;
    data->p[7] = u4;
    data->p[8] = alpha1;
    data->p[9] = alpha2;

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
    retval = CVodeSStolerances(cvode_mem, 1e-4, 1e-5);
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


    /* Sensitivity-related settings */
    if (sensi) {

        /* Set parameter scaling factor */
        int which = 1; // index of parameter for sensitivities
        pbar[0] = data->p[0];
        pbar[1] = data->p[1];
        pbar[2] = data->p[2];
        pbar[3] = data->p[3];
        pbar[4] = data->p[4];
        pbar[5] = data->p[5];
        pbar[6] = data->p[6];
        pbar[7] = data->p[7];
        pbar[8] = data->p[8];
        pbar[9] = data->p[9];

        for (int i=0; i<10; i++) plist[i] = i;

        /* Set sensitivity initial conditions */
        yS = N_VCloneVectorArray(NS, y);
        if (check_retval((void *) yS, "N_VCloneVectorArray", 0)) return (1);

        /* sets all sensitivity initial conditions to value */
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
        retval = CVodeSetSensParams(cvode_mem, data->p, pbar, plist);
        if (check_retval(&retval, "CVodeSetSensParams", 1)) return (1);

        printf("Sensitivity: YES ");
        if (sensi_meth == CV_SIMULTANEOUS){
                printf(" (SIMULTANEOUS +\n");
        }
        else if (sensi_meth == CV_STAGGERED) printf("( STAGGERED +");
        else printf("( STAGGERED1 +");
        if (err_con) printf(" FULL ERROR CONTROL )");
        else printf(" PARTIAL ERROR CONTROL )");

    } else {

        printf("Sensitivity: NO ");

    }

    /* In loop over output points, call CVode, print results, test for error */

    printf("\n\n");

//    for (iout = 1, tout = T1; iout <= NOUT; iout++, tout *= TMULT) {
    for (iout = 1, tout = T1; iout <= NOUT; iout++, tout += STEP) {

        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (check_retval(&retval, "CVode", 1)) break;

        printf("\t\tt=%f\n", t);

        PrintVectorArray(&y, 1, NEQ); // y is a 1xNEQ matrix or column vector

        /* Call CVodeGetSens to get the sensitivity solution vector after a
           successful return from CVode */
        if (sensi) {
            retval = CVodeGetSens(cvode_mem, &t, yS);
            if (check_retval(&retval, "CVodeGetSens", 1)) break;
            PrintVectorArray(yS, NS, NEQ);
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
 * Solution with keff1 = 0.0016829, i.e. theta - delta theta:
 *
 * Sensitivity: NO

		t=0.500000
0.0000000116	0.0199012863	19.7816748144	19.9799459792
		t=1.500000
0.0000000116	0.0116357950	19.7711481278	19.9780380800
		t=2.500000
0.0000000117	0.0071386582	19.7568718687	19.9761303716
		t=3.500000
0.0000000117	0.0046789354	19.7405769598	19.9742228539
		t=4.500000
0.0000000117	0.0033235157	19.7231964932	19.9723155269
		t=5.500000
0.0000000118	0.0025792079	19.7052236411	19.9704083907
		t=6.500000
0.0000000118	0.0021692445	19.6869351525	19.9685014452
		t=7.500000
0.0000000118	0.0019413777	19.6684832563	19.9665946904
		t=8.500000
0.0000000119	0.0018131983	19.6499503432	19.9646881262
		t=9.500000
0.0000000119	0.0017285186	19.6313925821	19.9627817527

Final Statistics

nst     =    25

nfe     =    39
netf    =     2    nsetups  =    10
nni     =    36    ncfn     =     0

nje    =     1    nfeLS     =     4


  * Solution with keff1 = 0.001717, i.e. theta - delta theta:


		t=0.500000
0.0000000116	0.0199012863	19.7816748144	19.9799459792
		t=1.500000
0.0000000116	0.0116357950	19.7711481278	19.9780380800
		t=2.500000
0.0000000117	0.0071386582	19.7568718687	19.9761303716
		t=3.500000
0.0000000117	0.0046789354	19.7405769598	19.9742228539
		t=4.500000
0.0000000117	0.0033235157	19.7231964932	19.9723155269
		t=5.500000
0.0000000118	0.0025792079	19.7052236411	19.9704083907
		t=6.500000
0.0000000118	0.0021692445	19.6869351525	19.9685014452
		t=7.500000
0.0000000118	0.0019413777	19.6684832563	19.9665946904
		t=8.500000
0.0000000119	0.0018131983	19.6499503432	19.9646881262
		t=9.500000
0.0000000119	0.0017285186	19.6313925821	19.9627817527

Final Statistics

nst     =    25

nfe     =    39
netf    =     2    nsetups  =    10
nni     =    36    ncfn     =     0

nje    =     1    nfeLS     =     4






 */