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

/* Accessor macros */

#define Ith(v, i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */
#define IJth(A, i, j) SM_ELEMENT_D(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Problem Constants */
#define Kin 1
#define Kf 0.1
#define Kout 0.2
#define S10 0.0
#define S20 0.0

#define NEQ 2
#define NP 3 // num parameters in the model
#define NS 3 // Number of parameters we want sensitivities for

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */


/* Type : UserData */

typedef struct {
    realtype p[3];           /* problem parameters */
} *UserData;

/* Prototypes of functions by CVODES */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

//static int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
//               void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

//static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
//              int iS, N_Vector yS, N_Vector ySdot,
//              void *user_data, N_Vector tmp1, N_Vector tmp2);

//static int ewt(N_Vector y, N_Vector w, void *user_data);

/* Prototypes of private functions */

static void ProcessArgs(int argc, char *argv[],
                        booleantype *sensi, int *sensi_meth,
                        booleantype *err_con);

static void WrongArgs(char *name);

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);

static void PrintOutputS(N_Vector *uS);

static void PrintFinalStats(void *cvode_mem, booleantype sensi);

static int check_retval(void *returnvalue, const char *funcname, int opt);

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

    /* Process arguments */
//  ProcessArgs(argc, argv, &sensi, &sensi_meth, &err_con);
    sensi_meth = CV_SIMULTANEOUS;
    err_con = 1;

    /* User data structure */
    data = (UserData) malloc(sizeof *data);
    if (check_retval((void *) data, "malloc", 2)) return (1);
    data->p[0] = Kin;
    data->p[1] = Kf;
    data->p[2] = Kout;

    /* Initial conditions */
    y = N_VNew_Serial(NEQ);
    if (check_retval((void *) y, "N_VNew_Serial", 0)) return (1);

    Ith(y, 1) = S10;
    Ith(y, 2) = S20;

    /* Create CVODES object */
    cvode_mem = CVodeCreate(CV_BDF);
    if (check_retval((void *) cvode_mem, "CVodeCreate", 0)) return (1);

    /* Allocate space for CVODES */
    retval = CVodeInit(cvode_mem, f, 0.0, y);
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

    /* Sensitivity-related settings */
    if (sensi) {

        /* Set parameter scaling factor */
        pbar[0] = data->p[0];
        pbar[1] = data->p[1];
        pbar[2] = data->p[2];

        /* Set sensitivity initial conditions */
        yS = N_VCloneVectorArray(NS, y);
        if (check_retval((void *) yS, "N_VCloneVectorArray", 0)) return (1);
        for (is = 0; is < NS; is++) N_VConst(0, yS[is]);

        /* Call CVodeSensInit1 to activate forward sensitivity computations
           and allocate internal memory for COVEDS related to sensitivity
           calculations. Computes the right-hand sides of the sensitivity
           ODE, one at a time */
        retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, yS);
        if (check_retval(&retval, "CVodeSensInit", 1)) return (1);

        /* Call CVodeSensEEtolerances to estimate tolerances for sensitivity
           variables based on the rolerances supplied for states variables and
           the scaling factor pbar */
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
    printf("===========================================");
    printf("============================\n");
    printf("     T     Q       H      NST           y1");
    printf("           y2           y3    \n");
    printf("===========================================");
    printf("============================\n");

    for (iout = 1, tout = 1.0; iout <= 10.0; iout++) {

        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        if (check_retval(&retval, "CVode", 1)) break;

        PrintOutput(cvode_mem, t, y);

        /* Call CVodeGetSens to get the sensitivity solution vector after a
           successful return from CVode */
        if (sensi) {
            retval = CVodeGetSens(cvode_mem, &t, yS);
            if (check_retval(&retval, "CVodeGetSens", 1)) break;
            PrintOutputS(yS);
        }
        printf("-----------------------------------------");
        printf("------------------------------\n");

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
    realtype s1, s2;

    s1 = NV_Ith_S(y, 0);
    s2 = NV_Ith_S(y, 1);

    NV_Ith_S(ydot, 0) = Kin - Kf * s1;
    NV_Ith_S(ydot, 1) = Kf * s1 - Kout * s2;

    return (0);
}

/*
 * Print current t, step count, order, stepsize, and solution.
 */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u) {
    long int nst;
    int qu, retval;
    realtype hu, *udata;

    udata = N_VGetArrayPointer(u);

    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1);
    retval = CVodeGetLastOrder(cvode_mem, &qu);
    check_retval(&retval, "CVodeGetLastOrder", 1);
    retval = CVodeGetLastStep(cvode_mem, &hu);
    check_retval(&retval, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#else
    printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

    printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#else
    printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#endif

}

/*
 * Print sensitivities.
*/

static void PrintOutputS(N_Vector *uS) {
    realtype *sdata;

    sdata = N_VGetArrayPointer(uS[0]);
    printf("                  Sensitivity 1  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

    sdata = N_VGetArrayPointer(uS[1]);
    printf("                  Sensitivity 2  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif

    sdata = N_VGetArrayPointer(uS[2]);
    printf("                  Sensitivity 3  ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf("%12.4Le %12.4Le %12.4Le \n", sdata[0], sdata[1], sdata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#else
    printf("%12.4e %12.4e %12.4e \n", sdata[0], sdata[1], sdata[2]);
#endif
}

/*
 * Print some final statistics from the CVODES memory.
 */

static void PrintFinalStats(void *cvode_mem, booleantype sensi) {
    long int nst;
    long int nfe, nsetups, nni, ncfn, netf;
    long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
    long int nje, nfeLS;
    int retval;

    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1);
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_retval(&retval, "CVodeGetNumRhsEvals", 1);
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_retval(&retval, "CVodeGetNumErrTestFails", 1);
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

    if (sensi) {
        retval = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
        check_retval(&retval, "CVodeGetSensNumRhsEvals", 1);
        retval = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
        check_retval(&retval, "CVodeGetNumRhsEvalsSens", 1);
        retval = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
        check_retval(&retval, "CVodeGetSensNumLinSolvSetups", 1);
        retval = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
        check_retval(&retval, "CVodeGetSensNumErrTestFails", 1);
        retval = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
        check_retval(&retval, "CVodeGetSensNumNonlinSolvIters", 1);
        retval = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
        check_retval(&retval, "CVodeGetSensNumNonlinSolvConvFails", 1);
    }

    retval = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_retval(&retval, "CVodeGetNumJacEvals", 1);
    retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
    check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

    printf("\nFinal Statistics\n\n");
    printf("nst     = %5ld\n\n", nst);
    printf("nfe     = %5ld\n", nfe);
    printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
    printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

    if (sensi) {
        printf("\n");
        printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
        printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
        printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
    }

    printf("\n");
    printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);

}

/*
 * Check function return value.
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns an integer value so check if
 *             retval < 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt) {
    int *retval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr,
                "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return (1);
    }

        /* Check if retval < 0 */
    else if (opt == 1) {
        retval = (int *) returnvalue;
        if (*retval < 0) {
            fprintf(stderr,
                    "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                    funcname, *retval);
            return (1);
        }
    }

        /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr,
                "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return (1);
    }

    return (0);
}
