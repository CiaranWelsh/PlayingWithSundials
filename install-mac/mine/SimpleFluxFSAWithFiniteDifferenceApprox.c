#include <cvodes/cvodes.h>                        /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include <math.h>
#include "helpers.h"

/**
 * Example not using user supplied SA rhs:
 * /Users/ciaranwelsh/Documents/sundials/install-mac/examples/cvodes/serial/cvsAdvDiff_FSA_non.c
 */

/**
 * Integrate, with sensitivities the following model:
 *
 * -> S1 -> S2 ->
 *
 * dS1/dt = k1 - k2*S1;
 * dS2/dt = +k2*S2 - k3*S2
 * k1 = 0.3
 * k2 = 0.6
 * k3 = 0.01
 * S10 = 0
 * S20 = 0
 *
 */


#define Kin 1
#define Kf 0.1
#define Kout 0.2
#define S10 0.0
#define S20 0.0

#define NEQ 2
#define NP 3 // num parameters in the model
#define NS 3 // Number of parameters we want sensitivities for

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */


typedef struct {
    double p[NS];
}* UserData;


/**
 *
 * -> S1 -> S2 ->
 *
 * dS1/dt = k1 - k2*S1;
 * dS2/dt = +k2*S2 - k3*S2
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data) {
    realtype s1, s2;

    s1 = NV_Ith_S(y, 0);
    s2 = NV_Ith_S(y, 1);

    NV_Ith_S(ydot, 0) = Kin - Kf * s1;
    NV_Ith_S(ydot, 1) = Kf * s1 - Kout * s2;

    return (0);
}


static void PrintOutput(realtype t, realtype y0, realtype y1, int qu, realtype hu) {
    printf("%10.5f    %12.5e   %12.5e   %2d    %6.4e\n", t, y0, y1, qu, hu);

    return;
}

static void PrintuS(N_Vector * nVectorArray){
    PrintVectorArray(nVectorArray, NS, NEQ);
}

static int ewt(N_Vector y, N_Vector w, void *user_data) {
    int i;
    realtype yy, ww, rtol, atol[2];

    rtol = 1e-3;
    atol[0] = 1e-6;
    atol[1] = 1e-6;

    for (i = 1; i <= 2; i++) {
        yy = Ith(y, i);
        ww = rtol * fabs(yy) + atol[i - 1];
        if (ww <= 0.0) return (-1);
        Ith(w, i) = 1.0 / ww;
    }

    return (0);
}

int main() {

    // use fixed point or newton solver, (mutually exclusive)
    int fixed_point = 0;
    int newton = 1;

    // with sensitivity or not
    int with_sens = 1;

    // 0 or 1. Use error control or not
    int err_con = 1;

    // place to store return value. err haandling
    int retval;

    UserData userData ;
    userData = (UserData) malloc(sizeof* userData);

    // sensitivity method
    int sensi_meth = CV_SIMULTANEOUS;

    if (fixed_point && newton) {
        printf("Can't have both fixed point and newton at the same time", stderr);
        return 1;
    }

    // solvers
    SUNLinearSolver LS;
    SUNNonlinearSolver NLS;
    SUNNonlinearSolver NLSsens;
    SUNMatrix A;

    N_Vector y = N_VNew_Serial(NEQ);
    CheckRetValNotNull(y, "N_VNewEmpty_Serial")

//    void *cvode_mem = CVodeCreate(CV_BDF);
    void *cvode_mem = CVodeCreate(CV_ADAMS);
    CheckRetValNotNull(cvode_mem, "CVodeCreate")

    /****************
     * Set initial conditions
     */
    NV_Ith_S(y, 0) = S10;
    NV_Ith_S(y, 1) = S20;

    // initialize cvode
    retval = CVodeInit(cvode_mem, f, 0.0, y);
    CheckRetValNotNegative(&retval, "CVodeInit")

    // set scalar tolerances
    retval = CVodeSStolerances(cvode_mem, 1e-6, 1e-12);
//    retval = CVodeWFtolerances(cvode_mem, ewt);
    CheckRetValNotNegative(&retval, "CVodeSStolerances")

    if (fixed_point) {
        NLS = SUNNonlinSol_FixedPoint(y, 0);
        CheckRetValNotNull(NLS, "SUNNonlinSol_FixedPoint");
    }
    if (newton) {
        NLS = SUNNonlinSol_Newton(y);
        CheckRetValNotNull(NLS, "SUNNonlinSol_Newton");
    }

    /* attach nonlinear solver object to CVode */
    retval = CVodeSetNonlinearSolver(cvode_mem, NLS);
    CheckRetValNotNegative(&retval, "CVodeSetNonlinearSolver");

    /* Create dense SUNMatrix for use in linear solves */
    A = SUNDenseMatrix(NEQ, NEQ);
    if (check_retval((void *) A, "SUNDenseMatrix", 0)) return (1);

    /* Create dense SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_Dense(y, A);
    CheckRetValNotNull(LS, "SUNLinSol_Dense");

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    CheckRetValNotNegative(&retval, "CVodeSetLinearSolver");

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(cvode_mem, NULL);
    CheckRetValNotNegative(&retval, "CVodeSetJacFn");

    int is;
    int plist[NP] = {0, 1, 2};
    double pbar[NP] = {Kin, Kf, Kout};
    userData->p[0] = Kin;
    userData->p[1] = Kf;
    userData->p[2] =  Kout;

    retval = CVodeSetUserData(cvode_mem, userData);
    CheckRetValNotNegative(&retval, "CVodeSetUserData")

    /**
     * y = state vector. shape = (num eqns)
     * uS = shape (num parameters x num eqns)
     *
     */
    N_Vector *uS;
    N_Vector *uDydk;
    N_Vector *uDky;
    if (with_sens) {
        // which sens

        // pointer to N_Vector - i.e. a matrix
        uS = N_VCloneVectorArray_Serial(NS, y);
        uDydk = N_VCloneVectorArray_Serial(NS, y);
        CheckRetValNotNull(uS, "N_VCloneVectorArray");

        setSensInitTo(uS, NS, NEQ, 5.0);


        printf("Initial uS:\n");
        PrintuS(uS);

        retval = CVodeSensInit1(cvode_mem, NS, sensi_meth, NULL, uS);

        double reltol = 1e-6;
        double abstol = 1e-8;
        retval = CVodeSensSStolerances(cvode_mem, reltol, &abstol);
        CheckRetValNotNegative(&retval, "CVodeSensSStolerances");


        retval = CVodeSetSensErrCon(cvode_mem, err_con);
        CheckRetValNotNegative(&retval, "CVodeSetSensErrCon");

        retval = CVodeSetSensDQMethod(cvode_mem, CV_FORWARD, 1);
        CheckRetValNotNegative(&retval, "CVodeSetSensDQMethod");

        // plist is only needed if we are being selective about which parameters
        // to compute sensitivities for.
        // typically we can use the parameter values for pbar
        retval = CVodeSetSensParams(cvode_mem, userData->p, NULL, NULL);
        printf("retval for CVodeSetSensParams: %d\n", retval);
        CheckRetValNotNegative(&retval, "CVodeSetSensParams");


        if (sensi_meth == CV_SIMULTANEOUS) {
            if (newton) {
                NLSsens = SUNNonlinSol_NewtonSens(NS + 1, y);
            }
            if (fixed_point) {
                NLSsens = SUNNonlinSol_FixedPointSens(NS + 1, y, 1);
            }
        } else if (sensi_meth == CV_STAGGERED) {
            if (newton) {
                NLSsens = SUNNonlinSol_NewtonSens(NS, y);
            } else {
                NLSsens = SUNNonlinSol_FixedPointSens(NS, y, 0);
            }
        }
        else {
            if (newton){
                NLSsens = SUNNonlinSol_Newton(y);
            } else {
                NLSsens = SUNNonlinSol_FixedPoint(y, 0);
            }
        }

        CheckRetValNotNull(NLSsens, "Creating NLSSens");

        /* attach nonlinear solver object to CVode */
        if (sensi_meth == CV_SIMULTANEOUS)
            retval = CVodeSetNonlinearSolverSensSim(cvode_mem, NLSsens);
        else if(sensi_meth == CV_STAGGERED)
            retval = CVodeSetNonlinearSolverSensStg(cvode_mem, NLSsens);
        else
            retval = CVodeSetNonlinearSolverSensStg1(cvode_mem, NLSsens);
        CheckRetValNotNegative(&retval, "Attaching nonlinear sens solver")

    }

    int i, temp_retval, qu;
    int nerr = 0;
    double t, tout, hu;
    t = 0.0;
    for (i = 0, tout = 1.0; i < 10; i++, tout += 1.0) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        CheckRetValNotNegative(&retval, "CVode");
        printf("retval for integration: %d\n", retval);
        temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
        if (check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
        temp_retval = CVodeGetLastStep(cvode_mem, &hu);
        if (check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
        PrintOutput(t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), qu, hu);
        if (retval != CV_SUCCESS) {
            nerr++;
            break;
        }

        if (with_sens) {
            retval = CVodeGetSens(cvode_mem, &t, uS);
            CheckRetValNotNegative(&retval, "CVodeGetSens");
            printf("retval for sens: %d\n", retval);

            retval = CVodeGetSensDky(cvode_mem, t, 1, uDydk);
            CheckRetValNotNegative(&retval, "CVodeGetSensDky");

            printf("uS: \n");
            PrintuS(uS);
            printf("uDydk: \n");
            PrintuS(uDydk);

        }
    }
    if (with_sens){
        long nfsevals = 0;
        // number of calls to the sensitivity function
        CVodeGetSensNumRhsEvals(cvode_mem, &nfsevals);
        printf("number of calls to the sensitivity funcition : %ld\n", nfsevals);

        CVodeGetNumRhsEvalsSens(cvode_mem,  &nfsevals);
        printf("number of calls to usres sens function due to finite differences : %ld\n", nfsevals);
    }

    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    NLS = NULL;
    LS = NULL;
    A = NULL;
}




/*
 * 0.2255 0.0741
 *
Initial uS:

			 S1 		 S2
 	 kin	 0.000000 	 0.000000
	 kf		 0.000000 	 0.000000
	 kout	 0.000000 	 0.000000

   1.00000     9.51626e-01    4.52796e-02    6    3.8757e-01
uS:

			 S1 		 S2
 	 kin	 0.951626 	 0.051796
	 kf		 -0.467884 	 0.501268
	 kout	 0.000000 	 -0.016291

uDydk:

			 S1 		 S2
 	 kin	 0.904837 	 0.105522
	 kf		 -0.904837 	 1.005091
	 kout	 0.000000 	 -0.048538

   2.00000     1.81269e+00    1.64293e-01    6    3.8757e-01
uS:

			 S1 		 S2
 	 kin	 1.812692 	 0.215477
	 kf		 -1.752310 	 2.020617
	 kout	 0.000000 	 -0.127961

uDydk:

			 S1 		 S2
 	 kin	 0.818731 	 0.224365
	 kf		 -1.637462 	 2.041585
	 kout	 0.000000 	 -0.189885

   3.00000     2.59182e+00    3.35876e-01    6    3.8757e-01
uS:

			 S1 		 S2
 	 kin	 2.591818 	 0.506259
	 kf		 -3.693631 	 4.606269
	 kout	 0.000000 	 -0.425957

uDydk:

			 S1 		 S2
 	 kin	 0.740818 	 0.360434
	 kf		 -2.222455 	 3.143708
	 kout	 0.000000 	 -0.421067

   4.00000     3.29680e+00    5.43444e-01    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 3.296800 	 0.943635
	 kf		 -6.155194 	 8.342631
	 kout	 0.000000 	 -1.000477

uDydk:

			 S1 		 S2
 	 kin	 0.670320 	 0.518407
	 kf		 -2.681280 	 4.349806
	 kout	 0.000000 	 -0.743540

   5.00000     3.93469e+00    7.74091e-01    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 3.934693 	 1.552239
	 kf		 -9.020401 	 13.355057
	 kout	 0.000000 	 -1.945370

uDydk:

			 S1 		 S2
 	 kin	 0.606531 	 0.703917
	 kf		 -3.032653 	 5.703665
	 kout	 0.000000 	 -1.163165

   6.00000     4.51188e+00    1.01785e+00    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 4.511884 	 2.362900
	 kf		 -12.190138 	 19.816048
	 kout	 0.000000 	 -3.362614

uDydk:

			 S1 		 S2
 	 kin	 0.548812 	 0.923768
	 kf		 -3.292870 	 7.256078
	 kout	 0.000000 	 -1.690377

   7.00000     5.03415e+00    1.26713e+00    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 5.034147 	 3.413951
	 kf		 -15.580498 	 27.953169
	 kout	 0.000000 	 -5.367047

uDydk:

			 S1 		 S2
 	 kin	 0.496585 	 1.186204
	 kf		 -3.476097 	 9.066724
	 kout	 0.000000 	 -2.340539

   8.00000     5.50671e+00    1.51619e+00    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 5.506710 	 4.752815
	 kf		 -19.120787 	 38.059029
	 kout	 0.000000 	 -8.091555

uDydk:

			 S1 		 S2
 	 kin	 0.449329 	 1.501232
	 kf		 -3.594632 	 11.206426
	 kout	 0.000000 	 -3.134499

   9.00000     5.93430e+00    1.76080e+00    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 5.934303 	 6.437974
	 kf		 -22.751765 	 50.503746
	 kout	 0.000000 	 -11.692939

uDydk:

			 S1 		 S2
 	 kin	 0.406570 	 1.881024
	 kf		 -3.659127 	 13.759867
	 kout	 0.000000 	 -4.099382

  10.00000     6.32121e+00    1.99788e+00    7    1.1056e+00
uS:

			 S1 		 S2
 	 kin	 6.321206 	 8.541352
	 kf		 -26.424113 	 65.750380
	 kout	 0.000000 	 -16.358673

uDydk:

			 S1 		 S2
 	 kin	 0.367879 	 2.340389
	 kf		 -3.678795 	 16.828858
	 kout	 0.000000 	 -5.269613

number of calls to the sensitivity funcition : 150
number of calls to usres sens function due to finite differences : 0

 */


















