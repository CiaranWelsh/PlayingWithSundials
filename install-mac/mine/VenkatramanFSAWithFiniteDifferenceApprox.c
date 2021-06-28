#include <cvodes/cvodes.h>                        /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
#include <math.h>

#include "helpers.h"



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

#define NEQ 4 // num equations.
#define NP 10 // num parameters in the model
#define NS 1  // Number of parameters we want sensitivities for

#define Ith(v, i)    NV_Ith_S(v,i-1)         /* i-th vector component i=1..NEQ */

static void PrintuS(N_Vector *nVectorArray) {
    PrintVectorArray(nVectorArray, NS, NEQ);
}

static void PrintOutput(realtype t, realtype y0, realtype y1, realtype y2, realtype y3, int qu, realtype hu) {
    if (t == 0.0) {
        printf("  t               scUPA         PLG          PLS              tcUPA       qu        hu \n");
        printf("================================================================================================\n");
    }
    printf("%10.5f    %12.5e   %12.5e   %12.5e  %12.5e   %2d    %6.4e\n", t, y0, y1, y2, y3, qu, hu);

    return;
}

typedef struct {
    double p[NP];
} *UserData;


/**
 *
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

static int computeInitSens(N_Vector y0, N_Vector y1, double step, N_Vector *yS) {
    int N1 = y0->ops->nvgetlength(y0);
    int N2 = y1->ops->nvgetlength(y1);
    if (N1 != N2) {
        printf("computeInitSens y0 dimensions are not the same as y2\n", stderr);
    }
    double *y0data = y0->ops->nvgetarraypointer(y0);
    double *y1data = y1->ops->nvgetarraypointer(y1);
    double *ySdata = y1->ops->nvgetarraypointer(*yS);
    for (int i = 0; i < y0->ops->nvgetlength(y0); i++) {
        ySdata[i] = (fabs(y1data[i] - y0data[i])) / step;
    }
    return 0;
}


static int integrate(void *cvode_mem, N_Vector y, N_Vector *yS, double start, double stop, int numIntegrationSteps,
                     int with_sens) {
    printf("yS on entry to integrate:\n");
    y->ops->nvprint(*yS);
    int retval = 0, nerr = 0, qu = 0;
    double hu = 0.0;
    double stepsize = (stop - start) / numIntegrationSteps;
    double cvodeInternalTime = 0.0;
    double tout = start + stepsize;
    for (int stepNum = 0; stepNum < numIntegrationSteps; stepNum++, start += stepsize, tout += stepsize) {
        printf("yS on entry to integration loop:\n");
        y->ops->nvprint(*yS);

        retval = CVode(cvode_mem, tout, y, &cvodeInternalTime, CV_NORMAL);

        printf("yS after call to cvode:\n");
        y->ops->nvprint(*yS);

        CheckRetValNotNegative(&retval, "CVode");
        retval = CVodeGetLastOrder(cvode_mem, &qu);
        if (check_retval(&retval, "CVodeGetLastOrder", 1)) ++nerr;
        retval = CVodeGetLastStep(cvode_mem, &hu);
        if (check_retval(&retval, "CVodeGetLastStep", 1)) ++nerr;
        PrintOutput(cvodeInternalTime, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3), qu, hu);
        if (retval != CV_SUCCESS) {
            printf("EERRRRRORRRRRR\n", stderr);
            nerr++;
            break;
        }

        if (with_sens) {
            printf("yS before call to get sens:\n");
            y->ops->nvprint(*yS);

            retval = CVodeGetSens(cvode_mem, &cvodeInternalTime, yS);
            CheckRetValNotNegative(&retval, "CVodeGetSens");
            printf("retval for CVodeGetSens : %d\n", retval);

            printf("yS after call to get sens:\n");
            y->ops->nvprint(*yS);

//            retval = CVodeGetSensDky(cvode_mem, t, 1, uDydk);
//            CheckRetValNotNegative(&retval, "CVodeGetSensDky");

            printf("\nyS for param at time %f: size: %lld\n", cvodeInternalTime, yS[0]->ops->nvgetlength(yS[0]));
            yS[0]->ops->nvprint(*yS);
//            PrintuS(yS);
//            printf("uDydk: \n");
//            PrintuS(uDydk);

        }
    }
    if (with_sens) {
        long nfsevals = 0;
        CVodeGetSensNumRhsEvals(cvode_mem, &nfsevals);
        printf("No. of calls to sensitivity r.h.s. function:     %ld \n", nfsevals);
        CVodeGetNumRhsEvalsSens(cvode_mem, &nfsevals);
        printf("No. of calls to r.h.s. function for sensitivity "
               "(only incremented when internal finite differences routine is used):"
               "     %ld \n", nfsevals);
        CVodeGetSensNumErrTestFails(cvode_mem, &nfsevals);
        printf("No. of sensitivity local error test failures:        %ld \n", nfsevals);
        CVodeGetSensNumLinSolvSetups(cvode_mem, &nfsevals);
        printf("No. of calls to lin. solv. setup routine for sens.:      %ld \n", nfsevals);
//        CVodeGetSensErrWeights(cvode_mem, &nfsevals);
//        printf("Error weight vector for sensitivity variables:       %f \n")
        CVodeGetSensNumNonlinSolvIters(cvode_mem, &nfsevals);
        printf("No. of sens. nonlinear solver iterations:        %ld \n", nfsevals);
        CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &nfsevals);
        printf("No. of sens. convergence failures:       %ld \n", nfsevals);
        CVodeGetStgrSensNumNonlinSolvIters(cvode_mem, &nfsevals);
        printf("No. of staggered nonlinear solver iterations %ld \n", nfsevals);
        long fails = 0;
        CVodeGetStgrSensNumNonlinSolvConvFails(cvode_mem, &fails);
        printf("No. of staggered convergence failures %ld \n", fails);
    }
    return 0;
}


int main() {

    // use fixed point or newton solver, (mutually exclusive)
    int fixed_point = 0;
    int newton = 1;

    // with sensitivity or not
    int with_sens = 1;

    // 0 or 1. Use error control or not
    int err_con = 0;

    // place to store return value. err haandling
    int retval;

    UserData userData;
    userData = (UserData) malloc(sizeof *userData);

    // sensitivity method
    int sensi_meth = CV_STAGGERED;

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

    void *cvode_mem = CVodeCreate(CV_BDF);
//    void *cvode_mem = CVodeCreate(CV_ADAMS);
    CheckRetValNotNull(cvode_mem, "CVodeCreate")

    /****************
     * Set initial conditions
     */
    double *yData = y->ops->nvgetarraypointer(y);
    yData[0] = scUPA0;
    yData[1] = PLG0;
    yData[2] = PLS0;
    yData[3] = tcUPA0;

    // initialize cvode
    retval = CVodeInit(cvode_mem, f, 0.0, y);
    CheckRetValNotNegative(&retval, "CVodeInit")

    // set tolerances
//    retval = CVodeSensEEtolerances(cvode_mem);
    retval = CVodeSStolerances(cvode_mem, 1e-6, 1e-12);
//    retval = CVodeWFtolerances(cvode_mem, ewt);
    CheckRetValNotNegative(&retval, "CVodeSStolerances")

    if (fixed_point) {
        NLS = SUNNonlinSol_FixedPoint(y, 1);
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

    retval = CVodeSetMaxNumSteps(cvode_mem, 100000);
    CheckRetValNotNegative(&retval, "CVodeSetMaxNumSteps");

    int is;
    int plistAll[NP] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    double pAll[NP] = {keff1, keff2, keff3, n, u1, u3, u2, u4, alpha1, alpha2};

    // which parameter to do sensitivities for?
    int which = 0;

    // we just do one parameter for now
    int plist[NS] = {plistAll[which]};
    double pbar[NS] = {pAll[which]};
    const char *pName[NP] = {
            "keff1",
            "keff2",
            "keff3",
            "n",
            "u1",
            "u3",
            "u2",
            "u4",
            "alpha1",
            "alpha2",
    };

//    userData->p[0] = pAll[which]; // just do one parameter for now
    userData->p[0] = keff1; // just do one parameter for now
    userData->p[1] = keff2;
    userData->p[2] = keff3;
    userData->p[3] = n;
    userData->p[4] = u1;
    userData->p[5] = u3;
    userData->p[6] = u2;
    userData->p[7] = u4;
    userData->p[8] = alpha1;
    userData->p[9] = alpha2;

    retval = CVodeSetUserData(cvode_mem, userData);
    CheckRetValNotNegative(&retval, "CVodeSetUserData")

    /**
     * y = state vector. shape = (num eqns)
     * uS = shape (num parameters x num eqns)
     *
     */
    N_Vector *yS;
    N_Vector *dySdk;
    if (with_sens) {
        // which sens

        // pointer to N_Vector - i.e. a matrix
        yS = N_VCloneVectorArray_Serial(NS, y);
        dySdk = N_VCloneVectorArray_Serial(NS, y); // kth derivative of sensitivites
        CheckRetValNotNull(yS, "N_VCloneVectorArray");


        retval = CVodeSensInit(cvode_mem, NS, sensi_meth, NULL, yS);

        retval = CVodeSensEEtolerances(cvode_mem);
//        double reltol = 1e-6;
//        double abstol = 1e-12;
//        retval = CVodeSensSStolerances(cvode_mem, reltol, &abstol);
        CheckRetValNotNegative(&retval, "CVodeSensSStolerances");

        retval = CVodeSetSensErrCon(cvode_mem, err_con);
        CheckRetValNotNegative(&retval, "CVodeSetSensErrCon");

        retval = CVodeSetSensDQMethod(cvode_mem, CV_FORWARD, 0);
        CheckRetValNotNegative(&retval, "CVodeSetSensDQMethod");

        retval = CVodeSetSensMaxNonlinIters(cvode_mem, 10000);
        CheckRetValNotNegative(&retval, "CVodeSetSensMaxNonlinIters");

        // plist is only needed if we are being selective about which parameters
        // to compute sensitivities for.
        // typically we can use the parameter values for pbar
        retval = CVodeSetSensParams(cvode_mem, userData->p, pbar, NULL);
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
        } else {
            if (newton) {
                NLSsens = SUNNonlinSol_Newton(y);
            } else {
                NLSsens = SUNNonlinSol_FixedPoint(y, 0);
            }
        }

        CheckRetValNotNull(NLSsens, "Creating NLSSens");

        /* attach nonlinear solver object to CVode */
        if (sensi_meth == CV_SIMULTANEOUS)
            retval = CVodeSetNonlinearSolverSensSim(cvode_mem, NLSsens);
        else if (sensi_meth == CV_STAGGERED)
            retval = CVodeSetNonlinearSolverSensStg(cvode_mem, NLSsens);
        else
            retval = CVodeSetNonlinearSolverSensStg1(cvode_mem, NLSsens);
        CheckRetValNotNegative(&retval, "Attaching nonlinear sens solver")

    }

    int numSteps, temp_retval, qu;
    int nerr = 0;
    double t, tout, hu;
    t = 0.0;

    // prints ICs
    PrintOutput(t, NV_Ith_S(y, 0), NV_Ith_S(y, 1), NV_Ith_S(y, 2), NV_Ith_S(y, 3), 0, 0);

    N_Vector y0 = N_VClone(y);
    N_Vector y1 = N_VClone(y);
    N_Vector yS0 = N_VClone(y);

    for (int i = 0; i < y->ops->nvgetlength(y); i++) {
        NV_Ith_S(y0, i) = NV_Ith_S(y, i);
        NV_Ith_S(y1, i) = NV_Ith_S(y, i);
        NV_Ith_S(yS0, i) = NV_Ith_S(y, i);
    }

    double firstStep = 0.001;
    integrate(cvode_mem, y1, yS, 0, firstStep, 1, 0);
    y0->ops->nvprint(y0);
    y1->ops->nvprint(y1);

    computeInitSens(y0, y1, firstStep, &yS0);

    printf("yS[0]\n");
    yS[0]->ops->nvprint(yS0);
    integrate(cvode_mem, y, &yS0, 0, 100, 10, 1);

    CVodeFree(&cvode_mem);
    SUNNonlinSolFree(NLS);
    NLS = NULL;
    LS = NULL;
    A = NULL;
}


















