#include <cvodes/cvodes.h>                        /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */
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

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
              int iS, N_Vector yS, N_Vector ySdot,
              void *user_data, N_Vector tmp1, N_Vector tmp2)
{
    UserData data;
    realtype p1, p2, p3;
    realtype y1, y2;
    realtype s1, s2;
    realtype sd1, sd2;

    data = (UserData) user_data;
    // unpack parameter values from user data
    p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2]; // kin kf and kout

    // unpack values of S1 and S2 at time t
    y1 = Ith(y,1);  y2 = Ith(y,2); // S1, S2

    // unpack the current values of the sensitivity vector Si
    s1 = Ith(yS,1); s2 = Ith(yS,2); //S1, S2

    // Jacobian * Si part of sensitivity equations.
    sd1 = -Kf*s1;
    sd2 = Kf*s1 + Kout*s2;

    switch (iS) {
        case 0:
            sd1 += 1;
            sd2 += 0;
            break;
        case 1:
            sd1 +=  -y1;
            sd2 +=   y1;
            break;
        case 2:
            sd1 +=  0;
            sd2 += -y2;
            break;
    }

    Ith(ySdot,1) = sd1;
    Ith(ySdot,2) = sd2;

    return(0);
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
    int is;
    int plist[NP] = {0, 1, 2};
    double pbar[NP] = {Kin, Kf, Kout};
    userData->p[0] = Kin;
    userData->p[1] = Kf;
    userData->p[2] =  Kout;

    // sensitivity method
//    int sensi_meth = CV_SIMULTANEOUS;
    int sensi_meth = CV_STAGGERED;

    if (fixed_point && newton) {
        printf("Can't have both fixed point and newton at the same time", stderr);
        return 1;
    }

    // solvers
    SUNNonlinearSolver NLS;
    SUNNonlinearSolver NLSsens;

    N_Vector y = N_VNew_Serial(NEQ);
    CheckRetValNotNull(y, "N_VNewEmpty_Serial")

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

    retval = CVodeSetUserData(cvode_mem, userData);
    CheckRetValNotNegative(&retval, "CVodeSetUserData")

    // set scalar tolerances
    retval = CVodeSStolerances(cvode_mem, 1e-6, 1e-12);
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
    SUNMatrix A = SUNDenseMatrix(NEQ, NEQ);
    if (check_retval((void *) A, "SUNDenseMatrix", 0)) return (1);

    /* Create dense SUNLinearSolver object for use by CVode */
    SUNLinearSolver LS = SUNLinSol_Dense(y, A);
    CheckRetValNotNull(LS, "SUNLinSol_Dense");

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    CheckRetValNotNegative(&retval, "CVodeSetLinearSolver");

    /* Set the user-supplied Jacobian routine Jac */
//    retval = CVodeSetJacFn(cvode_mem, NULL);
//    CheckRetValNotNegative(&retval, "CVodeSetJacFn");



    /**
     * y = state vector. shape = (1 x num eqns)
     * yS = shape (num parameters x num eqns)
     *
     */
    N_Vector *yS;
    N_Vector *uDydk;
    N_Vector *uDky;
    if (with_sens) {
        // which sens

        // pointer to N_Vector - i.e. a matrix
        yS = N_VCloneVectorArray_Serial(NS, y);
        uDydk = N_VCloneVectorArray_Serial(NS, y);
        CheckRetValNotNull(yS, "N_VCloneVectorArray");

        // set sens init cond
        setSensInitTo(yS, NS, NEQ, 0.0);

        printf("Initial yS:\n");
        PrintVectorArray(yS, NS, NEQ);

        retval = CVodeSensInit(cvode_mem, NS, sensi_meth, NULL, yS);

        retval = CVodeSensEEtolerances(cvode_mem);
        CheckRetValNotNegative(&retval, "CVodeSensEEtolerances");

        retval = CVodeSetSensErrCon(cvode_mem, err_con);
        CheckRetValNotNegative(&retval, "CVodeSetSensErrCon");

        retval = CVodeSetSensDQMethod(cvode_mem, CV_FORWARD, 1);
        CheckRetValNotNegative(&retval, "CVodeSetSensDQMethod");

        retval = CVodeSetSensParams(cvode_mem, userData->p, pbar, plist);
        CheckRetValNotNegative(&retval, "CVodeSetSensParams");

        SUNNonlinSol_NewtonSens(NS, y);

        if (sensi_meth == CV_SIMULTANEOUS) {
            if (newton) {
                NLSsens = SUNNonlinSol_NewtonSens(NS + 1, y);
            }
            if (fixed_point) {
                NLSsens = SUNNonlinSol_FixedPointSens(NS + 1, y, 0);
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
    N_Vector dky = N_VClone(y);
    for (i = 0, tout = 1.0; i < 10; i++, tout += 1.0) {
        retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
        CheckRetValNotNegative(&retval, "CVode");

        CVodeGetDky(cvode_mem, t, 1, dky);

        temp_retval = CVodeGetLastOrder(cvode_mem, &qu);
        if (check_retval(&temp_retval, "CVodeGetLastOrder", 1)) ++nerr;
        temp_retval = CVodeGetLastStep(cvode_mem, &hu);
        if (check_retval(&temp_retval, "CVodeGetLastStep", 1)) ++nerr;
        printf("\t\t solution at time t= %f\n", t);
        PrintVectorArray(&y, 1, y->ops->nvgetlength(y));
        if (retval != CV_SUCCESS) {
            nerr++;
            break;
        }

        if (with_sens) {
            retval = CVodeGetSens(cvode_mem, &t, yS);
            CheckRetValNotNegative(&retval, "CVodeGetSens");

            retval = CVodeGetSensDky(cvode_mem, t, 1, uDydk);
            CheckRetValNotNegative(&retval, "CVodeGetSensDky");

            printf("\t\tyS: \n");
            PrintVectorArray(yS, NS, NEQ);
            printf("\t\t first derivative: \n");
            PrintVectorArray(uDydk, NS, NEQ);

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
 * 1.00000     2.25594e-01    7.41466e-02    5    2.1578e-01
   2.00000     3.49403e-01    2.48784e-01    5    2.1578e-01
   3.00000     4.17350e-01    4.77238e-01    5    3.4978e-01
   4.00000     4.54640e-01    7.33909e-01    5    2.5686e-01
   5.00000     4.75105e-01    1.00476e+00    5    2.5686e-01
   6.00000     4.86338e-01    1.28209e+00    5    1.8582e-01
   7.00000     4.92502e-01    1.56171e+00    5    1.3618e-01
   8.00000     4.95885e-01    1.84131e+00    5    1.3618e-01
   9.00000     4.97742e-01    2.11965e+00    5    1.3618e-01
  10.00000     4.98760e-01    2.39605e+00    5    1.3618e-01
 */


















