#include <cvodes/cvodes.h>                        /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_serial.h>               /* access to serial N_Vector                    */
#include "sunnonlinsol/sunnonlinsol_newton.h"     /* access to the newton SUNNonlinearSolver      */


/**
 * Integrate, with sensitivities the following model:
 *
 * -> S1 -> S2 ->
 *
 * dS1/dt = k1*S1 - k2*S1;
 * dS2/dt = +k2*S2 - k3*S2
 * k1 = 0.3
 * k2 = 0.6
 * k3 = 0.01
 * S10 = 0
 * S20 = 0
 *
 *
 *
 */


#define K1 0.3
#define K2 0.6
#define K3 0.01
#define S10 0
#define S20 0

#define NEQ 2

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_retval(void *returnvalue, const char *funcname, int opt)
{
    int *retval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

        /* Check if retval < 0 */
    else if (opt == 1) {
        retval = (int *) returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                    funcname, *retval);
            return(1); }}

        /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    return(0);
}

#define CheckRetValNullPtr(thing, funcName)     if (check_retval((void*)thing, funcName, 0)) return 1;
#define CheckRetValNotNegative(thing, funcName)     if (check_retval((void*)thing, funcName, 1)) return 1;


int main(){

    N_Vector y = N_VNewEmpty_Serial(NEQ);
    CheckRetValNullPtr(y, "N_VNewEmpty_Serial")



















}























