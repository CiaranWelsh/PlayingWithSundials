//
// Created by Ciaran Welsh on 25/06/2021.
//

#ifndef SUNDIALS_HELPERS_H
#define SUNDIALS_HELPERS_H


/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns an integer value so check if
              retval < 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_retval(void *returnvalue, const char *funcname, int opt) {
    int *retval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return (1);
    }

        /* Check if retval < 0 */
    else if (opt == 1) {
        retval = (int *) returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                    funcname, *retval);
            return (1);
        }
    }

        /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return (1);
    }

    return (0);
}

#define CheckRetValNotNull(thing, funcName)     if (check_retval((void*)thing, funcName, 0)) return 1;
#define CheckRetValNotNegative(thing, funcName)     if (check_retval((void*)thing, funcName, 1)) return 1;





static void PrintVectorArray(N_Vector *nVectorArray, int nrow, int ncol) {
    realtype *sdata;
    for (int i=0; i<nrow; i++){
        sdata = N_VGetArrayPointer(nVectorArray[i]);
        for (int j=0; j<ncol; j++){
            printf("%.10f\t", sdata[j]);
        }
        printf("\n");
    }
}

static void setSensInitTo(N_Vector* nVectorArray, int numParamsYouWantSensFor, int numModelVariables, double value){
    for (int is = 0; is < numParamsYouWantSensFor ; is++) {
        double *data = nVectorArray[0]->ops->nvgetarraypointer(nVectorArray[is]);
        // then over num model species
        for (int species=0; species < numModelVariables; species++){
            // set sens init to 0.
            data[species] = value;
        }
    }
}


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


#endif //SUNDIALS_HELPERS_H
