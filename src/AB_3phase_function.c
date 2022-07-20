#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>

#include "../include/AB_3phase_function.h"


/* CHECK COST OF THE HEXAGON SIDE SOLUTIONS */
void checkSolution(float (*tempSol)[2], float (*tempOptSol)[2], float *costTempSol, float *costTempOptSol, int *flag, float (*H)[2][2], float (*f)[2])
{
 	*costTempSol            = costFunctionEval(H,f,tempSol);
	if(*flag == 1)
	{
		if(*costTempSol < *costTempOptSol)
		{
			/* The temporary optimal solution present a lower cost 
			 * compared to the previous one
			 */
			*costTempOptSol     = *costTempSol;
			/* Update temporary optimal solution vector */
			(*tempOptSol)[0]       = (*tempSol)[0];
			(*tempOptSol)[1]       = (*tempSol)[1];
		}
	}
	else
	{
		*costTempOptSol     = *costTempSol;
		/* Update temporary optimal solution vector */
		(*tempOptSol)[0]       = (*tempSol)[0];
		(*tempOptSol)[1]       = (*tempSol)[1];
		/* Rise up the "flagOptSol" flag*/
		*flag               = 1;                
	}
}




/* GET COEFFICIENTS FOR THE SIX SOLUTIONS IN ALPHA-BETA REFERENCE FRAME*/
void getCoefficientsSixSolutions(float (*H)[2][2], float (*f)[2], float *uBus_sca, float (*coeffSol_vec)[8])
{
    /* Get coefficients*/
    (*coeffSol_vec)[0]     = 2.0*(*uBus_sca)*(*H)[1][1];
    (*coeffSol_vec)[1]     = SQRT_3*(*f)[1];
    /* get a temporary value for the third coefficients: it spare one addition
     for calculating the fourth and fifth  coefficient*/
    (*coeffSol_vec)[2]     = ((*H)[0][1] + (*H)[1][0]);
    (*coeffSol_vec)[3]     = 1.0/((*H)[0][0] + 3.0*(*H)[1][1] - SQRT_3*(*coeffSol_vec)[2]);
    (*coeffSol_vec)[4]     = 1.0/((*H)[0][0] + 3.0*(*H)[1][1] + SQRT_3*(*coeffSol_vec)[2]);
    /* Complete the calculation of the third coefficient */
    (*coeffSol_vec)[2]     = (*uBus_sca)*(*coeffSol_vec)[2]*ONE_OVER_SQRT_3;
    (*coeffSol_vec)[5]     = 2.0*(*uBus_sca)*ONE_OVER_SQRT_3;
    (*coeffSol_vec)[6]     = (*uBus_sca)*ONE_OVER_SQRT_3;
    (*coeffSol_vec)[7]     = 1.0/((*H)[0][0]);
}



/* CHECK CONSTRAINTS */
char checkConstraints(float (*vec)[2], float *uBus_sca)
{
    char tempOut    = 0; /*Compliant with SESE concept*/
    
    /* Check all constraints on the input vector "vec"*/
    for(int ii=0; ii<6; ii++)
    {
        if((A_constr[ii][0]*(*vec)[0] + A_constr[ii][1]*(*vec)[1] - *uBus_sca*b_constr[ii]) > 0.0){
			tempOut = 1;
		}
    }
    
    return tempOut;
    
}

/* COST FUNCTION EVALUATION */
float costFunctionEval(float (*H)[2][2], float (*f)[2], float (*u)[2])
{
    /* No values of H, f nor u are modified */
    return (0.5*((*u)[0]*(*u)[0]*(*H)[0][0] + (*u)[1]*(*u)[1]*(*H)[1][1] + ((*H)[0][1]+(*H)[1][0])*(*u)[0]*(*u)[1]) + (*u)[0]*(*f)[0] + (*u)[1]*(*f)[1]);
}