#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>

#include "../include/AB_3phase_problem.h"
#include "../include/AB_3phase_function.h"

void aVsIs_AB3phase(float (*H)[2][2], float (*f)[2], float *uBus_sca, float (*out_vec)[2])
{
    /* Get the unconstrained solution */
    /* Calculate the inverse of the matrix H. It is suppose that H is positive*/
	
	float invDet        = 1.0/((*H)[0][0]*(*H)[1][1] - (*H)[0][1]*(*H)[1][0]);
      
	float invH[2][2]    = {{invDet*(*H)[1][1], -invDet*(*H)[1][0]}, {-invDet*(*H)[0][1], invDet*(*H)[0][0]}};
    	
    (*out_vec)[0]          = -(invH[0][0]*(*f)[0] + invH[0][1]*(*f)[1]);
    (*out_vec)[1]          = -(invH[1][0]*(*f)[0] + invH[1][1]*(*f)[1]); 
    /* Check constraints*/
    if(checkConstraints(out_vec,uBus_sca))
    {	
		/* Voltage saturation occurred. A new optimal solution has to be found*/
        /* First of all, get the coefficients for calculating the six
         different solutions and store them into the temporary vector 
         coeffSolTemp_vec */
        float coeffSolTemp_vec[8]   = {0.0};
        getCoefficientsSixSolutions(H,f,uBus_sca,&coeffSolTemp_vec);
        
        /* Store the temporary solution into the temp variable uTemp*/
        float uTemp[2]              = {0.0};
        /* Store the temporary optimal solution into the temp variable uOptTemp */
        float uOptTemp[2]           = {0.0};
        
        /* Calculate edge solutions, one for each side of the hexagon, and store
         the cost into the costTempComp. The temporary minimum cost is stored
         in the costMinTemp variable */
        float costOptTemp           = 0.0;
        float costOpt               = 0.0;
        
        /* use a flag to highlight that a possible and feasible solution has 
         * been calculated: when flagOptSol == 1, then at least one possible 
         * solution has been obtained
         */
        int flagOptSol              = 0;
        
        /* Get first solution */
        uTemp[0]                    = -((*f)[0]-coeffSolTemp_vec[0]-coeffSolTemp_vec[1]+coeffSolTemp_vec[2])*coeffSolTemp_vec[3];
        uTemp[1]                    = coeffSolTemp_vec[5] - SQRT_3*uTemp[0];
		if(!(uTemp[0] > (TWO_OVER_THREE*(*uBus_sca)) || uTemp[0] < (ONE_OVER_THREE*(*uBus_sca))))
        {
			checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Get second solution */
        uTemp[0]                    = -((*f)[0]+0.5*coeffSolTemp_vec[2])*coeffSolTemp_vec[7];
        uTemp[1]                    = *uBus_sca*ONE_OVER_SQRT_3;
		if(!(uTemp[0] > (ONE_OVER_THREE*(*uBus_sca)) || uTemp[0] < -(ONE_OVER_THREE*(*uBus_sca))))
        {
			checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Get third solution */
        uTemp[0]                    = -((*f)[0]+coeffSolTemp_vec[0]+coeffSolTemp_vec[1]+coeffSolTemp_vec[2])*coeffSolTemp_vec[4];
        uTemp[1]                    = coeffSolTemp_vec[5] + SQRT_3*uTemp[0];
		if(!(uTemp[0] > -(ONE_OVER_THREE*(*uBus_sca)) || uTemp[0] < -(TWO_OVER_THREE*(*uBus_sca))))
        {
			checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Get fourth solution */
        uTemp[0]                    = -((*f)[0]+coeffSolTemp_vec[0]-coeffSolTemp_vec[1]-coeffSolTemp_vec[2])*coeffSolTemp_vec[3];
        uTemp[1]                    = -coeffSolTemp_vec[5] - SQRT_3*uTemp[0];
		if(!(uTemp[0] > -(ONE_OVER_THREE*(*uBus_sca)) || uTemp[0] < -(TWO_OVER_THREE*(*uBus_sca))))
        {
			checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Get fifth solution */
        uTemp[0]                    = -((*f)[0]-0.5*coeffSolTemp_vec[2])*coeffSolTemp_vec[7];
        uTemp[1]                    = -*uBus_sca*ONE_OVER_SQRT_3;
		if(!(uTemp[0] > (ONE_OVER_THREE*(*uBus_sca)) || uTemp[0] < -(ONE_OVER_THREE*(*uBus_sca))))
        {
			checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Get sixth solution */
        uTemp[0]                 = -((*f)[0]-coeffSolTemp_vec[0]+coeffSolTemp_vec[1]-coeffSolTemp_vec[2])*coeffSolTemp_vec[4];
        uTemp[1]                 = -coeffSolTemp_vec[5] + SQRT_3*uTemp[0];
		if(!(uTemp[0] > (TWO_OVER_THREE*(*uBus_sca)) || uTemp[0] < (ONE_OVER_THREE*(*uBus_sca))))
        {
	        checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Check the vertex of the hexagon */
        for(int ii=0; ii<6; ii++)
        {
            uTemp[0]            = (*uBus_sca)*uAlphaVertex[ii];
            uTemp[1]            = (*uBus_sca)*uBetaVertex[ii];
	        /* being the vertex of the hexagon, no check on the feasibility is
             necessary. Go straight to the cost evaluation. In the worst case,
             * i.e. when none of the six solutions calculated on the edges of 
             * the hexagon are feasible, one solution on the vertex shall be
             * certainly applied */
            checkSolution(&uTemp,&uOptTemp, &costOptTemp, &costOpt, &flagOptSol, H, f);
        }
        
        /* Write the optimal solution to the output*/
        (*out_vec)[0]              = uOptTemp[0];
        (*out_vec)[1]              = uOptTemp[1];
    }


}