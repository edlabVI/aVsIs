#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <math.h>

#include "../include/dq_3phase_function.h"

/* PARK's DIRECT TRANSFORMATION */
void parkDirect(float (*vect)[2], float *Sthme, float *Cthme)
{
	float temp_ABVect[2] 	= {0.0,0.0};
	temp_ABVect[0]			= (*vect)[0];
	temp_ABVect[1]			= (*vect)[1];
	(*vect)[0]				= *Cthme*temp_ABVect[0] + *Sthme*temp_ABVect[1];
	(*vect)[1]				= -*Sthme*temp_ABVect[0] + *Cthme*temp_ABVect[1];
}

/* PARK's INVERSE TRANSFORMATION */
void parkInverse(float (*vect)[2], float *Sthme, float *Cthme)
{
	float temp_dqVect[2] 	= {0.0,0.0};
	temp_dqVect[0]			= (*vect)[0];
	temp_dqVect[1]			= (*vect)[1];
	(*vect)[0]				= *Cthme*temp_dqVect[0] - *Sthme*temp_dqVect[1];
	(*vect)[1]				= *Sthme*temp_dqVect[0] + *Cthme*temp_dqVect[1];
}

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


/* CONSTRAINT 1 AND 4 EVALUATION */
void constrONEandFOUR(float (*uTemp)[2],float (*uOptTemp)[2], float *uBus_sca, float *costOptTemp, float *costOpt, int *flagOptSol, float (*H)[2][2], float (*f)[2], float *Sthme, float *Cthme)
{
	float SplusSQRT3C 	= *Sthme+SQRT_3*(*Cthme);
	float CminusSQRT3S 	= *Cthme-SQRT_3*(*Sthme);
	
	if(!(fabsf(CminusSQRT3S)<MINIMUM_CHK)){
		float term[5] 		= {0.0};
		term[0]				= 1.0/CminusSQRT3S;
		term[1]				= term[0]*SplusSQRT3C;
		term[2]				= (*H)[0][1]+(*H)[1][0];
		term[3]				= *uBus_sca*ONE_OVER_SQRT_3;
		term[4]				= *uBus_sca*TWO_OVER_SQRT_3;
		term[5]				= term[2]*term[3];
		
		float coeff_1_4[5]  ={0.0};
		coeff_1_4[0]		= (*f)[1]*term[1];
		coeff_1_4[1]		= term[5]*term[0];
		coeff_1_4[2]		= term[4]*(*H)[1][1]*term[1]*term[0];
		coeff_1_4[3]		= term[4]*term[0];
		coeff_1_4[4]		= term[1];
		
		float den_1_4		= (*H)[0][0] - term[2]*term[1] + (*H)[1][1]*term[1]*term[1];
		
		
		/* Get first solution */
		(*uTemp)[0]                    = -((*f)[0]-coeff_1_4[0]+coeff_1_4[1]-coeff_1_4[2])/den_1_4;
		(*uTemp)[1]                    = coeff_1_4[3] - coeff_1_4[4]*(*uTemp)[0];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (TWO_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < (ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
		/* Get fourth solution */
		(*uTemp)[0]                    = -((*f)[0]-coeff_1_4[0]-coeff_1_4[1]+coeff_1_4[2])/den_1_4;
		(*uTemp)[1]                    = -coeff_1_4[3] - coeff_1_4[4]*(*uTemp)[0];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > -(ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(TWO_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
	}
	else
	{
		float term[5] 		= {0.0};
		term[0]     		= 1.0/SplusSQRT3C;
		term[1]				= term[0]*CminusSQRT3S;
		term[2]				= (*H)[0][1]+(*H)[1][0];
		term[3]				= *uBus_sca*ONE_OVER_SQRT_3;
		term[4]				= *uBus_sca*TWO_OVER_SQRT_3;
		term[5]				= term[2]*term[3];
		
		float coeff_1_4[5]  = {0.0};
		coeff_1_4[0]		= (*f)[0]*term[1];
		coeff_1_4[1]		= term[5]*term[0];
		coeff_1_4[2]		= term[4]*(*H)[0][0]*term[1]*term[0];
		
		coeff_1_4[3]		= term[4]*term[0];
		coeff_1_4[4]		= term[1];
		
		float den_1_4		= (*H)[1][1] - term[2]*term[1] + (*H)[0][0]*term[1]*term[1];
				
		/* Get first solution */
		(*uTemp)[1]                    = -((*f)[1]-coeff_1_4[0]+coeff_1_4[1]-coeff_1_4[2])/den_1_4;
		(*uTemp)[0]                    = coeff_1_4[3] - coeff_1_4[4]*(*uTemp)[1];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (TWO_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < (ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
		/* Get fourth solution */
		(*uTemp)[1]                    = -((*f)[1]-coeff_1_4[0]-coeff_1_4[1]+coeff_1_4[2])/den_1_4;
		(*uTemp)[0]                    = -coeff_1_4[3] - coeff_1_4[4]*(*uTemp)[1];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > -(ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(TWO_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
	}
}

/* CONSTRAINT 2 AND 5 EVALUATION */
void constrTWOandFIVE(float (*uTemp)[2],float (*uOptTemp)[2], float *uBus_sca, float *costOptTemp, float *costOpt, int *flagOptSol, float (*H)[2][2], float (*f)[2], float *Sthme, float *Cthme)
{
	if(!(fabsf(*Cthme)<MINIMUM_CHK)){
		
		float SdivC			= *Sthme/(*Cthme);
		float ONEdivC		= 1.0/(*Cthme);
		
		float term[3]		={0.0};
		term[0]				= (*H)[0][1]+(*H)[1][0];
		term[1]				= *uBus_sca*ONE_OVER_SQRT_3;
		term[2]				= *uBus_sca*ONE_OVER_2SQRT_3*term[0];
		
		float coeff_2_5[5]	={0.0};
		coeff_2_5[0]		= SdivC*(*f)[1];
		coeff_2_5[1]		= term[2]*ONEdivC;
		coeff_2_5[2]		= term[1]*(*H)[1][1]*SdivC*ONEdivC;
		coeff_2_5[3]		= term[1]*ONEdivC;
		coeff_2_5[4]		= SdivC;
		
		float den_2_5		= (*H)[0][0] - SdivC*term[0] + SdivC*SdivC*(*H)[1][1];
	
				
		/* Get second solution */
		(*uTemp)[0]                    = -((*f)[0]-coeff_2_5[0]+coeff_2_5[1]-coeff_2_5[2])/den_2_5;
		(*uTemp)[1]                    = coeff_2_5[3] - coeff_2_5[4]*(*uTemp)[0];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
		/* Get fifth solution */
		(*uTemp)[0]                    = -((*f)[0]-coeff_2_5[0]-coeff_2_5[1]+coeff_2_5[2])/den_2_5;
		(*uTemp)[1]                    = -coeff_2_5[3] - coeff_2_5[4]*(*uTemp)[0];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
	}
	else
	{
		float CdivS			= *Cthme/(*Sthme);
		float ONEdivS		= 1.0/(*Sthme);
		
		float term[3]		={0.0};
		term[0]				= (*H)[0][1]+(*H)[1][0];
		term[1]				= *uBus_sca*ONE_OVER_SQRT_3;
		term[2]				= *uBus_sca*ONE_OVER_2SQRT_3*term[0];
		
		float coeff_2_5[5]	={0.0};
		coeff_2_5[0]		= CdivS*(*f)[0];
		coeff_2_5[1]		= term[2]*ONEdivS;
		coeff_2_5[2]		= term[1]*(*H)[0][0]*CdivS*ONEdivS;
		coeff_2_5[3]		= term[1]*ONEdivS;
		coeff_2_5[4]		= CdivS;
		
		float den_2_5		= (*H)[1][1] - CdivS*term[0] + CdivS*CdivS*(*H)[0][0];

	
		/* Get second solution */
		(*uTemp)[1]                    = -((*f)[1]-coeff_2_5[0]+coeff_2_5[1]-coeff_2_5[2])/den_2_5;
		(*uTemp)[0]                    = coeff_2_5[3] - coeff_2_5[4]*(*uTemp)[1];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
		/* Get fifth solution */
		(*uTemp)[1]                    = -((*f)[1]-coeff_2_5[0]-coeff_2_5[1]+coeff_2_5[2])/den_2_5;
		(*uTemp)[0]                    = -coeff_2_5[3] - coeff_2_5[4]*(*uTemp)[1];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
	}
}

/* CONSTRAINT 3 AND 6 EVALUATION */
void constrTHREEandSIX(float (*uTemp)[2],float (*uOptTemp)[2], float *uBus_sca, float *costOptTemp, float *costOpt, int *flagOptSol, float (*H)[2][2], float (*f)[2], float *Sthme, float *Cthme)
{
	float SminusSQRT3C 	= *Sthme-SQRT_3*(*Cthme);
	float CplusSQRT3S 	= *Cthme+SQRT_3*(*Sthme);
	
	if(!(fabsf(CplusSQRT3S)<MINIMUM_CHK)){
		
		float term[5] 		= {0.0};
		term[0]				= 1.0/CplusSQRT3S;
		term[1]				= term[0]*SminusSQRT3C;
		term[2]				= (*H)[0][1]+(*H)[1][0];
		term[3]			= *uBus_sca*ONE_OVER_SQRT_3;
		term[4]			= *uBus_sca*TWO_OVER_SQRT_3;
		term[5]			= term[2]*term[3];
		
		float coeff_3_6[5]	={0.0};
		coeff_3_6[0]		= (*f)[1]*term[1];
		coeff_3_6[1]		= term[5]*term[0];
		coeff_3_6[2]		= term[4]*(*H)[1][1]*term[1]*term[0];
		
		coeff_3_6[3]		= term[4]*term[0];
		coeff_3_6[4]		= term[1];
		
		float den_3_6		= (*H)[0][0] - term[2]*term[1] + (*H)[1][1]*term[1]*term[1];
	
				
		
		/* Get third solution */
		(*uTemp)[0]                    = -((*f)[0]-coeff_3_6[0]+coeff_3_6[1]-coeff_3_6[2])/den_3_6;
		(*uTemp)[1]                    = coeff_3_6[3] - coeff_3_6[4]*(*uTemp)[0];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > -(ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(TWO_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
		/* Get sixth solution */
		(*uTemp)[0]                    = -((*f)[0]-coeff_3_6[0]-coeff_3_6[1]+coeff_3_6[2])/den_3_6;
		(*uTemp)[1]                    = -coeff_3_6[3] - coeff_3_6[4]*(*uTemp)[0];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (TWO_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < (ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
	}
	else
	{
		float SminusSQRT3C 	= *Sthme-SQRT_3*(*Cthme);
		float CplusSQRT3S 	= *Cthme+SQRT_3*(*Sthme);
		
		float term[6] 		= {0.0};
		term[0]				= 1.0/SminusSQRT3C;
		term[1]				= term[0]*CplusSQRT3S;
		
		term[2]				= (*H)[0][1]+(*H)[1][0];
		term[3]				= *uBus_sca*ONE_OVER_SQRT_3;
		term[4]				= *uBus_sca*TWO_OVER_SQRT_3;
		term[5]				= term[2]*term[3];
		
		float coeff_3_6[5]	={0.0};
		coeff_3_6[0]		= (*f)[0]*term[1];
		coeff_3_6[1]		= term[1]*term[0];
		coeff_3_6[2]		= term[4]*(*H)[0][0]*term[1]*term[0];
		
		coeff_3_6[3]		= term[4]*term[0];
		coeff_3_6[4]		= term[1];
		
		float den_3_6		= (*H)[1][1] - term[2]*term[1] + (*H)[0][0]*term[1]*term[4];
				
		/* Get third solution */
		(*uTemp)[1]                    = -((*f)[1]-coeff_3_6[0]+coeff_3_6[1]-coeff_3_6[2])/den_3_6;
		(*uTemp)[0]                    = coeff_3_6[3] - coeff_3_6[4]*(*uTemp)[1];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > -(ONE_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < -(TWO_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
		
		/* Get sixth solution */
		(*uTemp)[1]                    = -((*f)[1]-coeff_3_6[0]-coeff_3_6[1]+coeff_3_6[2])/den_3_6;
		(*uTemp)[0]                    = -coeff_3_6[3] - coeff_3_6[4]*(*uTemp)[1];
		parkInverse(uTemp,Sthme,Cthme);
		if(!((*uTemp)[0] > (TWO_OVER_THREE*(*uBus_sca)) || (*uTemp)[0] < (ONE_OVER_THREE*(*uBus_sca))))
		{
			parkDirect(uTemp,Sthme,Cthme);
			checkSolution(uTemp,uOptTemp,costOptTemp,costOpt,flagOptSol, H, f);
		}
	}
}

/* CHECK CONSTRAINTS */
char checkConstraints(float (*vec)[2], float *uBus_sca)
{
    char tempOut    = 0; /*Compliant with SESE concept*/
    
    /* Check all constraints on the input vector "vec"*/
    for(int ii=0; ii<6; ii++)
    {
        printf("conVer = %f \n", (A_constr[ii][0]*(*vec)[0] + A_constr[ii][1]*(*vec)[1] - *uBus_sca*b_constr[ii]));
		//printf("conVer = %f \n", (A_constr[ii][0]*(*vec)[0]));
		if((A_constr[ii][0]*(*vec)[0] + A_constr[ii][1]*(*vec)[1] - *uBus_sca*b_constr[ii]) > 0.0){
			tempOut = 1;
			printf("verificare le costrizioni \n");
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