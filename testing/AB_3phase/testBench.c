#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>

#include "../../include/AB_3phase_problem.h"


void main(){
	#define HEIGHT 640

	FILE *f_file;
	FILE *uAB_file;
	int i;

	float h1 = 0.00215312500000000;
	float h2 = 0.0;
	float h3 = 0.0;
	float h4 = 0.00215312500000000;

	float H[2][2] = {
		{h1, h2},
		{h3, h4}
	};

	float f1;
	float f2;

	float uA_ref_MAT;
	float uB_ref_MAT;

	float Ubus = 50.0;
	float uOpt[2] ={0.0, 0.0}; 

	f_file=fopen("data/f_txt.txt", "r");
	uAB_file=fopen("data/uAB_ref.txt", "r");

	for(int i = 0; i < HEIGHT; i++)
	{
		printf("----------------------- \n",uA_ref_MAT);
		
		fscanf(f_file,"%f %f",&f1,&f2);
		//printf("%.15f ",f1);
		//printf("%.15f \n",f2);
		float f[2] = {f1, f2};
		
		fscanf(uAB_file,"%f %f",&uA_ref_MAT,&uB_ref_MAT);
		
		aVsIs_AB3phase(&H,&f,&Ubus,&uOpt);
		
		printf("UA_FILE= %.15f \n",uA_ref_MAT);
		printf("UB_FILE= %.15f \n",uB_ref_MAT);
		
		printf("uA_OPT = %.15f \n", uOpt[0]);
		printf("uB_OPT = %.15f \n \n", uOpt[1]);
		
		printf("diff_uA = %f \n", uA_ref_MAT-uOpt[0]);
		printf("diff_uB = %f \n \n", uB_ref_MAT-uOpt[1]);
		
		if(abs(uA_ref_MAT-uOpt[0])>0.0001 || abs(uB_ref_MAT-uOpt[1])>0.0001){
			printf("Press Any Key to Continue\n");
			printf("index = %d \n", i);
			getchar(); 
		}
	}
	printf("Code works \n");
	printf("\n");


	fclose(f_file);
	fclose(uAB_file);

}