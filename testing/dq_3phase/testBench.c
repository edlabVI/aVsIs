#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <math.h>

#include "../../include/dq_3phase_problem.h"

void main(){
#define HEIGHT 640

FILE *thme_file;
FILE *H_file;
FILE *f_file;
FILE *udq_file;

int i;

float Ubus = 280.0;
float uOpt[2] ={0.0, 0.0}; 

thme_file	=fopen("data/thme_txt.txt", "r");
H_file		=fopen("data/H_txt.txt", "r");
f_file		=fopen("data/f_txt.txt", "r");
udq_file	=fopen("data/udq_ref.txt", "r");

float ud_ref_MAT;
float uq_ref_MAT;
float h1;
float h2;
float h3;
float h4;
float f1;
float f2;
float thme;

for(int i = 0; i < HEIGHT; i++)
{
	printf("----------------------- \n");
	
	fscanf(H_file,"%f %f %f %f",&h1,&h2,&h3,&h4);
	fscanf(f_file,"%f %f",&f1,&f2);
	fscanf(thme_file,"%f",&thme);
	printf("thme = %.15f \n",thme);
	fscanf(udq_file,"%f %f",&ud_ref_MAT,&uq_ref_MAT);
	
	float H[2][2] = {
    {h1, h2},
    {h3, h4}
	};
	float f[2] = {f1, f2};
	float Sthme = sin(thme);
	float Cthme = cos(thme);
	
	aVsIs_dq3phase(&H,&f,&Ubus,&uOpt,&Sthme,&Cthme);
	
	printf("Ud_FILE= %.15f \n",ud_ref_MAT);
	printf("Uq_FILE= %.15f \n",uq_ref_MAT);
	
	printf("ud_OPT = %.15f \n", uOpt[0]);
	printf("uq_OPT = %.15f \n \n", uOpt[1]);
	
	printf("diff_ud = %f \n", ud_ref_MAT-uOpt[0]);
	printf("diff_uq = %f \n \n", uq_ref_MAT-uOpt[1]);
	
	if(abs(ud_ref_MAT-uOpt[0])>0.0001 || abs(uq_ref_MAT-uOpt[1])>0.0001){
		printf("Press Any Key to Continue\n");
		printf("index = %d \n", i);
		getchar(); 
	}
}
printf("Code works \n");
printf("\n");

fclose(f_file);
fclose(udq_file);
fclose(thme_file);
fclose(H_file);
}