/*      PROTOTYPES      */
float costFunctionEval(float (*H)[2][2], float (*f)[2], float (*u)[2]);
char checkConstraints(float (*vec)[2], float *uBus_sca);
void getCoefficientsSixSolutions(float (*H)[2][2], float (*f)[2], float *uBus_sca, float (*coeffSol_vec)[8]);
void checkSolution(float (*tempSol)[2], float (*tempOptSol)[2], float *costTempSol, float *costTempOptSol, int *flag, float (*H)[2][2], float (*f)[2]);

/*      LIBRARIES       */

#include <stdlib.h>

/*      CONSTANTS       */
#define     SQRT_3                                   ((float)1.732050807568877)
#define     ONE_OVER_SQRT_3                          ((float)0.577350269189626)
#define     TWO_OVER_SQRT_3                          ((float)1.154700538379252)
#define     TWO_OVER_THREE                           ((float)0.666666666666667)
#define     ONE_OVER_THREE                           ((float)0.333333333333333)


static const float A_constr[6][2]  = {
                                        {SQRT_3, 1.0},
                                        {0.0, 1.0},
                                        {-SQRT_3, 1.0},
                                        {-SQRT_3, -1.0},
                                        {0.0, -1.0},
                                        {SQRT_3, -1.0}
                                    };
static const float b_constr[6]     = {
                                        TWO_OVER_SQRT_3,
                                        ONE_OVER_SQRT_3,
                                        TWO_OVER_SQRT_3,
                                        TWO_OVER_SQRT_3,
                                        ONE_OVER_SQRT_3,
                                        TWO_OVER_SQRT_3
                                    };

static const float uAlphaVertex[6]  = {
                                        TWO_OVER_THREE,
                                        ONE_OVER_THREE,
                                        -ONE_OVER_THREE,
                                        -TWO_OVER_THREE,
                                        -ONE_OVER_THREE,
                                        ONE_OVER_THREE
                                    };
static const float uBetaVertex[6]   = {
                                        0.0,
                                        ONE_OVER_SQRT_3,
                                        ONE_OVER_SQRT_3,
                                        0.0,
                                        -ONE_OVER_SQRT_3,
                                        -ONE_OVER_SQRT_3
                                    };