/* 
April 22, 2017. Header file for srkVasculature 2D header file.
*/

// 2D defines and globals.
static char help[] = "Synthetic 2D/3D geometry for human vasculature, 1 process job.\n";

#include<stdbool.h>
/*
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscsys.h>
*/

/* my standard sundials headers, constants, and macros that will fit Petsc.
   Petsc cannot handle cvodes, at least not my installations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
// #include <mpi.h>

//#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.  */
//#include <nvector/nvector_serial.h>       /* serial N_Vector types, fcts., macros */
//#include <cvode/cvode_dense.h>            /* prototype for CVDense                */
//#include <sundials/sundials_dense.h>      /* definitions DlsMat DENSE_ELEM        */
//#include <sundials/sundials_types.h>      /* definition of type realtype          */
//#include <sundials/sundials_math.h>

#define debug 0 // switch all printf statements on or off. some may not be guarded by debug

//#define Ith(v,i)    		NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
//#define IJth(A,i,j) 		DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */
//#define RTOL        	RCONST(1.0e-12)   	  /* scalar relative tolerance            */
//#define ATOL        	RCONST(1.0e-6)        /* scalar absolute tolerance components */
//#define MAXSTEPS     5000
//#define ZERO        	RCONST(0.0)

// have physical coordinates here for the 2D model.
#define usr_MXP      100.0
#define usr_MYP      100.0
#define usr_MZP		0.0

#define billion 1000000000

/*========================================================================*/
/*
Dimensions of the 3D ventricle model. These dimensions are taken from Comput Mech (2010) 45: 227-243, page 240, Figure 7.
A margin of 5 on each side of the model.
*/
#define LV_X_RADIUS 30
#define LV_Y_RADIUS 30
#define LV_Z_RADIUS 70
#define LV_THICKNESS 12

#define RV_X_RADIUS 30
#define RV_Y_RADIUS 51
#define RV_Z_RADIUS 60
#define RV_THICKNESS 6

#define MARGIN 2.0 // mm

#define EPI_PHI        1.0
#define ENDO_PHI  -1.0

#define MAXORDERNUM 12

// when you do this sort of this, the LHS is always a double - not an int.
#define TOTAL_X (MARGIN + LV_X_RADIUS  + LV_X_RADIUS  + MARGIN)    // 2 times radius of LV plus margins at both ends.
#define TOTAL_Y (MARGIN + LV_Y_RADIUS + RV_Y_RADIUS + MARGIN)
#define TOTAL_Z (MARGIN + LV_Z_RADIUS + MARGIN)

/*
EPI is the outer surface of the ventricles, as well as the suface of the LV facing the RV in the septum. Septum is not segmented as of now.
ENDO is the inner surface of the ventricles. The Septum has inner surface on the LV side and outer surface on the RV side.
*/
#define EPI      1000
#define ENDO 2000

#define NBS    7            // 3D arrays have 7 neighbours: itself, and 6 surrounding it in a standard 1st order FD stencil.

#define EMPTY_SPACE -1

/*========================================================================*/

// taking 0.01 for DX makes the run longer than desirable.
#define DX          0.1   /* X internode spacing in mm. My model has the same dx in each direction. */

/*
#define usr_MX      ((usr_MXP)/(DX))
#define usr_MY      ((usr_MYP)/(DX))
#define usr_MZ		0
*/

// when you do this sort of this, the LHS is always a double - not an int.
#define usr_MX      (int)((TOTAL_X)/(DX))
#define usr_MY      (int)((TOTAL_Y)/(DX))
#define usr_MZ	   (int)((TOTAL_Z)/(DX))


// physical length and breadth of rectangle.
//#define usr_M		80.0 // mm, x direction
//#define usr_P		81.0 // mm, y direction

#define EMPTY_SPACE -1

#define DIM 2   /* Dimension, in 2D it is 2, in 3D it is 3. Do you want to play with nd problems? Maybe later. */

// 21 June 2017 revision: we now take the number of leaves and the nuNodes is determined based on how the tree work itself out.
// number of segments = (number of nodes - 1)
// number ofleaves = (number of nodes + 1)/2
#define	nuLeaves 6000
#define nuNodes  (2*nuLeaves) // total number of nodes in the model. This has to be at least 3. 1 is root, 2 is second.

#define MB_mm 	0.6413 // Units: kg / (seconds^3-mm).

#define Q_perf		4160.0 // mm^3/s, total perfusion at inlet.		
#define MURRAY	2.1
// #define MU_mm 	0.0000036 // // Pa-s in SI. When meters to mm is done, this is what we get.
#define Pin			100.0 * 133.3 // mm Hg.: Convert mmHg to Pa (=N/m^2) by multiplying by 133.3. 1 mmHg = 133.3 Pa/m^2
// #define Pout			20.0   // mm Hg.

#define qrequired 0.02 // per s in annealing.
// #define Q0 4160 // mm^3/s. in annealing.
#define			zeta 2.0

#define mu_0 0.0036 // hematocrit with units: grams / (mm-s).

/* min max macros. these are already in PetSc headers. */
// #define MAX(a,b) ((a) > (b) ? a : b)
// #define MIN(a,b) ((a) < (b) ? a : b)
	
#define TESTINGSEGS 	1 // 0 is for full run, 1 is when I was/am developing the exchange segments part of the code.

// the binary tree structure.
struct bin_tree {
int 		data;					// data is just an int identifier of the node.
int 		parent_data;			// to identify who is the parent.
int 		strahler;
int			parent_srahler;
double 	radius, length, perfusion, pressure, resistance, resistanceSegment, downResistance; // resistanceSegment is the resistance of me and my parent segment. downResistance is the resistance in me subtree.
double	x_cood, y_cood, z_cood; // this may also be in an array.
double 	theta1, theta2, theta;
	int assignedOrNot, onOff;
	int whichCase;
struct 		bin_tree *left, *right; // parent is also the data in the struct.
};

// my struct array to hold all this for calculation where the next node goes.
typedef struct {
	int me_data, parent_data, me_order;
	int assignedOrNot; // 0 1 array, 0 if not assigned, 1 if done and finished.
	int onOff; // this is for use when calculating coordinates.
	double length, x_cod, y_cod, z_cod, radius, resistance, resistanceSegment, pressure, perfusion, theta1, theta2, theta;
	int whichCase;
} nodeInfo;

typedef struct bin_tree node;
typedef struct { double x, y, z; } XYZ;
typedef struct {double total_perfusion, distanceToEndo; } arr;

#define NBS    7            // 3D arrays have 7 neighbours: itself, and 6 surrounding it in a standard 1st order FD stencil.

/**************************************************Random number generator.****************************************************************************/
/************************************************************************************************************************************************************************/
/*
Notes on the random number generator:
generates a random number on [0,1]-real-interval 
genrand_real1() gives you a double between 0 and 1. This is probably what you will use.
*/
/* For Random Numbers */
/* Period parameters */  
#define N 				624
#define M 				397
#define MATRIX_A 		0x9908b0dfUL   	/* constant vector a 		*/
#define UPPER_MASK 	0x80000000UL 	/* most significant w-r bits 	*/
#define LOWER_MASK 	0x7fffffffUL 		/* least significant r bits 	*/

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

int solution_broke = 0, somerandomnumber = 23334556;

/*
The random number generator.
*/

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */



/**************************** segment intersection function.***********************/
// this function is not used, I dont know if it works or not.
// int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,  double p2_x, double p2_y, double p3_x, double p3_y, double *i_x, double *i_y)
int get_line_intersection(double p0_x, double p0_y, double p1_x, double p1_y,  double p2_x, double p2_y, double p3_x, double p3_y){
    double s02_x, s02_y, s10_x, s10_y, s32_x, s32_y, s_numer, t_numer, denom, t;
    s10_x = p1_x - p0_x;
    s10_y = p1_y - p0_y;
    s32_x = p3_x - p2_x;
    s32_y = p3_y - p2_y;

    denom = s10_x * s32_y - s32_x * s10_y;
    if (fabs(denom)<0.001 ) 													        return 0; // Collinear

    bool denomPositive = denom > 0;

    s02_x = p0_x - p2_x;
    s02_y = p0_y - p2_y;
    s_numer = s10_x * s02_y - s10_y * s02_x;
    if ((s_numer < 0) == denomPositive)        											return 0; // No collision

    t_numer = s32_x * s02_y - s32_y * s02_x;
    if ((t_numer < 0) == denomPositive)        											return 0; // No collision
    if (((s_numer > denom) == denomPositive) || ((t_numer > denom) == denomPositive)) 	return 0; // No collision
    // Collision detected
    t = t_numer / denom;
/* actual intersection.
    if (i_x != NULL)
        *i_x = p0_x + (t * s10_x);
    if (i_y != NULL)
        *i_y = p0_y + (t * s10_y);
*/
    return 1;
}

/* *************************** segment intersection function.********************** */

/******************* VTK of segments. ***************************************************/
void writeRadiusVTK(double **coodsMe, int **srkNodes, double *radiusSegment, int accepted_cood){

	char str[1000];
	FILE *output;
	int temp_x;
	
	sprintf(str,"segmentsRadius%06d.vtk", accepted_cood);
	output = fopen(str,"w");
	fprintf(output,"# vtk DataFile Version 3.0\n");
	fprintf(output,"Segments\n");
	fprintf(output,"ASCII\n");
	fprintf(output,"\n");
	fprintf(output,"DATASET POLYDATA\n");
	fprintf(output,"POINTS %d float\n", nuNodes);
	for(temp_x = 0; temp_x < nuNodes; temp_x++)
	fprintf(output,"%f %f %f\n",coodsMe[temp_x][0], coodsMe[temp_x][1], coodsMe[temp_x][2]);
	fprintf(output,"\n");
	fprintf(output,"LINES %d %d\n", nuNodes-1, 3*(nuNodes-1) );
	for(temp_x = 0; temp_x < nuNodes-1; temp_x++)
	fprintf(output,"2 %d %d\n",srkNodes[temp_x][0], srkNodes[temp_x][1]);
	fprintf(output,"\n");
	fprintf(output, "CELL_DATA %d\n",nuNodes-1);
	fprintf(output, "SCALARS radius double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < nuNodes-1; temp_x++)
	fprintf(output,"%20.20f \n", radiusSegment[temp_x]);
	fclose(output);

}


/******************* VTK of segments. ***************************************************/

/* Cost function. Calculates total cost. */
double totalCost(double *radiusSegment, double *lengthSegment, double *perfusion, double *resistivity, double **coodsMe, int *leaves, double radius_V_mm ){

double CF_tissue_volume, CF_power, CF_tissue_supply, bCost, X1, X2, Y1, Y2, Z1, Z2, length, total_cost;
int temp_x, temp_y, mecoodsIdx, terminals;

// volume and power and spheres of influence.
CF_tissue_volume = 0.0; CF_power = 0.0;

// volume and power.
for(temp_x = 0; temp_x < nuNodes; temp_x++){
CF_tissue_volume 	= CF_tissue_volume 	+ MB_mm * radiusSegment[temp_x] * lengthSegment[temp_x];
CF_power 			= CF_power 			+ perfusion[temp_x] * perfusion[temp_x] * resistivity[temp_x];
}

// Tissue supply.
CF_tissue_supply = 0.0;
	for(temp_x=0; temp_x<usr_MX; temp_x++)
	for(temp_y=0; temp_y<usr_MY; temp_y++){
	bCost = 0; // this is just before the temp_z loop start.
	X2 = (double)temp_x * DX;    Y2 = (double)temp_y * DX;    Z2 = 0.0;
	for(terminals=0; terminals<nuLeaves; terminals++){
		mecoodsIdx 						= leaves[terminals];
		X1 = coodsMe[mecoodsIdx][0];     Y1 = coodsMe[mecoodsIdx][1]; Z1 = 0.0;
		length 							= (X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2);
		if(length<=radius_V_mm) bCost 	= bCost + 1.0;
	} // end of terminals loop.
	if(bCost>0.0)    bCost = (bCost - 1.0)*(bCost - 1.0);
	else		  bCost = 10.0;
	CF_tissue_supply      = CF_tissue_supply + bCost;
	} // end of X-Y loop.


total_cost = CF_tissue_supply * pow(10, 30) + (CF_tissue_volume + CF_power)  * pow(10, 4);

	return total_cost;
}
/* Cost function. Calculates total cost. */

/* Radius function ************************************************************************************************************************************/
// this function needs to become recursion.
void RadiusSegments(int ** srkNodes, int *ThisNodeIsALeaf, int **srkStrahlerNumberOfNode, double *radiusSegment, int *parent){

double ratio, radiusLeft, radiusRight, ratioTemp1, ratioTemp2;
int temp_x, temp_y, numberofChildren, locLeft, locRight, leftStrahler, rightStrahler, allRadiiDefined;
	
temp_x 					= 0;
radiusSegment[temp_x] 	= 2.1; // mm. temp_x = 0 is always the first/root segment. 

// look for children of srkNodes[0][0]
// identify the 2 children, put their data into locLeft and locRight.
numberofChildren = 0; locLeft = locRight = leftStrahler = rightStrahler = -1;
for(temp_y=0; temp_y < nuNodes; temp_y++)
		if(srkNodes[temp_x][0]==srkNodes[temp_y][1]){
			numberofChildren++;
			if(numberofChildren==1){
			locLeft 				= temp_y; 	leftStrahler 			= srkStrahlerNumberOfNode[temp_y][0];
			}
			if(numberofChildren==2){
			locRight 			= temp_y; 	rightStrahler 		= srkStrahlerNumberOfNode[temp_y][0];
			}		
if(numberofChildren==2) break;	
		}	

if(debug)
printf("Loc, Str (L, R), first: %d\t%d\t%d\t%d\n", locLeft, leftStrahler, locRight, rightStrahler);

// calculate the radii of locLeft and locRight.

ratio 				= (double)leftStrahler/(double)rightStrahler;
ratioTemp1			= pow(ratio, (double)MURRAY);
ratioTemp2 			= pow( (1.0 + ratioTemp1), 1.0/(double)MURRAY);

radiusSegment[locRight]	= radiusSegment[temp_x] / ratioTemp2;
radiusSegment[locLeft]	= radiusSegment[locRight] * ratio;

parent[temp_x]								= 2; // segment 0 is done, it has radius = 2.1
if(ThisNodeIsALeaf[locLeft]<0)  parent[locLeft] 	= 1; // this goes live.
if(ThisNodeIsALeaf[locRight]<0) parent[locRight] 	= 1; // this goes live also.

if(debug)
printf("Strs: %f %f %f\n", ratio, ratioTemp1, ratioTemp2);

// now look in parent array and come up with the first location that needs dealing with.
temp_x = -1;
for(temp_y=0; temp_y < nuNodes; temp_y++){ if(parent[temp_y]==1) temp_x = temp_y; if(temp_x>0) break; }

if(debug)
printf("next node at location in srkNodes: %d %d\n", temp_x, locLeft);

// do the next part.
// look for children of srkNodes[temp_x][0]
// identify the 2 children, put their data into locLeft and locRight.

do{

numberofChildren = 0; locLeft = locRight = leftStrahler = rightStrahler = -1;
for(temp_y=0; temp_y < nuNodes; temp_y++)
		if(srkNodes[temp_x][0]==srkNodes[temp_y][1]){
			numberofChildren++;
			if(numberofChildren==1){
			locLeft 				= temp_y; 	leftStrahler 			= srkStrahlerNumberOfNode[temp_y][0];
			}
			if(numberofChildren==2){
			locRight 			= temp_y; 	rightStrahler 		= srkStrahlerNumberOfNode[temp_y][0];
			}		
if(numberofChildren==2) break;	
		}	

if(debug)
printf("Loc, Str (L, R), second: %d\t%d\t%d\t%d\n", locLeft, leftStrahler, locRight, rightStrahler);

ratio 				= (double)leftStrahler/(double)rightStrahler;
ratioTemp1			= pow(ratio, (double)MURRAY);
ratioTemp2 			= pow( (1.0 + ratioTemp1), 1.0/(double)MURRAY);

radiusSegment[locRight]	= radiusSegment[temp_x] / ratioTemp2;
radiusSegment[locLeft]	= radiusSegment[locRight] * ratio;

parent[temp_x]								= 2; // segment 0 is done, it has radius = 2.1
if(ThisNodeIsALeaf[locLeft]<0)  parent[locLeft] 	= 1; // this goes live.
if(ThisNodeIsALeaf[locRight]<0) parent[locRight] 	= 1; // this goes live also.

if(debug)
printf("Strs: %f %f %f\n", ratio, ratioTemp1, ratioTemp2);

// now look in parent array and come up with the first location that needs dealing with.
temp_x = -1;
for(temp_y=0; temp_y < nuNodes; temp_y++){ if(parent[temp_y]==1) temp_x = temp_y; if(temp_x>0) break; }

if(debug)
printf("LOOP, next node at location in srkNodes: %d %d\n", temp_x, locLeft);

allRadiiDefined = 0;
for(temp_y=0; temp_y < nuNodes-1; temp_y++){ if(radiusSegment[temp_y]>=0) allRadiiDefined++; }

} while(allRadiiDefined<nuNodes-1);

} // end of radius function.

/* end radius iterative function ****************************************************************************************************************************/

/* Bounding box function ************************************************************************************************************************************/
/* This function only does internal nodes, not leaves. */
void boundingBoxx(int **srkNodes, double **boundingBox, double **coodsMe, int *leaves){

	int temp_x, temp_y, locLeft, locRight, terminals, cont, cn;
	double x_min, x_max, y_min, y_max, z_min, z_max;

for(temp_x=0; temp_x < nuNodes; temp_x++){
	y_min = x_min =  100000000.0; y_max = x_max = -100000000.0; z_min = z_max = 0.0; cont = 0;
	for(terminals=0; terminals<nuLeaves;terminals++) if(leaves[terminals]==srkNodes[temp_x][0]) cont = 1;
	if(cont==1) 					continue;
	if(srkNodes[temp_x][0]==0) 	continue;

	// do me.
	cn = srkNodes[temp_x][0];
	if(x_min > coodsMe[cn][0]) x_min = coodsMe[cn][0];
	if(y_min > coodsMe[cn][1]) y_min = coodsMe[cn][1];
	if(x_max < coodsMe[cn][0]) x_max = coodsMe[cn][0];
	if(y_max < coodsMe[cn][1]) y_max = coodsMe[cn][1];

	// do parent.
	cn = srkNodes[temp_x][1];
	if(x_min > coodsMe[cn][0]) x_min = coodsMe[cn][0];
	if(y_min > coodsMe[cn][1]) y_min = coodsMe[cn][1];
	if(x_max < coodsMe[cn][0]) x_max = coodsMe[cn][0];
	if(y_max < coodsMe[cn][1]) y_max = coodsMe[cn][1];

	locLeft = locRight = -1;
	for(temp_y=0; temp_y<nuNodes; temp_y++){
		if(srkNodes[temp_x][0]==srkNodes[temp_y][1]){
		if(locLeft>=0&& locRight<0) locRight = srkNodes[temp_y][0];
		if(locLeft<0 && locRight<0) locLeft  = srkNodes[temp_y][0];
		if(locLeft>=0&&locRight>=0) break;
		}
	}

	// do child left
	cn = locLeft;
	if(x_min > coodsMe[cn][0]) x_min = coodsMe[cn][0];
	if(y_min > coodsMe[cn][1]) y_min = coodsMe[cn][1];
	if(x_max < coodsMe[cn][0]) x_max = coodsMe[cn][0];
	if(y_max < coodsMe[cn][1]) y_max = coodsMe[cn][1];

	// do child right
	cn = locRight;
	if(x_min > coodsMe[cn][0]) x_min = coodsMe[cn][0];
	if(y_min > coodsMe[cn][1]) y_min = coodsMe[cn][1];
	if(x_max < coodsMe[cn][0]) x_max = coodsMe[cn][0];
	if(y_max < coodsMe[cn][1]) y_max = coodsMe[cn][1];

	boundingBox[temp_x][0] = x_min; boundingBox[temp_x][1] = y_min; 
	boundingBox[temp_x][2] = x_max; boundingBox[temp_x][3] = y_max; 
} // end of boundingBox loop.


} // end of boundingBoxx


/************************************************************************************************************************************************************************/
/************************************************************************************************************************************************************************/

// my srk recursion functions.

// new constrct tree functions from:http://www.geeksforgeeks.org/construct-a-binary-tree-from-parent-array-representation/
/*
other sources of same code:
http://www.techiedelight.com/build-binary-tree-given-parent-array/
http://www.ideserve.co.in/learn/construct-binary-tree-from-parent-array
https://stackoverflow.com/questions/37941318/how-to-build-an-incomplete-binary-tree-from-array-representation
*/

// Utility function to create new Node
node *newNode(int data){
    node *temp = (node*)malloc(sizeof(node));
    temp->data  = data;
    temp->left  = temp->right = NULL;
    return (temp);
}

// Creates a node with data as 'i'.  If i is root, then it changes
// root.  If parent of i is not created, then it creates parent first
void createNode(int *parent, int i, node **created, node **root){
    // If this node is already created
    if (created[i] != NULL)         return;
    // Create a new node and set created[i]
    created[i] = newNode(i);
    // If 'i' is root, change root pointer and return
    if (parent[i] == -1){ *root = created[i]; return; }
 
    // If parent is not created, then create parent first
    if (created[parent[i]] == NULL) createNode(parent, parent[i], created, root);
 
    // Find parent pointer
    node *p = created[parent[i]];
 
    // If this is first child of parent
    if (p->left == NULL) 	p->left  = created[i];
    else 			        p->right = created[i]; // If second child
} // end of createNode.



// Creates tree from parent[0..n-1] and returns root of the created tree
node *createTree(int *parent) {
	int i;
    // Create an array created[] to keep track
    // of created nodes, initialize all entries
    // as NULL
    node *created[nuNodes]; // array of pointers of size nuNodes.
    for (i=0; i<nuNodes; i++)      created[i] = NULL;
 
    node *root = NULL;
    for (i=0; i<nuNodes; i++)     createNode(parent, i, created, &root);
 
    return root;
}


/* If target is present in tree, then prints the ancestors and returns true, otherwise returns false. */
bool printAncestors(node *root, int target, int *TheParent) {
  /* base cases */
  if (root == NULL) 			     return false;   if (root->data == target)	     return true;
 
  /* If target is present in either left or right subtree of this node,
     then print this node */
  if ( printAncestors(root->left, target, TheParent) || printAncestors(root->right, target, TheParent) )  {
//	  printf("Ancestor: %d \n", root->data);
	  if( *TheParent < 0) *TheParent = root->data;
	  return true;
  }
  /* Else return false */
  return false;
}
 
/* Insert function for the binary tree. */
void insert(node ** tree, int val)
{
    node *temp = NULL;
    if(!(*tree)){ // this is for any new node without children.
        temp = (node *)malloc(sizeof(node));
        temp->data = val;
        temp->left = temp->right = NULL;
        *tree = temp;
        return;
    }
         if(val < (*tree)->data)    insert(&(*tree)->left,  val);
    else if(val > (*tree)->data) 	insert(&(*tree)->right, val);
} // end of insert function. This is a void function.

void print_preorder(node * tree)
{
    if (tree){
        printf("%d ",tree->data);
//        tree->data					; // whats this doing?
        print_preorder(tree->left)		;
        print_preorder(tree->right)	;
    }
}

/* Compute the "height" of a tree -- the number of
    nodes along the longest path from the root node
    down to the farthest leaf node.*/
int height(node *node){
    if (node==NULL)        return 0;
    else{
        /* compute the height of each subtree */
        int lheight = height(node->left);
        int rheight = height(node->right);
        /* use the larger one */
        if (lheight > rheight)
            return(lheight+1);
        else return(rheight+1);
    }
}

/* Print nodes at a given level, write to a file so that I can read it. */
void printGivenLevel2(node *tree, int level){
	FILE *print_level;

    if (tree == NULL)         return;
    if (level == 1){
	print_level = fopen("level.txt","a");	    
//	    printf("%d ", tree->data);
	    fprintf(print_level, "%d ", tree->data);
	    fclose(print_level);	    
    }
    else if (level > 1){
				        printGivenLevel2(tree->left, level-1);
				        printGivenLevel2(tree->right, level-1);
    }
}



/* Print nodes at a given level */
void printGivenLevel(node *tree, int level){
    if (tree == NULL)         return;
    if (level == 1) 	        printf("(Radius: %f, length: %f\n ", tree->radius, tree->length);
    else if (level > 1){
				        printGivenLevel(tree->left, level-1);
				        printGivenLevel(tree->right, level-1);
    }
}

/* Function to print level order traversal a tree*/
void printLevelOrder(node *tree, int val){ // val at 1st call is not used.
    int h = height(tree);
    int i;
    for (i=1; i<=h; i++){
        printGivenLevel(tree, i);
        printf("\n");
    }
}

//get the leaves count
unsigned int getLeafCount(node *node){
 if(node == NULL) 									  	return 0;
 if(node->left == NULL && node->right == NULL)  			return 1;
 else   												return getLeafCount(node->left)+getLeafCount(node->right);
}

/*  This function traverses tree in post order to 
    to delete each and every node of the tree */
void deleteTree(node *node){
    if (node == NULL) return;
 
    /* first delete both subtrees */
    deleteTree(node->left);
    deleteTree(node->right);
   
    /* then delete the node */
//    printf("\n Deleting node: %d", node->data);
    free(node);
    node = NULL;
} 

/* Helper function for getLevel().  It returns level of the data if data is
   present in tree, otherwise returns 0.*/
int getLevelUtil(node *node, int data, int level)
{
    if (node == NULL)
        return 0;
 
    if (node->data == data)
        return level;
 
    int downlevel = getLevelUtil(node->left, data, level+1);
    if (downlevel != 0)
        return downlevel;
 
    downlevel = getLevelUtil(node->right, data, level+1);
    return downlevel;
}
 
/* Returns level of given data value */
int getLevel(node *node, int data)
{
    return getLevelUtil(node,data,1);
}

// the next 2 functions are untested. Take from:
// http://www.techcrashcourse.com/2016/06/c-program-to-create-duplicate-binary-tree.html?m=1
//
node* getNewNode(int data) {
  /* dynamically allocate memory for a new node */
  node* newNode = (node*)malloc(sizeof(node));
  
  /* populate data in new Node */
  newNode->data = data;
  newNode->left = NULL;
  newNode->right = NULL;
   
  return newNode;
}
 
/* Returns a tree which is exact copy of passed tree */
node *cloneBinaryTree(node *root){
    if(root == NULL)
        return NULL;
    /* create a copy of root node */
    node *newNode = getNewNode(root->data);
    /* Recursively create clone of left and right sub tree */
    newNode->left = cloneBinaryTree(root->left);
    newNode->right = cloneBinaryTree(root->right);
    /* Return root of cloned tree */
    return newNode;
}

    /*_________________________________________________________________________________
____________________________________________________________________________________*/
    
    /*Function which helps the print_path to recursively print all the nodes*/ 
    void print_paths_recur(node *node, int path[], int path_len, int n1, int n2){ // this needs to return a value to print paths.
      FILE *parentTest;
      int i, count, return_value;      
      return_value = -1;

      if (node == NULL){
	      return; 
      }
      path[path_len] = node->data;
      path_len++;
      if (node->left == NULL && node->right == NULL) {
// this is where you test if n1 and n2 are in the same path or not.	      
//      for (i = 0; i < path_len; i++)   printf("%d -> ", path[i]); printf("\n");    
	      count = 0;
	      for (i = 0; i < path_len; i++){ if(path[i]==n1) count++; if(path[i]==n2) count++; }
              if(count > return_value) return_value = count;
              parentTest = fopen("p.txt","a+"); // temporary file.
              fprintf(parentTest, "%d\n", return_value);
              fclose(parentTest);
      }
      else{
        print_paths_recur(node->left,    path, path_len, n1, n2);    //recursively calls the left node of the tree
        print_paths_recur(node->right, path, path_len, n1, n2);    //recursively calls the right node of the tree
      }
    } // end of print_paths_recur.

    /*Function to store all the paths from the root node to all leaf nodes in  a array*/
    void print_paths(node *node, int n1, int n2){  // this needs to return a value to the main. 0 if n1 and n2 are parents
      int path[100000];
      print_paths_recur(node, path, 0, n1, n2);
    }


/* Function to be used in assignment of Strahler number to each node. This is the same as printPostorder, but changed for Strahler assignment. */
/* This function writes a file with node data and the Strahler number. */
int binaryStrahler(int NN, node *node){ // the N is the node data value. the *node is always root.

FILE *strahler;
int binStr, binStrL, binStrR;

     if (node == NULL)         return -1; // error trapping.
		if(node->left==NULL && node->right==NULL){ 	
//			printf("Leaf: %d Strahler: %d\n", node->data, 1); // leaf strahler. 
			binStr = 1;
			strahler = fopen("strahler.txt", "a+");
			fprintf(strahler, "%d %d %d\n", node->data, 1, 1); // data, strahler, ifnode is a leaf.
			fclose(strahler);
			node->strahler = binStr;
			return           binStr;
		}
		else if(node->left==NULL && node->right!=NULL){ 	
//			printf("%d Strahler: %d\n", node->data,  binaryStrahler(node->data, node->right) ); // at a node which has only 1 child.
			binStr = binaryStrahler(node->data, node->right);
			strahler = fopen("strahler.txt", "a+");
			fprintf(strahler, "%d %d %d\n", node->data, binStr, -1 );
			fclose(strahler);
			node->strahler = binStr;
			return           binStr;
		}
		else if(node->left!=NULL && node->right==NULL){ 	
//			printf("%d Strahler: %d\n", node->data,  binaryStrahler(node->data, node->left) ); // at a node which has only 1 child.
			binStr = binaryStrahler(node->data, node->left);
			strahler = fopen("strahler.txt", "a+");
			fprintf(strahler, "%d %d %d\n", node->data, binStr, -1);
			fclose(strahler);
			node->strahler = binStr;
			return           binStr;
		}
		else if(node->left!=NULL && node->right!=NULL){
			binStrL = binaryStrahler(node->data, node->left);
			binStrR = binaryStrahler(node->data, node->right);

			if(binStrL==binStrR) binStr = binStrL + 1; // it could also be binStrR + 1, which is the same.
			if(binStrL< binStrR) binStr = binStrR;
			if(binStrL> binStrR) binStr = binStrL;

			strahler = fopen("strahler.txt", "a+");
			fprintf(strahler, "%d %d %d\n", node->data, binStr, -1);
			fclose(strahler);
			node->strahler = binStr;
			return           binStr;
		}


return -10000; // this return value has exit meaning.

} // end of binaryStrahler.


// for the 3D model.
double binaryDownResistance3D(node *node){

double rR, rL;
	
if(node==NULL){ printf("got a null node in down resistance of 3D, exiting.\n"); exit(0); } // error trapping, must never happen.

	 if(node->left==NULL && node->right==NULL){ // this is a leaf, base case.
		return node->downResistance = 0.0;
	}
else if(node->left!=NULL && node->right==NULL){
//		printf("left, node->right %d", node->left->data);
		return node->downResistance = node->left->resistanceSegment + binaryDownResistance3D(node->left);
	}
else if(node->left==NULL && node->right!=NULL){
//	printf("right, node->right %d", node->right->data);
		return node->downResistance = node->right->resistanceSegment + binaryDownResistance3D(node->right);
	}
else if(node->left!=NULL && node->right!=NULL){
	rR = node->right->resistanceSegment + binaryDownResistance3D(node->right  );
	rL = node->left   ->resistanceSegment + binaryDownResistance3D(node->left     );
//	printf("two sides: %20.20f\n", rR * rL / (rR + rL));
		return node->downResistance = rR * rL / (rR + rL);
	}
	return -1000000; // this return value has exit meaning, must not happen.
} // end of binaryDownResistance


// this binarydownResistance is working.
/* definition of binarydownResistance();*/
//  NN is the child node to which resistivty is attached to. Node is pointer to this child node.
// I am not passing NN any more, dont need it.
double binarydownResistance(node *node, double *resistivity, int **srkNodes){
	
	FILE *resist;
	
	int tmp_x, tmpx_hold, tmp_xL, tmpxL_hold, tmp_xR, tmpxR_hold;
	double rL, rR;
// resistanceSegment, downResistance	
	if(node==NULL) return -1.0; // error trapping, this must never ever happen.

	if(node->left==NULL && node->right==NULL){ // this is a leaf.
	// at leaves, the total downstream resistance is 0.
		node->downResistance = 0.0;
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++) if(node->data==srkNodes[tmp_x][0]){ tmpx_hold = tmp_x; }
		node->resistanceSegment = resistivity[tmpx_hold];
		resist = fopen("resistance.txt","a+"); fprintf(resist, "%d %f\n", node->data, node->downResistance); fclose(resist);
	return node->downResistance;
	}
	else if(node->left!=NULL && node->right==NULL){ // only left node. In that case, look for the resistivity at node->left->data and add that to binaryDownR
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++) if(node->left->data==srkNodes[tmp_x][0]){ tmpx_hold = tmp_x; }
		node->downResistance = resistivity[tmpx_hold] + binarydownResistance(node->left, resistivity, srkNodes);
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++) if(node->data==srkNodes[tmp_x][0]){ tmpx_hold = tmp_x; }
		node->resistanceSegment = resistivity[tmpx_hold];
		resist = fopen("resistance.txt","a+");fprintf(resist, "%d %f\n", node->data, node->downResistance); fclose(resist);
		return node->downResistance;		
	}
	else if(node->left==NULL && node->right!=NULL){ // only right node. In that case, look for the resistivity at node->left->data and add that to binaryDownR
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++) if(node->right->data==srkNodes[tmp_x][0]){ tmpx_hold = tmp_x; }
		node->downResistance = resistivity[tmpx_hold] + binarydownResistance(node->right, resistivity, srkNodes);
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++) if(node->data==srkNodes[tmp_x][0]){ tmpx_hold = tmp_x; }
		node->resistanceSegment = resistivity[tmpx_hold];
		resist = fopen("resistance.txt","a+");fprintf(resist, "%d %f\n", node->data, node->downResistance); fclose(resist);
	return node->downResistance;		
	}
	else if(node->left!=NULL && node->right!=NULL){ 
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++){
			 if(node->left->data==srkNodes[tmp_x][0]){ tmp_xL = tmp_x; }
			 if(node->right->data==srkNodes[tmp_x][0]){ tmp_xR = tmp_x; }
		}
		rL =  resistivity[tmp_xL]  + binarydownResistance(node->left,  resistivity, srkNodes);
		rR =  resistivity[tmp_xR]  + binarydownResistance(node->right, resistivity, srkNodes);
		node->downResistance = rL * rR / (rL + rR);
		
		tmpx_hold = -1;
		for(tmp_x=0; tmp_x < nuNodes; tmp_x++) if(node->data==srkNodes[tmp_x][0]){ tmpx_hold = tmp_x; }
/*		if(tmpx_hold > 0) */	node->resistanceSegment = resistivity[tmpx_hold];
//		else					node->resistanceSegment = 0.0;			
		resist = fopen("resistance.txt","a+"); fprintf(resist, "%d %f\n", node->data, node->downResistance); fclose(resist);
		return node->downResistance;
	}

	return -1000000; // this return value has exit meaning.
} // end binaryResistance.
    
void writeFile(node *tNode){

FILE *nodeNumbas;	
	
	  if(tNode==NULL) return; // error trapping, this must never ever happen.	

	
	      if(tNode->left!=NULL&&tNode->right!=NULL){
	nodeNumbas = fopen("meandParent.txt", "a+");
	fprintf(nodeNumbas, "%d %d %d %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %d %10.10f %10.10f %10.10f %d %d %d\n", tNode->data, tNode->parent_data, tNode->strahler, tNode->length,  tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->theta1, tNode->theta2, tNode->theta, tNode->assignedOrNot, tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->onOff, tNode->assignedOrNot, tNode->whichCase);
	fclose(nodeNumbas);
			writeFile(tNode->left); writeFile(tNode->right);
		}
else	if(tNode->left!=NULL&&tNode->right==NULL){
	nodeNumbas = fopen("meandParent.txt", "a+");
	fprintf(nodeNumbas, "%d %d %d %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %d %10.10f %10.10f %10.10f %d %d %d\n", tNode->data, tNode->parent_data, tNode->strahler, tNode->length,  tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->theta1, tNode->theta2, tNode->theta, tNode->assignedOrNot, tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->onOff, tNode->assignedOrNot, tNode->whichCase);
	fclose(nodeNumbas);	
			writeFile(tNode->left);
		}
else	if(tNode->left==NULL&&tNode->right!=NULL){
	nodeNumbas = fopen("meandParent.txt", "a+");
	fprintf(nodeNumbas, "%d %d %d %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %d %10.10f %10.10f %10.10f %d %d %d\n", tNode->data, tNode->parent_data, tNode->strahler, tNode->length,  tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->theta1, tNode->theta2, tNode->theta, tNode->assignedOrNot, tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->onOff, tNode->assignedOrNot, tNode->whichCase);
	fclose(nodeNumbas);	
			writeFile(tNode->right);			
		}		
else	if(tNode->left==NULL&&tNode->right==NULL){ 		
	nodeNumbas = fopen("meandParent.txt", "a+");
	fprintf(nodeNumbas, "%d %d %d %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %10.10f %d %10.10f %10.10f %10.10f %d %d %d\n", tNode->data, tNode->parent_data, tNode->strahler, tNode->length,  tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->theta1, tNode->theta2, tNode->theta, tNode->assignedOrNot, tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->onOff, tNode->assignedOrNot, tNode->whichCase);
	fclose(nodeNumbas);	
		}
}

// recursion used in the 3D code.
void perfusionPressure3D(node *node, double Pout){

double pme, pin, pout, qin, qL, qR, pL, pR, dr, rme, rL, rR, drL, drR, x, x1, x2, x_norm, dRd;
int tmp_x, tmpx_hold, temp_x, NNL, NNR;

// These are variables for calculating the bifurcation angles theta1 and theta2.
// By convention, d1 > d2, which makes theta1 < theta2.
double d1, d2, q1, q2, d0, q0;
double theta1, theta2;
double term1, term2, term3, term4;
double temp_p;
	
	  if(node==NULL) return; 
	  
else  if(node->left==NULL&&node->right==NULL){
if(debug)
 printf("Pressure: %f %f\n", node->pressure, node->perfusion);
	 return;
	}
else{
		pin 	= node->pressure;
		qin 	= node->perfusion;
		dRd = node->downResistance;
		pout 	= Pout;
		// temp init: error trapping.
		pL = pR = -1.0;
		qL = qR = -1.0;
		temp_p = pin - qin * dRd;

if(node->left!=NULL&&node->right!=NULL){
			drL 		= node->left ->downResistance;
			drR 	= node->right->downResistance;
			rL 		= node->left ->resistanceSegment;
			rR 		= node->right->resistanceSegment;

			// perfusion and pressure of next nodes.
			qL = qin * (rR+drR)/(rL+drL) / (1.0 + (rR+drR)/(rL+drL)  ); 
			qR = qin * (rL+drL)/(rR+drR) / (1.0 + (rL+drL)/(rR+drR)  ); 

			pL = pin - qL * (rL);  // poissule at left node for pressure at left node.
			pR = pin - qR * (rR);
 if(debug)
printf("Perfusions: %f %f %f %f %f %f\n", qL, qR, qL+qR, qin, pL, pR);

			// error trapping?
			if(fabs(qL+qR - qin) > 0.1) printf("conservation not, see perfusionPressure3D in the header file, bifurcation %f %f %f!!\n", qL, qR, qin);
			// assign.
			node->left->pressure   = pL;
			node->left->perfusion  = qL;
			
			node->right->pressure  = pR;
			node->right->perfusion = qR;
			
//			if( fabs(pL - Pout) < 0.1 || fabs(pR - Pout) < 0.1){ printf("The pressure is wrong, tree goes both ways.\n"); exit(0); }
			
// this is where you assign the theta1 and theta2.
// there are definitions by: a) Yang and Wang; b) Fung (a la Kassab), Nic Smith; and c) Tamaddon et al.
// I am using the one based on approximate calculation of perfusion. Nic Smith et al, and Yang and Wang.

q0 = node->perfusion;   d0 = 2.0*(node->radius);
term1 = q0 * q0 / (d0 * d0 * d0 * d0);

if( (node->left->radius) > (node->right->radius) ){
	d1 = 2.0*(node->left->radius)		;
	q1 = (node->left->perfusion)			;
	d2 = 2.0*(node->right->radius)		;
	q2 = (node->right->perfusion)		;
}
else
{
	d2 = 2.0*(node->left->radius)		;
	q2 = (node->left->perfusion)			;
	d1 = 2.0*(node->right->radius)		;
	q1 = (node->right->perfusion)		;
}

term2 = q1 * q1 / (d1 * d1 * d1 * d1)			;
term3 = q2 * q2 / (d2 * d2 * d2 * d2)		;
term4 = 2.0 * q0 * q1 / (d0 * d0 * d1 *d1)	;
x1 = (term1 + term2 - term3)/term4		;
x2 = (term1 - term2 + term3)/term4		;

// now normalise the two.
x_norm = sqrt(x1 * x1 + x2 * x2)			;
x1 = x1/x_norm						;
x2 = x2/x_norm						;

	theta1 =   acos(x1 );
	theta2 =  -acos(x2);


if(theta1!=theta1) {printf("theta1 is not defined. back to the recursion.\n") ; exit(0); }
if(theta2!=theta2) {printf("theta2 is not defined. back to the recursion.\n") ; exit(0); }

node->theta1 = theta1;
node->theta2 = theta2;

node->left->theta = theta1;
node->right->theta = theta2;

			perfusionPressure3D(node->left, Pout);
			perfusionPressure3D(node->right, Pout);

			if(debug)
			printf("pressure, perfusion: %f %f\n", node->pressure, node->perfusion);
		}
	else if(node->left!=NULL&&node->right==NULL){

			drL 		= node->left ->downResistance;
			rL 		= node->left ->resistanceSegment;
		
			qL = qin;			
			pL = pin - qL * (rL);  // poissule at left node for pressure at left node.

			// error trapping?
			if(fabs(qL - qin) > 0.1) printf("conservation not, see perfusionPressure3D in the header file, left!! %f %f %f!!\n", qL, qR, qin);
			// assign.
			node->left->pressure   = pL;
			node->left->perfusion  = qL;			

//			if( fabs(pL - Pout) < 0.1 ){ printf("The pressure is wrong, Left %f %50.50f %50.50f %f %f %f\n%50.50f\n", pin, pL, Pout, rL, qL, rL * qL, pL - Pout); exit(0); }

			perfusionPressure3D(node->left, Pout);
//			perfusionPressure3D(node->right, Pout);

// printf("%f %f %f %f\n", node->perfusion, node->radius, node->left->perfusion, node->left->radius );

// this is where you assign the theta1 and theta2.
theta1 = theta2 = 0.0; // theta2 is never used.
theta1 = 37.0 / 360.0 * (2.0 * M_PI);
node->theta1 = theta1;
node->theta2 = theta2;

			if(debug)
			printf("pressure, perfusion: %f %f\n", node->pressure, node->perfusion);
		}
	else if(node->left==NULL&&node->right!=NULL){
			rR 		= node->right->resistanceSegment;
			drR 	= node->right->downResistance;

			qR = qin;
			pR = pin - qR * (rR);
			// error trapping?
//			if(fabs(qR - qin) > 0.01) printf("conservation not, see perfusionPressure3D in the header file, right!! %f %f %f!!\n", qL, qR, qin);
			// assign.
//			node->left->pressure   = pL;
//			node->left->perfusion  = qL;
			node->right->pressure  = pR;			
			node->right->perfusion = qR;
			
//			if(fabs(pR - Pout) < 0.1){ printf("The pressure is wrong, Right %f %f\n", pR, Pout); exit(0); }			
			
//			perfusionPressure3D(node->left, Pout);
			perfusionPressure3D(node->right, Pout);

// this is where you assign the theta1 and theta2.		RIGHT NEVER HAPPENS if there is no bifurcation!
theta1 = theta2 = 0.0; // both theta1 and theta2 are never used.
node->theta1 = theta1;
node->theta2 = theta2;

			if(debug)			
			printf("pressure, perfusion: %f %f\n", node->pressure, node->perfusion);
		}	
} // end of big else.

} // end of perfusionPressure3D.


/*
This recursion is used in the 2D vasculature I was making for simulated annealing. For the 3D, I may have written another function. See the function
call in the main code.
*/
void perfusionPressure(node *node, double *resistivity, int **srkNodes, double *perfusion, double *pressure, double Pout){

double pme, pin, pout, qin, qL, qR, pL, pR, dr, rme, rL, rR, drL, drR, x;
int tmp_x, tmpx_hold, temp_x, NNL, NNR;

if(node==NULL) return; // error trapping, this must never ever happen.
if(node->left!=NULL&&node->right!=NULL){

	pin 	= node->pressure;
	qin 	= node->perfusion;
	pout 	= Pout;
	rL 		= node->left ->resistanceSegment;
	rR 		= node->right->resistanceSegment;
	drL 		= node->left ->downResistance;
	drR 	= node->right->downResistance;
	// perfusion and pressure of next nodes.
	pL = pR = -1.0;
	qL = qR = -1.0;
	// perfusion first.
	x = (rR + drR)/(rL + drL);
	qR = qin/(1.0 + x);
	qL = qin - qR;
	pL = pin - qL * rL;
	pR = pin - qR * rR;
	// assign.
	node->left->pressure   = pL;
	node->right->pressure  = pR;
	node->left->perfusion  = qL;
	node->right->perfusion = qR;
	NNL = NNR = -1;
	NNL = node->left ->data;
	NNR = node->right->data;
	for(temp_x = 0; temp_x < nuNodes; temp_x++){
		if(srkNodes[temp_x][0]==NNL){
			pressure[temp_x] = pL; perfusion[temp_x] = qL;
		}
		if(srkNodes[temp_x][0]==NNR){
			pressure[temp_x] = pR; perfusion[temp_x] = qR;
		}
	}
// printf("%f %f %f %f %f\n", qin, qL, qR, pL, pR);
	// call left and right.
perfusionPressure(node->left,  resistivity, srkNodes, perfusion, pressure, Pout);
perfusionPressure(node->right, resistivity, srkNodes, perfusion, pressure, Pout);
}
else if(node->left==NULL&&node->right!=NULL){

	pin 	= node->pressure;
	qin 	= node->perfusion;
	pout 	= Pout;
//	rL 		= node->left ->resistanceSegment;
	rR 		= node->right->resistanceSegment;
//	drL 	= node->left ->downResistance;
	drR 	= node->right->downResistance;
	// perfusion and pressure of next nodes.
	pL = pR = -1.0;
	qL = qR = -1.0;
	// perfusion first.
//	qL = (pin - pout)/(rL + drL);
	qR = qin; // (pin - pout)/(rR + drR);
//	pL = pin - qL * rL;
	pR = pin - qR * rR;
	// assign.
//	node->left->pressure   = pL;
	node->right->pressure  = pR;
//	node->left->perfusion  = qL;
	node->right->perfusion = qR;
	NNL = NNR = -1;
//	NNL = node->left ->data;
	NNR = node->right->data;
	for(temp_x = 0; temp_x < nuNodes; temp_x++){
		if(srkNodes[temp_x][0]==NNL){
			pressure[temp_x] = pL; perfusion[temp_x] = qL;
		}
		if(srkNodes[temp_x][0]==NNR){
			pressure[temp_x] = pR; perfusion[temp_x] = qR;
		}
	}
// printf("%f %f %f %f %f\n", qin, qL, qR, pL, pR);
	// call right.
perfusionPressure(node->right, resistivity, srkNodes, perfusion, pressure, Pout);
}
else if(node->left!=NULL&&node->right==NULL){

	pin 	= node->pressure;
	qin 	= node->perfusion;
	pout 	= Pout;
	rL 		= node->left ->resistanceSegment;
//	rR 		= node->right->resistanceSegment;
	drL 	= node->left ->downResistance;
//	drR 	= node->right->downResistance;
	// perfusion and pressure of next nodes.
	pL = pR = -1.0;
	qL = qR = -1.0;
	// perfusion first.
	qL = qin; // (pin - pout)/(rL + drL);
//	qR = (pin - pout)/(rR + drR);
	pL = pin - qL * rL;
//	pR = pin - qR * rR;
	// assign.
	node->left->pressure   = pL;
//	node->right->pressure  = pR;
	node->left->perfusion  = qL;
//	node->right->perfusion = qR;
	NNL = NNR = -1;
	NNL = node->left ->data;
//	NNR = node->right->data;
if(debug)
printf("NNL is %d %f %f\n", NNL, pL, qL);
	for(temp_x = 0; temp_x < nuNodes; temp_x++){
		if(srkNodes[temp_x][0]==NNL){
			pressure[temp_x] = pL; perfusion[temp_x] = qL;
		}
		if(srkNodes[temp_x][0]==NNR){
			pressure[temp_x] = pR; perfusion[temp_x] = qR;
		}
	}
// printf("In vasc: %f %f %f %f %f\n", qin, qL, qR, pL, pR);
	// call left.
perfusionPressure(node->left,  resistivity, srkNodes, perfusion, pressure, Pout);
}
else if(node->left==NULL&&node->right==NULL){
// do nothing: by now, it is all done.
return;
}

} // end of perfusionPressure


/*_________________________________________________________________________*/
// now work out the pressures and perfusions, given pressure and perfusion at inlet, and pressures at the leaves.
// dont use this function, see above.
/*_________________________________________________________________________*/

// return value is pressure.
void pressurePerfusion(node *node){

	FILE *pressure, *perfusion;
	
	double finalPressure = 2.0; // mmHg
	double deltaP;
	double perfL, perfR;
	double rL, rR;
	
	if(node==NULL) return; // error trapping, this must never ever happen.

	if(node->left!=NULL && node->right!=NULL){ 
		node->pressure = finalPressure + (node->perfusion)*(node->downResistance); // mmHg, given		
		perfL  = (node->pressure-finalPressure)/( node->left->resistanceSegment + node->left->downResistance );
		perfR = (node->pressure-finalPressure)/( node->right->resistanceSegment + node->right->downResistance );		
		node->left->perfusion = perfL;
		node->right->perfusion = perfR;
		pressurePerfusion(node->left);
		pressurePerfusion(node->right);	

		pressure = fopen("pressure.txt","a+");
		fprintf(pressure, "%d %f\n", node->data, node->pressure);
		fclose(pressure);				
		perfusion = fopen("perfusion.txt","a+");		
		fprintf(perfusion, "%d %f\n", node->data, node->perfusion);
		fclose(perfusion);
		
	if(debug)
	printf("Root? %d %f %f %f %f\n", node->data, node->perfusion, perfL, perfR, perfL+perfR);		
	}
	else if(node->left!=NULL && node->right==NULL){ 
		node->pressure = finalPressure + (node->perfusion)*(node->downResistance); // mmHg, given		
		perfL  = (node->pressure-finalPressure)/( node->left->resistanceSegment    + node->left->downResistance );
		perfR = 0.0; 
		node->left->perfusion = perfL;
//		node->right->perfusion = perfR;
		pressurePerfusion(node->left);
//		pressurePerfusion(node->right);		

		pressure = fopen("pressure.txt","a+");
		fprintf(pressure, "%d %f\n", node->data, node->pressure);
		fclose(pressure);				
		perfusion = fopen("perfusion.txt","a+");		
		fprintf(perfusion, "%d %f\n", node->data, node->perfusion);
		fclose(perfusion);


	if(debug)
	printf("Root? %d %f %f %f %f\n", node->data, node->perfusion, perfL, perfR, perfL+perfR);
	}
	else if(node->left==NULL && node->right!=NULL){ 
		node->pressure = finalPressure + (node->perfusion)*(node->downResistance); // mmHg, given		
		perfL  = 0.0; 
		perfR = (node->pressure-finalPressure)/( node->right->resistanceSegment + node->right->downResistance );		
//		node->left->perfusion = perfL;
		node->right->perfusion = perfR;				
//		pressurePerfusion(node->left);
		pressurePerfusion(node->right);
		
		pressure = fopen("pressure.txt","a+");
		fprintf(pressure, "%d %f\n", node->data, node->pressure);
		fclose(pressure);				
		perfusion = fopen("perfusion.txt","a+");		
		fprintf(perfusion, "%d %f\n", node->data, node->perfusion);
		fclose(perfusion);
		
	if(debug)		
	printf("Root? %d %f %f %f %f\n", node->data, node->perfusion, perfL, perfR, perfL+perfR);		
	}
	else if(node->left==NULL && node->right==NULL){ 
		node->pressure = 2.0 + (node->perfusion)*(node->downResistance); // mmHg, given		
		perfL  = 0.0; 
		perfR = 0.0; 

		pressure = fopen("pressure.txt","a+");
		fprintf(pressure, "%d %f\n", node->data, node->pressure);
		fclose(pressure);				
		perfusion = fopen("perfusion.txt","a+");		
		fprintf(perfusion, "%d %f\n", node->data, node->perfusion);
		fclose(perfusion);

	if(debug)
	printf("Leaves. %d %f %f %f %f\n", node->data, node->perfusion, perfL, perfR, perfL+perfR);
	}

} // end of void pressurePerfusion function.


/*
DO NOT USE THIS FUNCTION ANY MORE. the job is done in pressurePerfusion function.
---------------------------------------------------------------------------------------------------------------------
Function for perfusion calculation.
1 May 2017:
Perfusion out of a node is not exactly a spatial constant. 
Smaller leaf vessels will have less perfusion through them.
Larger leaf vessals will have larger perfusion through them.
Revise this function to make sure that:
a) The perfusion through a leaf is q0/N if it is the only child.
b) if there are 2 children, then the perfusion is set up proportional to the radius/diameter of the segment.
*/
double binaryPerfusion(int NN, double leafPerf, node *node){ 
	FILE *perfusion;
	double binPerf;
	if(node==NULL) return -1.0; // error trapping, this must never ever happen.
	if(node->left==NULL && node->right==NULL){
			binPerf = leafPerf;
			node->perfusion = binPerf;
			perfusion = fopen("perfusion.txt", "a+");
			fprintf(perfusion, "%d %f\n", node->data, binPerf);
			fclose(perfusion);
	return binPerf;
	}
	else if(node->left==NULL && node->right!=NULL){ 	
			binPerf = binaryPerfusion(node->data, leafPerf, node->right);
			node->perfusion = binPerf;
			perfusion = fopen("perfusion.txt", "a+");
			fprintf(perfusion, "%d %f\n", node->data, binPerf );
			fclose(perfusion);
		return binPerf;
	}
	else if(node->left!=NULL && node->right==NULL){ 	
			binPerf = binaryPerfusion(node->data, leafPerf, node->left);
			node->perfusion = binPerf;
			perfusion = fopen("perfusion.txt", "a+");
			fprintf(perfusion, "%d %f\n", node->data,  binPerf );
			fclose(perfusion);

	return binPerf;
	}
	else if(node->left!=NULL && node->right!=NULL){ 	
			binPerf = binaryPerfusion(node->data, leafPerf, node->left) + binaryPerfusion(node->data, leafPerf, node->right);
			node->perfusion = binPerf;
			perfusion = fopen("perfusion.txt", "a+");
			fprintf(perfusion, "%d %f\n", node->data, binPerf );
			fclose(perfusion);
	return binPerf ;
	}

return -10000; // this return value has exit meaning.	
} // end of binaryPerfusion.


/*_________________________________________________________________________*/
/*_________________________________________________________________________*/
// these functions are used in the rearrangement algorithm of the annealing process.
/* Given a binary tree, print its nodes according to the
  "bottom-up" postorder traversal. */
void printPostorder(node *node){
     if (node == NULL)         return;
	if(node->left!=NULL)     printPostorder(node->left);
     if(node->right!=NULL)     printPostorder(node->right);
     // now deal with the node
     FILE *out;  out = fopen("postorder.txt","a+"); fprintf(out,"%d ", node->data); fclose(out);
}

// two functions used in the 3D code.
/* Given a binary tree, print its nodes in inorder*/
void printInorder(node *node, char *str){
     if (node == NULL)       	    return;
     if(node->left!=NULL)  	     printInorder(node->left, str);
     FILE *out; out = fopen(str,"a+"); fprintf(out,"%d ", node->data); fclose(out);

// write all the other data for this node struct.
/*
int 			parent_data
int 			strahler;
int			parent_srahler;
double 		radius, length, perfusion, pressure, resistance, resistanceSegment, downResistance;
double		x_cood, y_cood, z_cood; // this may also be in an array.
double 		theta1, theta2, theta;
int 			assignedOrNot, onOff;
int 			whichCase;
*/
     char filename[1000] = {0};
     if(node->whichCase==1)
     sprintf(filename, "NodeData.data.%d", 1);
     else
     sprintf(filename, "NodeData.data.%d", 2);	     
     out = fopen(filename, "a+");
     fprintf(out, "%d %d %d %d %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %50.50f %d %d %d\n", node->data, node->parent_data, node->strahler, node->parent_srahler, node->radius, node->length, node->perfusion, node->pressure, node->resistance, node->resistanceSegment, node->downResistance, node->x_cood, node->y_cood, node->z_cood, node->theta1, node->theta2, node->theta, node->assignedOrNot, node->onOff, node->whichCase);
     fclose(out);

     /* now recur on right child */
if(node->right!=NULL)      printInorder(node->right, str);
}
 
/* Given a binary tree, print its nodes in preorder*/
void printPreorder(node *node, char *str){
     if (node == NULL)           return;
     /* first print data of node */
//     printf("%d ", node->data);  
     FILE *out; 
     out = fopen(str,"a+"); // add because you are in a recursion.
     fprintf(out,"%d ", node->data);
     fclose(out);

     /* then recur on left sutree */
     printPreorder(node->left, str);  
     /* now recur on right subtree */
     printPreorder(node->right, str);
}    

/*_________________________________________________________________________*/
/*_________________________________________________________________________*/


/* UTILITY FUNCTIONS USED IN RE-ARRAGEMENT */
/* Function to find index of value in arr[start...end]
   The function assumes that value is present in in[] */
int searchR(int arr[], int strt, int end, int value){
  int i;
  int answer;
  answer = -1;
  for(i = strt; i <= end; i++)  if(arr[i] == value){ answer = i; };
  return answer;
}

/* Helper function that allocates a new node with the
   given data and NULL left and right pointers. */
  node *newNodeR(int data){
  node *temp 	= NULL;
  temp 			= (node *)malloc(sizeof(node));
  temp->data 	= data;
  temp->left 		= NULL;
  temp->right 	= NULL;
  return(temp);
}
 
/* This funtcion is here just to test buildTree() */
void printInorderR(node *node) {
  if (node == NULL)      return;
  /* first recur on left child */
  printInorderR(node->left);
  /* then print the data of node */
  printf("%d ", node->data);
  /* now recur on right child */
  printInorderR(node->right);
}
 
/* Recursive function to construct binary of size len from
   Inorder traversal in[] and Preorder traversal pre[].  Initial values
   of inStrt and inEnd should be 0 and len -1.  The function doesn't
   do any error checking for cases where inorder and preorder
   do not form a tree */
/*
NOTE:
The static attribute tells the compiler that you only want to initialize the variable once, and afterwards it should retain it's data. This will fix your problem because after each recursive call, preIndex wont be reset to zero. Instead, it will contain the value of the last call (in this case, the calling function). This is a problem if
you are iterating over and over as done in simulated annealing.

At the end of the function, you should reset preIndex to 0 so that in the next call, it doesn't still contain the value of the call previous. This is also why a temporary variable is made - so that it can return its value before it's set to 0.
*/

/*
Special variable used by buildTreeR. This HAS to be global static for now, till I work out other more robust ways of doing this.
*/

 static int preIndex;

node *buildTreeR(int in[], int pre[], int inStrt, int inEnd){
// static int preIndex = 0; // I made this a global variable.

  if(inStrt > inEnd){ /* printf("this did happen and this printf gave out something\n"); */ return NULL; } 
  
  /* Pick current node from Preorder traversal using preIndex and increment preIndex */
//  printf("preindex, inStrt, and inEnd: %d %d %d\n", preIndex, inStrt, inEnd);
  node *tNode = newNodeR(pre[preIndex++]); // no use of preIdx.
//  printf("After preindex and inIndex %d %d\n", preIndex, inIndex);
  /* If this node has no children then return */
  if(inStrt == inEnd){ return tNode; } // where are inStrt and inEnd changing?
  /* Else find the index of this node in Inorder traversal */
// 185087107 185087107
// printf("%d \n", tNode->data);
//  exit(0);

  int inIndex = searchR(in, inStrt, inEnd, tNode->data); // inIndex is calculated accurately. So that is okay.
  /* Using index in Inorder traversal, construct left and right subtrees */
  
  tNode->left 		= buildTreeR(in, pre, inStrt, inIndex-1);
  tNode->right 	= buildTreeR(in, pre, inIndex+1, inEnd);

  return tNode;
}

 /*___________________________________________________________________________________
____________________________________________________________________________________*/

void NEWTON(double theta, double phi, double a, double b, double c, double p_x, double p_y, double p_z, double c_x, double c_y, double c_z, double *thetaphi){
double cost, sint, cosp, sinp, a11, a12, a21, a22, f1, f2;
double b11, b12, b21, b22, det, dtheta, dphi, tx, ty, tz, tdist;
int iters;

iters = 0;
do{
	cost 	= cos(theta);  sint 	= sin(theta);
	cosp 	= cos(phi);     	sinp 	= sin(phi);
/*
	f1 = - sint * cost *( a * a * cosp * cosp + b * b * sinp * sinp - c * c) + a * p_x * sint * cosp + b * p_y * sint * sinp - c * p_z * cost;
	f2 =  cost * cost * sinp * sinp *(b * b - a * a) + a * p_x * cost * sinp - b * p_y * cost * cosp;

	a11 = (sint * sint - cost * cost) * ( a * a * cosp * cosp + b * b * sinp * sinp - c * c) + a * p_x * cost * cosp + b * p_y * cost * sinp + c * p_z * sint;
	a12 = 2.0 * sint * cost * cosp * sinp * ( a * a  - b * b ) - a * p_x * sint * sinp + b * p_y * sint * cosp;
	a21 = ( a * a  - b * b ) * sint * sinp * cosp;
	a22 = ( a * a  - b * b ) * cost * (sinp * sinp - cosp * cosp) + a * p_x * cosp + b * p_y * sinp;
*/

	f1 = (a*a - b*b)* cost * sint * cosp - p_x * a * sint + p_y * b * cost;
	f2 = sinp * cosp * ( a * a * cost * cost + b * b * sint * sint - c * c) -p_x * a * sinp * cost - p_y * b * sinp * sint + p_z * c * cosp;
	
	a11 = (a * a - b * b) * (cost * cost - sint * sint) * cosp - p_x * a * cost - p_y * b * sint;
	a12 = (b * b - a * a) * cost * sint * sinp;
	a21 = 2.0 * (b * b - a * a) * cost * sint * sinp * cosp + p_x * a * sinp * sint - p_y * b * sinp * cost;
	a22 = (a * a * cost * cost + b * b * sint * sint - c * c) * (cosp * cosp - sinp * sinp) - p_x * a * cosp * cost - p_y * b * cosp * sint - p_z * c * sinp;
	det = a11 * a22 - a21 * a12;

	if(fabs(det) > 0.0){ b11 = a22 / det;     b12 = -a12 / det; b21 = -a21/det;      b22 = a11/det; }
	
	dtheta 	= b11  * f1 + b12 * f2; dphi 	= b21 * f1 + b22 * f2;

	// revise your approximation.
	theta = theta + dtheta; phi	   = phi + dphi; iters++;
}while(fabs(dtheta) > 0.01 && fabs(dphi) > 0.01);

	thetaphi[0] = theta;
	thetaphi[1] = phi;

	if(debug)
	printf("number of iterations: %d\n", iters);
}

 /*___________________________________________________________________________________
____________________________________________________________________________________*/

double closestAVBDR(XYZ xu_vector, XYZ *avbdr_closest, XYZ *avbdr_normal, double fractionforOrder){
// analytical solution. It was found to be more accurate than the search in grid solution at DX = 0.1. Always take numIters > 1000
// analytical method for the same distance.

	double LV_x_c, LV_y_c, LV_z_c;	
	double tempdistance, distToSurface_avbdr;
	double theta, phi, cost, sint, cosp, sinp, a, b, c, a11, a12, a21, a22, f1, f2, p_x, p_y, p_z, nx, ny, nz;
	double b11, b12, b21, b22, det, dtheta, dphi, tx, ty, tz, tdist, distt, my_tx, my_ty, my_tz;
	int iters;
	double numIters;
	
	numIters = 100.0; 
	LV_x_c = MARGIN + LV_X_RADIUS;  LV_y_c = MARGIN + LV_Y_RADIUS; LV_z_c = MARGIN;
	tdist = 10000000000000000000000.0; distToSurface_avbdr = 1000000000000000.0; // a large impossible value.

// use parametric equation of ellipse: a * cos(theta), b * sin(theta) with theta in [0 , 2PI) and phi in 0, M_PI
p_x = xu_vector.x;   p_y = xu_vector.y;  p_z = xu_vector.z; // you have to shift origin to Lx, Ly, Lz.

// if projection is inside AV Border, then the distance is just the z distance between p_z and AV Border.
if( (((LV_x_c - p_x) * (LV_x_c - p_x)/(LV_X_RADIUS * LV_X_RADIUS) + (LV_y_c - p_y) * (LV_y_c - p_y)/(LV_Y_RADIUS * LV_Y_RADIUS)  <= 1.0)
&& ((LV_x_c-p_x)*(LV_x_c-p_x)/( (LV_X_RADIUS-fractionforOrder*LV_THICKNESS)*(LV_X_RADIUS-fractionforOrder*LV_THICKNESS))+(LV_y_c-p_y)*(LV_y_c-p_y)/( (LV_Y_RADIUS-fractionforOrder*LV_THICKNESS) * (LV_Y_RADIUS-fractionforOrder*LV_THICKNESS))>= 1.0))
|| 
    (((LV_x_c - p_x) * (LV_x_c - p_x)/(RV_X_RADIUS * RV_X_RADIUS)+  (LV_y_c - p_y) * (LV_y_c - p_y)/(RV_Y_RADIUS * RV_Y_RADIUS) <= 1.0)
&&((LV_x_c-p_x)*(LV_x_c-p_x)/( (RV_X_RADIUS-fractionforOrder*RV_THICKNESS)*(RV_X_RADIUS-fractionforOrder*RV_THICKNESS) )+(LV_y_c-p_y)*(LV_y_c-p_y)/( (RV_Y_RADIUS-fractionforOrder*RV_THICKNESS) * (RV_Y_RADIUS-fractionforOrder*RV_THICKNESS))>=1.0)
&&( (LV_x_c - p_x) * (LV_x_c - p_x)/(LV_X_RADIUS * LV_X_RADIUS)+ (LV_y_c - p_y) * (LV_y_c - p_y)/(LV_Y_RADIUS * LV_Y_RADIUS) > 1.0)
&&(p_y > RV_X_RADIUS+(double)9.0) 
)
){
	tdist 	= (p_z - (double)MARGIN)*(p_z - (double)MARGIN); // final distance.
	my_tx 	= p_x; my_ty = p_y; my_tz = (double)MARGIN ;
}
else if( // inside LV cavity is easy, no need to do anything else.
((LV_x_c - p_x) * (LV_x_c - p_x)/( (LV_X_RADIUS -fractionforOrder*LV_THICKNESS) * (LV_X_RADIUS-fractionforOrder*LV_THICKNESS) ) 
    + (LV_y_c - p_y) * (LV_y_c - p_y)/( (LV_Y_RADIUS-fractionforOrder*LV_THICKNESS) * (LV_Y_RADIUS-fractionforOrder*LV_THICKNESS) ) 
    < 1.0)
) {
		// inner LV edge.
a = LV_X_RADIUS -fractionforOrder*LV_THICKNESS; b = LV_Y_RADIUS-fractionforOrder*LV_THICKNESS; c=LV_Z_RADIUS-fractionforOrder*LV_THICKNESS;
		for(theta = 0; theta < 2.0 * M_PI; theta = theta + 2.0 * M_PI / numIters){
			tx = a * cos(theta) + LV_x_c;  ty = b * sin(theta) + LV_y_c;
			distt = (tx - p_x) * (tx - p_x) + (ty - p_y) * (ty - p_y);
			if(distt < tdist ){ tdist = distt; my_tx = tx; my_ty = ty; my_tz = (double)MARGIN ; }
		}
} else {
		// in RV cavity, LV outside.
		a = LV_X_RADIUS;     b = LV_Y_RADIUS;    c = LV_Z_RADIUS;
		for(theta = 0; theta < 2.0 * M_PI; theta = theta + 2.0 * M_PI / numIters){
			tx = a * cos(theta) + LV_x_c;  ty = b * sin(theta) + LV_y_c;
			distt = (tx - p_x) * (tx - p_x) + (ty - p_y) * (ty - p_y);
			if(distt < tdist && ty > 53.3 ){ tdist = distt; my_tx = tx; my_ty = ty; my_tz = (double)MARGIN ; } // the ty value was found by looking at geometry in paravu.
		}
		//  in RV cavity, RV inside.
a=RV_X_RADIUS-fractionforOrder*RV_THICKNESS; b=RV_Y_RADIUS-fractionforOrder*RV_THICKNESS; c=RV_Z_RADIUS-fractionforOrder*RV_THICKNESS;
		for(theta = 0; theta < 2.0 * M_PI; theta = theta + 2.0 * M_PI / numIters){
			tx = a * cos(theta) + LV_x_c; ty = b * sin(theta) + LV_y_c;
			distt = (tx - p_x) * (tx - p_x) + (ty - p_y) * (ty - p_y);
			if(distt < tdist && ty > 53.3 ){ tdist = distt; my_tx = tx; my_ty = ty; my_tz = (double)MARGIN ; } // the ty value was found by looking at geometry in paravu.
		}
				tdist = tdist + (p_z - (double)MARGIN)*(p_z - (double)MARGIN);
} // end of big if.

// outside surface condition does not happen because AV border is the broadest part of all the ellipsoids in Z direction.

		avbdr_normal->x = 0.0; 	avbdr_normal->y = 0.0; 	avbdr_normal->z = -1.0;
		avbdr_closest->x = my_tx; 	avbdr_closest->y = my_ty; 	avbdr_closest->z = my_tz; 
		distToSurface_avbdr = tdist;
if(debug){
		printf("In AV border, distance: %f\n", sqrt(distToSurface_avbdr) );
		printf("Coordinate at border: %f %f %f\n", avbdr_closest->x, avbdr_closest->y, avbdr_closest->z );
}

return sqrt(distToSurface_avbdr);
}

/*
double closestLVCHAMBER(XYZ xu_vector, XYZ *LVCHAMBER_closest, XYZ *LVCHAMBER_normal, double fractionforOrder){
	double LV_x_c, LV_y_c, LV_z_c;	
// analytical method.
	double theta, phi, cost, sint, cosp, sinp, a, b, c, a11, a12, a21, a22, f1, f2, p_x, p_y, p_z, nx, ny, nz;
	double b11, b12, b21, b22, det, dtheta, dphi, tx, ty, tz, tdist, distToSurface_LVCHAMBER;
	double thetaphi[2];
	
	LV_x_c = MARGIN + LV_X_RADIUS;
	LV_y_c = MARGIN + LV_Y_RADIUS;
	LV_z_c = MARGIN;
	distToSurface_LVCHAMBER = 1000000000000000.0; // a large value.

	a = LV_X_RADIUS-fractionforOrder*LV_THICKNESS;  b = LV_Y_RADIUS-fractionforOrder*LV_THICKNESS; c = LV_Z_RADIUS-fractionforOrder*LV_THICKNESS;
	p_x 	= xu_vector.x - LV_x_c;  	p_y = xu_vector.y - LV_y_c;  					p_z  = xu_vector.z - LV_z_c; // you have to shift origin to Lx, Ly, Lz.
// initial theta and phi. 
	theta 	= atan2(a * p_y, b * p_x ) ;
	phi 		= atan2(p_z, c * sqrt( p_x * p_x / (a * a) + p_y* p_y / (b * b) ) ); 
	
//	theta = atan2( p_z * sqrt(a*a+b*b), c * sqrt(p_x * p_x + p_y * p_y)  );
//	phi = atan2( a * p_y, b * p_x);

	// thetaphi is the return value.
	NEWTON(theta, phi, a, b, c, p_x, p_y, p_z, LV_x_c, LV_y_c, LV_z_c, thetaphi);
	theta 	= thetaphi[0]						;
	phi 		= thetaphi[1]						;
	cost 	= cos(theta);  sint 	= sin(theta)		;
	cosp 	= cos(phi);     	sinp 	= sin(phi)	;

// calculate your x location on the ellipse. This needs to be translated to my models coordinates before returning to main.
	tx 						= a * cosp * cost + LV_x_c; 		ty			     = b * cosp * sint + LV_y_c; 	tz 				= c * sinp + LV_z_c;
	p_x = p_x + LV_x_c; p_y = p_y + LV_y_c; p_z = p_z + LV_z_c;	
	
	nx 						= -2.0 * ( tx -  LV_x_c )/( a * a); 	ny 					     = -2.0 * ( ty - LV_y_c)/(  b *b); 	nz 				= -2.0*(tz-LV_z_c)/(c *c);
	
	
	LVCHAMBER_closest->x 	= tx; 						LVCHAMBER_closest->y  = ty; 						LVCHAMBER_closest->z = tz;
	LVCHAMBER_normal->x 	= nx; 						LVCHAMBER_normal->y = ny; 						LVCHAMBER_normal->z = nz;

	distToSurface_LVCHAMBER = (p_x - tx) * (p_x - tx) + (p_y - ty) * (p_y - ty) + (p_z - tz) * (p_z - tz);
if(debug){
	printf("In LV chamber, distance: %f\n", sqrt(distToSurface_LVCHAMBER) );
	printf("Coordinate at LV distance: %f %f %f\n", LVCHAMBER_closest->x, LVCHAMBER_closest->y, LVCHAMBER_closest->z = tz );
}

return sqrt(distToSurface_LVCHAMBER);
}

*/

double closestLVCHAMBER(XYZ xu_vector, XYZ *LVCHAMBER_closest, XYZ *LVCHAMBER_normal, double fractionforOrder){
	double LV_x_c, LV_y_c, LV_z_c;	
// analytical method didnt do it for me. I am now going simple.
	double theta, phi, cost, sint, cosp, sinp, a, b, c, a11, a12, a21, a22, f1, f2, p_x, p_y, p_z, nx, ny, nz;
	double b11, b12, b21, b22, det, dtheta, dphi, tx, ty, tz, tdist, distToSurface_LVCHAMBER;
	double numIters, tempdistance;
	
	LV_x_c = MARGIN + LV_X_RADIUS;
	LV_y_c = MARGIN + LV_Y_RADIUS;
	LV_z_c = MARGIN;
	distToSurface_LVCHAMBER = 1000000000000000.0; // a large value.
	numIters = 1000.0;
	a = LV_X_RADIUS-fractionforOrder*LV_THICKNESS;  b = LV_Y_RADIUS-fractionforOrder*LV_THICKNESS; c = LV_Z_RADIUS-fractionforOrder*LV_THICKNESS;
	p_x  = xu_vector.x;								   p_y  = xu_vector.y;  							       p_z   = xu_vector.z; 	
	
/* Rather than do an i,j,k grid iteration, do a theta phi iteration. That may be a few 1000 times faster. */
// condition for inside RV wall.
for(theta = 0; theta < 2 * M_PI; theta = theta + 2 * M_PI /numIters )
for(phi = 0;     phi < M_PI;             phi    = phi + M_PI /numIters        ) {
	tx = a * cos(phi) * cos(theta) + LV_x_c; ty = b * cos(phi) * sin(theta) + LV_y_c; tz  = c * sin(phi) + LV_z_c;	

	tempdistance = ( p_x - tx ) * (  p_x - tx )  +	( p_y - ty ) * (  p_y - ty ) + ( p_z -  tz ) * (  p_z -  tz );
	 if(tempdistance < distToSurface_LVCHAMBER ){
	distToSurface_LVCHAMBER = tempdistance;
	// note the location.
	LVCHAMBER_closest->x  = tx ; LVCHAMBER_closest->y = ty ; LVCHAMBER_closest->z = tz ;
//	nx 						= -2.0 * ( tx -  LV_x_c )/( a * a); 	ny 					     = -2.0 * ( ty - LV_y_c)/(  b *b); 	nz 				= -2.0*(tz-LV_z_c)/(c *c);	
	LVCHAMBER_normal->x = -2.0 * ( tx - LV_x_c )/( a*a ); LVCHAMBER_normal->y = -2.0 * ( ty - LV_y_c )/( b*b ); LVCHAMBER_normal->z = -2.0 * ( tz - LV_z_c )/( c*c );
         }
}
	return sqrt(distToSurface_LVCHAMBER);
}



double closestRVCHAMBER(XYZ xu_vector, XYZ *RVCHAMBER_closest, XYZ *RVCHAMBER_normal, double fractionforOrder){
	int usr_i, usr_j, usr_k;
	double LV_x_c, LV_y_c, LV_z_c;	
	double tempdistance, distToSurface_RVCHAMBER;
	double  my_tx, my_ty, my_tz, tx, ty, tz, a, b, c, p_x, p_y, p_z, aa, bb, cc;
	double theta, phi; // theta is in 0, 2PI and phi is in 0, PI.
	double numIters;
	double stheta, sphi;
	
	LV_x_c = MARGIN + LV_X_RADIUS;
	LV_y_c = MARGIN + LV_Y_RADIUS;
	LV_z_c = MARGIN;
	distToSurface_RVCHAMBER = 1000000000000000.0; // a large value.

// use parametric equation of ellipse: a * cos(theta), b * sin(theta) with theta in [0 , 2PI)     
// test if answer is inside, it should not be.
numIters = 1000.0;
a=RV_X_RADIUS-fractionforOrder*RV_THICKNESS; b = RV_Y_RADIUS-fractionforOrder*RV_THICKNESS; c = RV_Z_RADIUS - fractionforOrder*RV_THICKNESS;
p_x  = xu_vector.x;								   p_y  = xu_vector.y;  							       p_z   = xu_vector.z; 

/* Rather than do an i,j,k grid iteration, do a theta phi iteration. That may be a few 1000 times faster. */
// condition for inside RV wall.
for(theta = 0; theta < 2 * M_PI; theta = theta + 2 * M_PI /numIters )
for(phi = 0;     phi < M_PI;             phi    = phi + M_PI /numIters        ) {
	tx = a * cos(phi) * cos(theta) + LV_x_c; ty = b * cos(phi) * sin(theta) + LV_y_c; tz  = c * sin(phi) + LV_z_c;	
//	if(tz>2.0) 	printf("%f %f %f\n", tx, ty, tz);
	tempdistance = ( p_x - tx ) * (  p_x - tx )  +	( p_y - ty ) * (  p_y - ty ) + ( p_z -  tz ) * (  p_z -  tz );
	 if(tempdistance < distToSurface_RVCHAMBER && (ty > 53.3) ){
	distToSurface_RVCHAMBER = tempdistance;
	// note the location.
	RVCHAMBER_closest->x  = tx ; RVCHAMBER_closest->y = ty ; RVCHAMBER_closest->z = tz ;
	RVCHAMBER_normal->x = -2.0 * ( tx - LV_x_c )/( a*a ); RVCHAMBER_normal->y = -2.0 * ( ty - LV_y_c )/( b*b ); RVCHAMBER_normal->z = -2.0 * ( tz - LV_z_c )/( c*c );
	stheta = theta; sphi = phi;
         }
}

// see if the RV sides LV outer surface is any closer.
a = LV_X_RADIUS;     b = LV_Y_RADIUS;    c = LV_Z_RADIUS;
for(theta = 0; theta < 2.0 * M_PI; theta = theta + 2.0 * M_PI /numIters )
for(phi = 0;     phi < M_PI;             phi    = phi + M_PI /numIters              ){
	tx = a * cos(phi) * cos(theta) + LV_x_c; ty = b * cos(phi) * sin(theta) + LV_y_c; tz  = c * sin(phi) + LV_z_c;
	tempdistance = ( p_x - tx ) * (  p_x - tx )  +	( p_y - ty ) * (  p_y - ty ) + ( p_z -  tz ) * (  p_z -  tz );	
	                if(tempdistance < distToSurface_RVCHAMBER && ( ty > 53.3) ){	
	distToSurface_RVCHAMBER = tempdistance;
	// note the location.
	RVCHAMBER_closest->x = tx ; RVCHAMBER_closest->y = ty ; RVCHAMBER_closest->z = tz ;
	RVCHAMBER_normal->x = -2.0 * ( tx - LV_x_c )/( a*a ); RVCHAMBER_normal->y = -2.0 * ( ty - LV_y_c )/( b*b ); RVCHAMBER_normal->z = -2.0 * ( tz - LV_z_c )/( c*c );
		stheta = theta; sphi = phi;
	                }
}
if(debug){
		printf("In RV chamber, distance: %f\n", sqrt(distToSurface_RVCHAMBER) );
		printf("Coordinate at RV distance: %f %f %f\n", RVCHAMBER_closest->x, RVCHAMBER_closest->y, RVCHAMBER_closest->z );
		printf("Angles: %f %f\n", stheta, sphi);
}
return sqrt(distToSurface_RVCHAMBER);
}

double closestOUTSIDE(XYZ xu_vector, XYZ *OUTSIDE_closest, XYZ *OUTSIDE_normal, double fractionforOrder){
	
	double tempdistance, distToSurface_OUTSIDE;
	double  my_tx, my_ty, my_tz, tx, ty, tz, a, b, c, p_x, p_y, p_z, aa, bb, cc;
	double theta, phi; // theta is in 0, 2PI and phi is in 0, PI.
	double numIters;
	double LV_x_c, LV_y_c, LV_z_c;		
	
	numIters = 1000.0;		
	LV_x_c = MARGIN + LV_X_RADIUS; 	LV_y_c = MARGIN + LV_Y_RADIUS; 	LV_z_c = MARGIN;
	distToSurface_OUTSIDE = 1000000000000000.0; // a large value.
	a = LV_X_RADIUS;     b = LV_Y_RADIUS;    c = LV_Z_RADIUS;  aa = RV_X_RADIUS;     bb = RV_Y_RADIUS;    cc = RV_Z_RADIUS;
	p_x = xu_vector.x; p_y = xu_vector.y; p_z = xu_vector.z;

// LV part.
for(theta = 0; theta < 2.0 * M_PI; theta = theta + 2.0 * M_PI /numIters )
for(phi = 0;     phi < M_PI;             phi    = phi + M_PI /numIters )               {
	tx = a * cos(phi) * cos(theta) + LV_x_c; ty = b * cos(phi) * sin(theta) + LV_y_c; tz  = c * sin(phi) + LV_z_c;
	if( ( ( (LV_x_c - tx) * (LV_x_c - tx)/( a * a ) + (LV_y_c - ty) * (LV_y_c - ty)/( b * b ) + (LV_z_c - tz) * (LV_z_c - tz)/( c * c ) - 1.0) <= 0.01 )  // square of DX.
	       && ( ty < RV_X_RADIUS+(double)9.0 )
	 ){
			tempdistance = ( p_x - tx ) * ( p_x - tx ) + ( p_y - ty ) * ( p_y - ty ) + ( p_z - tz ) * ( p_z - tz ) ;
			// normal on LV surface
			if(tempdistance < distToSurface_OUTSIDE){
				distToSurface_OUTSIDE = tempdistance;
				OUTSIDE_closest->x = tx; OUTSIDE_closest->y = ty; OUTSIDE_closest->z = tz;
				OUTSIDE_normal->x = 2.0 * ( tx - LV_x_c) /( a*a ); OUTSIDE_normal->y = 2.0 * ( ty - LV_y_c) /( b*b ); OUTSIDE_normal->z = 2.0 * ( tz - LV_z_c) /( c*c );
			}
	 }
}

// RV part.
for(theta = 0; theta < 2.0 * M_PI; theta = theta + 2.0 * M_PI /numIters )
for(phi = 0;     phi < M_PI;             phi    = phi + M_PI /numIters ){
// the abc are for RV.	
tx = aa * cos(phi) * cos(theta) + LV_x_c; ty = bb * cos(phi) * sin(theta) + LV_y_c; tz  = cc * sin(phi) + LV_z_c;	
	if( ( (LV_x_c - tx) * (LV_x_c - tx)/( a * a ) + (LV_y_c - ty) * (LV_y_c - ty)/( b * b ) + (LV_z_c - tz) * (LV_z_c - tz)/( c * c ) > 1.0) &&( ty > RV_X_RADIUS+(double)9.0 ) ) {
		tempdistance = ( p_x - tx ) * ( p_x - tx )  + ( p_y - ty ) * ( p_y - ty )  + ( p_z -  tz ) * ( p_z - tz );			
		// normal on RV surface.
		if(tempdistance < distToSurface_OUTSIDE){
			distToSurface_OUTSIDE = tempdistance;
			OUTSIDE_closest->x = tx ; OUTSIDE_closest->y = ty ; OUTSIDE_closest->z = tz ;
			OUTSIDE_normal->x=2.0*( tx-LV_x_c) /(aa*aa ); OUTSIDE_normal->y=2.0*(ty-LV_y_c) /(bb*bb); OUTSIDE_normal->z=2.0*(tz-LV_z_c) /(cc*cc);
		}
	}
}

if(debug){
		printf("Outside, distance: %f\n", sqrt(distToSurface_OUTSIDE) );
		printf("Coordinate outside surface: %f %f %f\n", OUTSIDE_closest->x, OUTSIDE_closest->y, OUTSIDE_closest->z );
}

return sqrt(distToSurface_OUTSIDE);
}

int isNodeInTissue(double x, double y, double z, double fractionforOrder){

	int nodeInside;
	double LV_x_c, LV_y_c, LV_z_c;
	
	// this is the centre of all 4 ellipsoids.
	LV_x_c = MARGIN + LV_X_RADIUS;
	LV_y_c = MARGIN + LV_Y_RADIUS;
	LV_z_c = MARGIN;

			if ( ( ( ( (LV_x_c - x) * (LV_x_c - x)/((LV_X_RADIUS+0.02) * (LV_X_RADIUS+0.02)) 
			       + (LV_y_c - y) * (LV_y_c - y)/((LV_Y_RADIUS+0.02) * (LV_Y_RADIUS+0.02)) 
			       + (LV_z_c - z) * (LV_z_c - z)/((LV_Z_RADIUS+0.02) * (LV_Z_RADIUS+0.02)) <= 1.05)
			&&
				( (LV_x_c - x) * (LV_x_c - x)/( (LV_X_RADIUS -fractionforOrder*LV_THICKNESS) * (LV_X_RADIUS-fractionforOrder*LV_THICKNESS) ) 
			       + (LV_y_c - y) * (LV_y_c - y)/( (LV_Y_RADIUS-fractionforOrder*LV_THICKNESS) * (LV_Y_RADIUS-fractionforOrder*LV_THICKNESS) ) 
			       + (LV_z_c - z) * (LV_z_c - z)/( (LV_Z_RADIUS-fractionforOrder*LV_THICKNESS) * (LV_Z_RADIUS-fractionforOrder*LV_THICKNESS) ) > 0.98)
				) 
||
		(		( (LV_x_c - x) * (LV_x_c - x)/((LV_X_RADIUS-0.02) * (LV_X_RADIUS-0.02)) 
			       + (LV_y_c - y) * (LV_y_c - y)/((LV_Y_RADIUS-0.02) * (LV_Y_RADIUS-0.02)) 
			       + (LV_z_c - z) * (LV_z_c - z)/((LV_Z_RADIUS-0.02) * (LV_Z_RADIUS-0.02)) > 0.98)
			&&
				( (LV_x_c - x) * (LV_x_c - x)/((RV_X_RADIUS+0.01) * (RV_X_RADIUS+0.01)) 
			       + (LV_y_c - y) * (LV_y_c - y)/((RV_Y_RADIUS+0.01) * (RV_Y_RADIUS+0.01)) 
			       + (LV_z_c - z) * (LV_z_c - z)/((RV_Z_RADIUS+0.01) * (RV_Z_RADIUS+0.01)) <= 1.05)
			&&
				( (LV_x_c - x) * (LV_x_c - x)/( (RV_X_RADIUS -fractionforOrder*RV_THICKNESS) * (RV_X_RADIUS-fractionforOrder*RV_THICKNESS) ) 
			       + (LV_y_c - y) * (LV_y_c - y)/( (RV_Y_RADIUS-fractionforOrder*RV_THICKNESS) * (RV_Y_RADIUS-fractionforOrder*RV_THICKNESS) ) 
			       + (LV_z_c - z) * (LV_z_c - z)/( (RV_Z_RADIUS-fractionforOrder*RV_THICKNESS) * (RV_Z_RADIUS-fractionforOrder*RV_THICKNESS) ) > 0.98)
			&&  ( y >= RV_X_RADIUS+(double)6.0 )
				) )			
				&& x>=0 && y>=0 && z>= LV_z_c )
					nodeInside = 1;  else  nodeInside = 0;

	return nodeInside;
}

 /*________________________________________________________________________________________________________________________*/

// check the AV border condition: 8 August 2017.
void writeVentricleInCoordinates(){

	char str[10000];
	FILE *output;
	int usr_i, usr_j, usr_k;
	float LV_x_c, LV_y_c, LV_z_c;
	float xtmp, ytmp, ztmp;
	
	// this is the centre of all 4 ellipsoids.
	LV_x_c = MARGIN + LV_X_RADIUS;
	LV_y_c = MARGIN + LV_Y_RADIUS;
	LV_z_c = MARGIN;
// header
printf("%d %d %d\n", usr_MX, usr_MY, usr_MZ);
  		sprintf(str,"vascVentricle.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING %f %f %f\n", DX, DX, DX);
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n", usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS huamnVentricle int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
// do some calculation.		
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
			if(	( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/(LV_X_RADIUS * LV_X_RADIUS) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/(LV_Y_RADIUS * LV_Y_RADIUS) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/(LV_Z_RADIUS * LV_Z_RADIUS) <= 1.0)
			&&
				( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/( (LV_X_RADIUS -LV_THICKNESS) * (LV_X_RADIUS-LV_THICKNESS) ) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/( (LV_Y_RADIUS-LV_THICKNESS) * (LV_Y_RADIUS-LV_THICKNESS) ) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/( (LV_Z_RADIUS-LV_THICKNESS) * (LV_Z_RADIUS-LV_THICKNESS) ) >= 1.0)
				) 
				fprintf(output," 100 ");
		else if(	
				( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/(LV_X_RADIUS * LV_X_RADIUS) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/(LV_Y_RADIUS * LV_Y_RADIUS) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/(LV_Z_RADIUS * LV_Z_RADIUS) > 1.0)
			&&
				( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/(RV_X_RADIUS * RV_X_RADIUS) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/(RV_Y_RADIUS * RV_Y_RADIUS) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/(RV_Z_RADIUS * RV_Z_RADIUS) <= 1.0)
			&&
				( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/( (RV_X_RADIUS -RV_THICKNESS) * (RV_X_RADIUS-RV_THICKNESS) ) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/( (RV_Y_RADIUS-RV_THICKNESS) * (RV_Y_RADIUS-RV_THICKNESS) ) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/( (RV_Z_RADIUS-RV_THICKNESS) * (RV_Z_RADIUS-RV_THICKNESS) ) >= 1.0)
			&&  ( (double)usr_j * DX > RV_X_RADIUS+(double)6.0 )
				) 
				fprintf(output," 200 ");				
				else
				fprintf(output," 0 ");					
			}
				fprintf(output,"\n");			
			}
			
		fclose(output);
		
// epicardial surface or a small volume.
  		sprintf(str,"vascVentricleEpiSurface.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING %f %f %f\n", DX, DX, DX);
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n", usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS huamnVentricle int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
// do some calculation.		
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
			if(	( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/(LV_X_RADIUS * LV_X_RADIUS) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/(LV_Y_RADIUS * LV_Y_RADIUS) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/(LV_Z_RADIUS * LV_Z_RADIUS) <= 1.0)
			&&
				( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/( (LV_X_RADIUS -0.05*LV_THICKNESS) * (LV_X_RADIUS-0.05*LV_THICKNESS) ) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/( (LV_Y_RADIUS-0.05*LV_THICKNESS) * (LV_Y_RADIUS-0.05*LV_THICKNESS) ) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/( (LV_Z_RADIUS-0.05*LV_THICKNESS) * (LV_Z_RADIUS-0.05*LV_THICKNESS) ) >= 1.0)
				) 
				fprintf(output," 100 ");
		else if( 	( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/(LV_X_RADIUS * LV_X_RADIUS) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/(LV_Y_RADIUS * LV_Y_RADIUS) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/(LV_Z_RADIUS * LV_Z_RADIUS) > 1.0)
			&&
				  ( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/(RV_X_RADIUS * RV_X_RADIUS) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/(RV_Y_RADIUS * RV_Y_RADIUS) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/(RV_Z_RADIUS * RV_Z_RADIUS) <= 1.0)
			&&
				( (LV_x_c - (double)usr_i * DX) * (LV_x_c - (double)usr_i * DX)/( (RV_X_RADIUS-0.05*RV_THICKNESS) * (RV_X_RADIUS-0.05*RV_THICKNESS) ) 
			       + (LV_y_c - (double)usr_j * DX) * (LV_y_c - (double)usr_j * DX)/( (RV_Y_RADIUS-0.05*RV_THICKNESS) * (RV_Y_RADIUS-0.05*RV_THICKNESS) ) 
			       + (LV_z_c - (double)usr_k * DX) * (LV_z_c - (double)usr_k * DX)/( (RV_Z_RADIUS-0.05*RV_THICKNESS) * (RV_Z_RADIUS-0.05*RV_THICKNESS) ) >= 1.0)
			&&  ( (double)usr_j * DX > RV_X_RADIUS+(double)6.0 )
				) 
				fprintf(output," 200 ");				
				else
				fprintf(output," 0 ");					
			}
				fprintf(output,"\n");			
			}
			
		fclose(output);
}


// this function could not count the number of nodes in the new Kassab tree.
int count(node *tree){
    int c =  1;             //Node itself should be counted
         if (tree==NULL) return 1;  
    else{
        c += count(tree->left);  c += count(tree->right);
        return c;
    }
}

int binarytree_count_recursive(node *root)
{
    int count = 1;
    if (root->left != NULL) {
       count += binarytree_count_recursive(root->left);
    }
    if (root->right != NULL) {
        count += binarytree_count_recursive(root->right);
    }
    return count;
}
 
int binarytree_count(node *root) {
    int count = 0;
        count = binarytree_count_recursive(root);
    return count;
}


// count number of nodes of order X coming off element of order Y.
int binarytree_count_recursive_atOrder(node *root, int rootOrder, int childOrder)
{
    int count = 0;

    if(root->strahler==rootOrder){
	// the element goes along left, if there is a childOrder its at the end or on the right.    
	    if (root->right != NULL  && root->right->strahler==childOrder ) 	count ++; 
	    if (root->left != NULL && root->left->strahler > childOrder) 	       	count+=binarytree_count_recursive_atOrder(root->left, rootOrder, childOrder);
	    else if(root->left != NULL && root->left->strahler == childOrder) 	count++;
    }
  
	    return count;
}
 
int binarytree_count_atOrder(node *root, int rootOrder, int childOrder) {
    int count = 0;
        count = binarytree_count_recursive_atOrder(root, rootOrder, childOrder);
    return count;
}

double lengthh(double mean, double sigma){
	double X1, X2, U1, U2;

		do{
			// generate 2 random numbers in (0, 1), uniformly distributed.
			U1 = genrand_real1();  							U2 = genrand_real2();
			X1 = 										X2  = -1.0;
			// from normal dbn with mean 0 and std. deviation 1.
			X1 = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2); 		X2 = sqrt(-2.0 * log(U1)) * sin(2.0 * M_PI * U2);
			// from normal dbn with mean = given mean and std. deviation as given.
			X1 = sigma * X1 + mean; 						X2 = sigma * X2 + mean;
			U1 = 										U2 = -1.0;
		}while(X1 < 0.0); // X1 should be more than 0.5 so that you get at least 1 numberOfSegments.
			return X1;
}

int numSegsF(double mean, double sigma){
	double X1, X2, U1, U2;
	int numberOfSegmentsInElement2;
		do{
			// generate 2 random numbers in (0, 1), uniformly distributed.
			U1 = genrand_real1();  							U2 = genrand_real2();
			X1 = 										X2  = -1.0;
			// from normal dbn with mean 0 and std. deviation 1.
			X1 = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2); 		X2 = sqrt(-2.0 * log(U1)) * sin(2.0 * M_PI * U2);
			// from normal dbn with mean = given mean and std. deviation as given.
			X1 = sigma * X1 + mean; 						X2 = sigma * X2 + mean;
			U1 = 										U2 = -1.0;
		}while(X1 < 0.5); // X1 should be more than 0.5 so that you get at least 1 numberOfSegments.			
			// just take X1, ignore X2.
			numberOfSegmentsInElement2 = (int)(X1);
//			if(numberOfSegmentsInElement2 < 1)  numberOfSegmentsInElement2 = 1; // there is no 0 segment case here. Generalise after.

			return numberOfSegmentsInElement2;
}

// tree generator using kassab data. makeKassabTree. This is the engine of the topology generator.
// some large vessels are assigned, the rest are generated.
node *makeKassabTree(int numberOfSegmentsInElement ,unsigned int *nodeData, int order, double *segmentElementRatio_mean, double *segmentElementRatio_SD, double **cuprob_conn , double *segmentLengths_mean,  double *segmentLengths_SD, double *elementDiameters_mean, int parent_data, int whichCase, int stopOrder, int startOrder, double *amx, double x0, double y0, double z0, int onoff, int assigned, double t0, int parent_order, int segLabel ){

	int order2, numberOfSegmentsInElement2;
	double X1, X2, U1, U2, sigma, mean;
	int tem, row;
	double eta_rel, gama, mu, temp_term, a;
	int mainArterySegs, maxOrder;
	int maxOrderNum;
	
	int callSegLevel1, callSegLevel2, callSegLevel3, callSegLevel4, callSegLevel5, callSegLevel6, callSegLevel7, callSegLevel8, callSegLevel9, callSegLevel10, callSegLevel11, callSegLevel12;
callSegLevel1= callSegLevel2= callSegLevel3= callSegLevel4= callSegLevel5= callSegLevel6= callSegLevel7= callSegLevel8= callSegLevel9= callSegLevel10= callSegLevel11= callSegLevel12= 1;

	double n, assym;
	
	double LV_x_c, LV_y_c, LV_z_c;		
	LV_x_c = MARGIN + LV_X_RADIUS; 	LV_y_c = MARGIN + LV_Y_RADIUS; 	LV_z_c = MARGIN;
	n = (segmentElementRatio_mean[order] - (double)numberOfSegmentsInElement);
	assym = pow(amx[order], n);
	if(debug) 	printf("%f\n", assym);

	// create the node and its data.	This applies throughout the function, so do it once.
	node *tNode 		= 	(node *)malloc(sizeof(node));
	do{ tem 		=      (int)(  (double)billion * genrand_real1() + 0.5 ); } while(nodeData[tem]<=0);
	tNode->data 	=      nodeData[tem];
	tNode->parent_data = parent_data;
	nodeData[tem] 	=      -1; // this node is now assigned, untouchable.
	tNode->strahler 	=      order;
	tNode->parent_srahler = parent_order;
	tNode->radius    	= elementDiameters_mean[order]/2.0;  // ***
	mean 			= segmentLengths_mean[order];
	sigma 			= segmentLengths_SD[order];
	tNode->length 	=         lengthh(mean, sigma);
	tNode->resistanceSegment 	= 8.0 * mu_0 * (tNode->length)/ ( M_PI * (tNode->radius) * (tNode->radius) * (tNode->radius) * (tNode->radius) );	
	tNode->whichCase = whichCase;	
 	tNode->assignedOrNot = 0;  tNode->x_cood = -1; tNode->y_cood = -1; tNode->z_cood = -1;  tNode->onOff = 0;		
	tNode->left = NULL; tNode->right = NULL;

	
// there are 3 cases: numberOfSegmentsInElement > 0, or when num segs is 0 and order > stopOrder, or when num segs is 0 and order is stop order.	
// 00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
if(numberOfSegmentsInElement > 0){

// all coordinates for in element,
#include "RCAO11.c" // RCA O11 main artery.
#include "RCAsegLabel1.c"	
#include "RCAsegLabel3.c"
#include "RCAsegLabel5.c"
#include "RCAsegLabel7.c"
#include "LAO11.c"
#include "LCXO10.c"
#include "LAsegLabel2.c"
#include "LAsegLabel4.c"
#include "LAsegLabel6.c"

	
// coordinates must be assigned to this node before the following new nodes are generated.	
/* left node is the same element continuing. */
#include "callLeft.c"

// right branch. This is where you add O10 elements. To do this, you have to know the O11 coordinates - the conditions are hard coded.
// a whole bunch of if else statements: which right statement are you calling.
		if( (fabs(tNode->x_cood-32.0)<0.001) && (fabs(tNode->y_cood-59.0)<0.001) && (fabs(tNode->z_cood-2.0)<0.001)  ){ // RCA septal O(10)
			numberOfSegmentsInElement2 = 9	; order2 = 10	;  whichCase = 1; segLabel = 1;
			#include "callRight.c"
		}
		else if( (fabs(tNode->x_cood-50.5)<0.001) && (fabs(tNode->y_cood-72.0)<0.001) && (fabs(tNode->z_cood-7.0)<0.001)  ){ 
			numberOfSegmentsInElement2 	= 9	; order2 = 10	; whichCase = 1; segLabel = 3	; 
			#include "callRight.c"
		}
		else if( (fabs(tNode->x_cood-13)<0.001) && (fabs(tNode->y_cood-65)<0.001) && (fabs(tNode->z_cood-28)<0.001)  ){ 
			numberOfSegmentsInElement2 	= 9	; order2 = 10	; whichCase = 1	; segLabel = 5	; 
			#include "callRight.c"
		}
		else if( (fabs(tNode->x_cood-4.2)<0.001) && (fabs(tNode->y_cood-29)<0.001) && (fabs(tNode->z_cood-26.0)<0.001)  ){ 
			numberOfSegmentsInElement2 	= 9	; order2 = 10	; whichCase = 1; segLabel = 7	;
			#include "callRight.c"
		}
		else if( (fabs(tNode->x_cood-37.0)<0.001) && (fabs(tNode->y_cood-49.0)<0.001) && (fabs(tNode->z_cood-2.0)<0.001)  ){ // LA septal O(10)
		numberOfSegmentsInElement2 = 9	; order2 	= 10	; whichCase = 2	; segLabel = 2	; // LA septal.
		#include "callRight.c"
		}
		else if( (fabs(tNode->x_cood-61.911180)<0.001) && (fabs(tNode->y_cood-36)<0.001) && (fabs(tNode->z_cood-7.382513)<0.001)  ){ 
		numberOfSegmentsInElement2 	= 9	; order2 = 10	; whichCase = 2; segLabel = 4; 
		#include "callRight.c"
		}
		else if( (fabs(tNode->x_cood-52.937464)<0.001) && (fabs(tNode->y_cood-36)<0.001) && (fabs(tNode->z_cood-52.132610)<0.001)  ){ 
		numberOfSegmentsInElement2 	= 9	; order2 = 10	; whichCase = 2; segLabel = 6; 
		#include "callRight.c"
		} else { // general case.
			segLabel = -1;
			// right branch. This is where you add O10 elements. To do this, you have to know the O11 coordinates - the conditions are hard coded.
			#include "generateNumSegs.c"
			if(order2<stopOrder)   tNode->right = NULL; 	
			else{ 	
			// apparently the # has to be the first character on the line.
			#include "callRight.c" 	
			}
		}


	
}
// 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
// 1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
else if(numberOfSegmentsInElement == 0 && order>stopOrder){ // add one subtree at the end of this node on the right branch.

			#include "endNodes.c"

			// coordinates must be assigned to this node before the following new nodes are generated.
			/* left node is where you attach a suitable new element. */
			#include "generateNumSegs2.c"
			if(order2<stopOrder)   tNode->left = NULL; 
			else
			{
			#include "callLeft2.c"			
			}
			/* right node is where you attach a suitable new element. */
			#include "generateNumSegs2.c"
			if(order2<stopOrder)   tNode->right = NULL; 
			else
			{
			#include "callRight.c"			
			}
}
// 222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
// 222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
// 222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
else if(numberOfSegmentsInElement == 0 && order==stopOrder){ // when number of segments is 0, just make a node and return. done.

// no coordinate assignment here.
	tNode->left 		=      NULL;
	tNode->right 	=      NULL;

} // end of numberOfSegmentsInElement == 0 && order==stopOrder
// 3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
// 3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
// 3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
// 3333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
	

return tNode;
} // end of makeKassabTree.


    /*Function which helps the print_path to recursively print all the nodes*/ 
    void print_paths_recurL(node *node, double path[], int path_len, int startOrder, int stopOrder){ // this needs to return a value to print paths.
      FILE *pathLength;
      int i, count, return_value;      
      double lengthofPath;
      return_value = -1;

      if (node == NULL){
	      return; 
      }
      
      if(node->strahler<=startOrder && node->strahler>=stopOrder)
      path[path_len] = node->length;
      else
      path[path_len] = 0.0;
      
      path_len++;
      if (node->left == NULL && node->right == NULL) {
	      lengthofPath = 0.0;
      for (i = 0; i < path_len; i++) lengthofPath = lengthofPath + path[i];
//      printf("Length of path %f between orders %d and %d\n", lengthofPath, stopOrder, startOrder);
	pathLength = fopen("pathLength.data","a+");
	fprintf(pathLength, "%f\n", lengthofPath);
	fclose(pathLength);
      }
      else{
        print_paths_recurL(node->left,    path, path_len, startOrder, stopOrder);    //recursively calls the left node of the tree
        print_paths_recurL(node->right, path, path_len, startOrder, stopOrder);    //recursively calls the right node of the tree
      }
    } // end of print_paths_recur.

    /*Function to store all the paths from the root node to all leaf nodes in  a array*/
    void print_pathsL(node *node, int startOrder, int stopOrder){  // this needs to return a value to the main. 0 if n1 and n2 are parents
      double path[100000];
      print_paths_recurL(node, path, 0, startOrder, stopOrder);
    }

    
void onOffAssignedRecurse(node *root){

if(root->left!=NULL && root->right!=NULL){
	if(root->assignedOrNot==1 && root->onOff==1)	{ root->left->onOff = 1; root->right->onOff = 1;   }
	else 											{ root->left->onOff = 0; root->right->onOff = 0; }
	onOffAssignedRecurse(root->left);    onOffAssignedRecurse(root->right);
}
else if(root->left!=NULL && root->right==NULL){
	if(root->assignedOrNot==1 && root->onOff==1)	{ root->left->onOff = 1;  }
	else 									{ root->left->onOff = 0; }
	onOffAssignedRecurse(root->left);    
}
else if(root->left==NULL && root->right!=NULL){
	if(root->assignedOrNot==1 && root->onOff==1)	{ root->right->onOff = 1;   }
	else 									{ root->right->onOff = 0; }
	onOffAssignedRecurse(root->right);
}
else if(root->left==NULL && root->right==NULL){ }


} // end of onOffAssignedRecurse.

//***********************************************************************************************************************************************
//***********************************************************************************************************************************************
//***********************************************************************************************************************************************

void boundaryVector(XYZ xu_vector, XYZ nusunit_vector, XYZ *direction_vector, double fractionforOrder, double Ls, int myStrahler, int parentStrahler){

XYZ  avbdr_closest, avbdr_normal, LVCHAMBER_closest, LVCHAMBER_normal, RVCHAMBER_closest, RVCHAMBER_normal, OUTSIDE_closest, OUTSIDE_normal, wns_vector;
XYZ projection;
double distToSurface_avbdr, distToSurface_LVCHAMBER, distToSurface_RVCHAMBER, distToSurface_OUTSIDE, nus_length;

wns_vector.x = wns_vector.y = wns_vector.z = -1.0;

// on the AV border plane.
distToSurface_avbdr 	= closestAVBDR(xu_vector, &avbdr_closest, &avbdr_normal, fractionforOrder);
// on the inside LV cavity surface. This is not an inequality, but should be made so for the discrete case testing of distances.
distToSurface_LVCHAMBER 	= closestLVCHAMBER(xu_vector, &LVCHAMBER_closest, &LVCHAMBER_normal, fractionforOrder);
// RV cavity. It is a bit of LV outside wall and the RV inside wall.
distToSurface_RVCHAMBER 	= closestRVCHAMBER(xu_vector, &RVCHAMBER_closest, &RVCHAMBER_normal, fractionforOrder);
// Outside.
distToSurface_OUTSIDE = closestOUTSIDE(xu_vector, &OUTSIDE_closest, &OUTSIDE_normal, fractionforOrder);
	
//
nus_length=sqrt( avbdr_normal.x * avbdr_normal.x + avbdr_normal.y * avbdr_normal.y + avbdr_normal.z * avbdr_normal.z );
avbdr_normal.x = avbdr_normal.x/nus_length; avbdr_normal.y = avbdr_normal.y/nus_length; avbdr_normal.z = avbdr_normal.z/nus_length;
nus_length= sqrt(LVCHAMBER_normal.x*LVCHAMBER_normal.x+LVCHAMBER_normal.y*LVCHAMBER_normal.y+LVCHAMBER_normal.z*LVCHAMBER_normal.z);
LVCHAMBER_normal.x=LVCHAMBER_normal.x/nus_length; LVCHAMBER_normal.y=LVCHAMBER_normal.y/nus_length; LVCHAMBER_normal.z= LVCHAMBER_normal.z/nus_length;
nus_length=sqrt( RVCHAMBER_normal.x*RVCHAMBER_normal.x+RVCHAMBER_normal.y*RVCHAMBER_normal.y+RVCHAMBER_normal.z*RVCHAMBER_normal.z);
RVCHAMBER_normal.x=RVCHAMBER_normal.x/nus_length; RVCHAMBER_normal.y=RVCHAMBER_normal.y/nus_length; RVCHAMBER_normal.z= RVCHAMBER_normal.z/nus_length;
nus_length 			= sqrt( OUTSIDE_normal.x * OUTSIDE_normal.x + OUTSIDE_normal.y * OUTSIDE_normal.y + OUTSIDE_normal.z * OUTSIDE_normal.z );
OUTSIDE_normal.x = OUTSIDE_normal.x/nus_length; OUTSIDE_normal.y = OUTSIDE_normal.y/nus_length; OUTSIDE_normal.z = OUTSIDE_normal.z/nus_length;

// the weighted sums of the normal vectors, eqn 3 of page 286.
wns_vector.x = avbdr_normal.x * exp(-distToSurface_avbdr/(2.0 * Ls) ) 				+ LVCHAMBER_normal.x * exp(-distToSurface_LVCHAMBER/(2.0 * Ls) ) 
			+ RVCHAMBER_normal.x * exp(-distToSurface_RVCHAMBER/(2.0 * Ls) ) 	+ OUTSIDE_normal.x * exp(-distToSurface_OUTSIDE/(2.0 * Ls) )			;
wns_vector.y = avbdr_normal.y * exp(-distToSurface_avbdr/(2.0 * Ls) ) 				+ LVCHAMBER_normal.y * exp(-distToSurface_LVCHAMBER/(2.0 * Ls) ) 
			+ RVCHAMBER_normal.y * exp(-distToSurface_RVCHAMBER/(2.0 * Ls) ) 	+ OUTSIDE_normal.y * exp(-distToSurface_OUTSIDE/(2.0 * Ls) )			;
wns_vector.z = avbdr_normal.z * exp(-distToSurface_avbdr/(2.0 * Ls) ) 				+ LVCHAMBER_normal.z * exp(-distToSurface_LVCHAMBER/(2.0 * Ls) ) 
			+ RVCHAMBER_normal.z * exp(-distToSurface_RVCHAMBER/(2.0 * Ls) ) 	+ OUTSIDE_normal.z * exp(-distToSurface_OUTSIDE/(2.0 * Ls) )			;

nus_length = sqrt(wns_vector.x*wns_vector.x+wns_vector.y*wns_vector.y+wns_vector.z*wns_vector.z);
if(nus_length > 0.0){
wns_vector.x = wns_vector.x / nus_length; wns_vector.y = wns_vector.y / nus_length;  wns_vector.z = wns_vector.z / nus_length;
}
else
{
wns_vector.x = 0.0; wns_vector.y = 0.0;  wns_vector.z = 1.0;	// bug, but for now.
}

			if(debug){
			printf("xu vector: %f %f %f\n", xu_vector.x, xu_vector.y, xu_vector.z);
			printf("wns_vector: %f %f %f %f %f %f\n", wns_vector.x, wns_vector.y, wns_vector.z, wns_vector.x+xu_vector.x, wns_vector.y+xu_vector.y, wns_vector.z+xu_vector.z);
			printf("nusunit_vector: %f %f %f %f\n", nusunit_vector.x, nusunit_vector.y, nusunit_vector.z, nusunit_vector.x * nusunit_vector.x+ nusunit_vector.y*nusunit_vector.y+ nusunit_vector.z* nusunit_vector.z);
			}
			
direction_vector->x = nusunit_vector.x - wns_vector.x; direction_vector->y = nusunit_vector.y - wns_vector.y; direction_vector->z = nusunit_vector.z - wns_vector.z;
nus_length = sqrt(direction_vector->x * direction_vector->x + direction_vector->y * direction_vector->y + direction_vector->z * direction_vector->z );
// normalised direction vector:
direction_vector->x = direction_vector->x / nus_length; direction_vector->y = direction_vector->y / nus_length;   direction_vector->z = direction_vector->z / nus_length;
			
// try the projection.
// projection of direction vector along wns_vector.
/*
if(myStrahler != 8){
projection.x = ( (direction_vector->x)*(wns_vector.x) + (direction_vector->y)*(wns_vector.y) + (direction_vector->z)*(wns_vector.z) ) * wns_vector.x;
projection.y = ( (direction_vector->x)*(wns_vector.x) + (direction_vector->y)*(wns_vector.y) + (direction_vector->z)*(wns_vector.z) ) * wns_vector.y;
projection.z = ( (direction_vector->x)*(wns_vector.x) + (direction_vector->y)*(wns_vector.y) + (direction_vector->z)*(wns_vector.z) ) * wns_vector.z;

// substract the projection from direction vector.
direction_vector->x = direction_vector->x - projection.x;
direction_vector->y = direction_vector->y - projection.y;
direction_vector->z = direction_vector->z - projection.z;

// normalise this.
nus_length = sqrt(direction_vector->x * direction_vector->x + direction_vector->y * direction_vector->y + direction_vector->z * direction_vector->z );
direction_vector->x = direction_vector->x / nus_length; direction_vector->y = direction_vector->y / nus_length;   direction_vector->z = direction_vector->z / nus_length;
}
else 
*/	
if(myStrahler==8 && parentStrahler>8){
	direction_vector->x = -wns_vector.x; direction_vector->y = -wns_vector.y; direction_vector->z = -wns_vector.z;
}

} // end of direction vector.

//***********************************************************************************************************************************************
//***********************************************************************************************************************************************
//***********************************************************************************************************************************************

void supplyVector(XYZ xu_vector, XYZ *nus_vector, node *tNode, int my_order, int my_data, double Ls){

	XYZ di_vector, yihat_vector;
	double di_length, nus_length, temp_termLsdi;
	
	if(tNode->data!=my_data && my_order <= tNode->strahler && tNode->assignedOrNot==1){
		di_vector.x = xu_vector.x - tNode->x_cood;
		di_vector.y = xu_vector.y - tNode->y_cood;
		di_vector.z = xu_vector.z - tNode->z_cood;		
		
	di_length   = sqrt ( di_vector.x * di_vector.x + di_vector.y * di_vector.y + di_vector.z * di_vector.z ); // di_length in principle could end up being 0.
//	if(di_length<=0) { printf("The length is unrealistic or negative, exiting.\n"); exit(0); }
	if(di_length>0.0){
	// normalise the vector, this is the y_i^ on page 285, para 1 of Beard and Bassingthwaighte.
	yihat_vector.x = di_vector.x/di_length;  yihat_vector.y = di_vector.y/di_length; yihat_vector.z = di_vector.z/di_length;
	temp_termLsdi = pow( Ls / di_length, zeta ) / (1.0 + pow( Ls / di_length, zeta ) );
		nus_vector->x = nus_vector->x  + temp_termLsdi * yihat_vector.x;
		nus_vector->y = nus_vector->y + temp_termLsdi * yihat_vector.y;
		nus_vector->z = nus_vector->z  + temp_termLsdi * yihat_vector.z;
	}
	else
	{
		nus_vector->x = nus_vector->x  + 0.0001; // hack.
		nus_vector->y = nus_vector->y + 0.0001;
		nus_vector->z = nus_vector->z  + 0.0001;	
	}



		if(debug)
	printf("%f %f %f\n", Ls, di_length, temp_termLsdi );
	} // end of nus vector construction.

	if(tNode->left!=NULL&tNode->right!=NULL){
		supplyVector(xu_vector, nus_vector, tNode->left, my_order, my_data, Ls);
		supplyVector(xu_vector, nus_vector, tNode->right, my_order, my_data, Ls);
	}
	else if(tNode->left!=NULL&tNode->right==NULL){
		supplyVector(xu_vector, nus_vector, tNode->left, my_order, my_data, Ls);
	}
	else if(tNode->left==NULL&tNode->right!=NULL){
		supplyVector(xu_vector, nus_vector, tNode->left, my_order, my_data, Ls);
	}
	else if(tNode->left==NULL&tNode->right==NULL){}

} // end of supply vector function.


//***********************************************************************************************************************************************
//***********************************************************************************************************************************************
//***********************************************************************************************************************************************

/* Inputs to set coordinates of tree: which node to consider, root of RCA, root of LA.
   By going to me node, you are setting up coordinates of left node and right node.
*/

void setCoordinatesofTree(node *tNode, node *rootRCA, node *rootLA, XYZ xu_vectorp, int stopOrder, int thisOrder){

XYZ 	     xu_vector, nus_vector, nusunit_vector, direction_vector, new_loc, sp, vd, vd1, vd2, nb;
XYZ 	     avbdr_closest, avbdr_normal, LVCHAMBER_closest, LVCHAMBER_normal, RVCHAMBER_closest, RVCHAMBER_normal, OUTSIDE_closest, OUTSIDE_normal;
double 	nus_length, di_length, cs, cb, spdotvd, ux, uy, uz, thetaT, sint, cost, Ls;
double 	rotation[3][3], distToSurface_avbdr, distToSurface_LVCHAMBER, distToSurface_RVCHAMBER, distToSurface_OUTSIDE, fractionforOrder;
double 	epsilon; 
int 		stopIter, doLeft, doRight, maxStopIter, dontDoIt;
FILE 	*output;

stopIter = 0; epsilon = 0.1; // mm.

// by the time you get here, the node is definitely defined with assignedOrNot = 1.
if(tNode->strahler>=stopOrder){
	if(tNode->whichCase==1){
output = fopen("rawSegsRCA.data","a+");
fprintf(output, "%50.50lf %50.50lf %50.50lf %d %d %50.50lf %50.50lf %50.50lf %50.50lf %d %d\n", tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->data, tNode->parent_data, tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->strahler, tNode->parent_srahler ); // the last one was thisOrder in pre- 22 Sept versions.
fclose(output);
	}
	else{
output = fopen("rawSegsLA.data","a+");
fprintf(output, "%50.50lf %50.50lf %50.50lf %d %d %50.50lf %50.50lf %50.50lf %50.50lf %d %d\n", tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->data, tNode->parent_data, tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->strahler, tNode->parent_srahler ); // the last one was thisOrder in pre- 22 Sept versions.
fclose(output);
	}		
}

// go down tree until you have left and right kids who have thisOrder or more.

// common parts, do this each time. *********************************************************************************************************************
dontDoIt = 1;
// get supply vector, xu is now me node.
xu_vector.x = tNode->x_cood; xu_vector.y = tNode->y_cood; xu_vector.z = tNode->z_cood;
nus_vector.x = nus_vector.y = nus_vector.z = 0.0;

// Fung sp vector, common to both branches.
// sp is the vector of parent segment from the end of which the children start.
if(tNode->parent_data<0){
sp.x = 1000000.0; sp.y = 1000.0; sp.z = -2.34;	// this does not matter, it is never used.
dontDoIt = 0;
}
else
{
sp.x = xu_vector.x - xu_vectorp.x; sp.y = xu_vector.y - xu_vectorp.y; sp.z = xu_vector.z - xu_vectorp.z;
dontDoIt = 1;
}

if(debug)
printf("XU and XUP, %f %f %f %f %f %f, %d\n", xu_vector.x , xu_vectorp.x , xu_vector.y , xu_vectorp.y, xu_vector.z, xu_vectorp.z, tNode->whichCase);
// normalise sp.
di_length = sqrt ( sp.x * sp.x + sp.y * sp.y + sp.z * sp.z );
// normalise if di_length is non-zero finite or exit.
if(di_length>0){ // do the solution only if di_length is okay. this is a big if.
	sp.x = sp.x /di_length;  sp.y = sp.y / di_length; sp.z = sp.z / di_length;
	dontDoIt = 1;
}
else
{
	sp.x = 1.0; sp.y = 1.0; sp.y = 1.0; // hack.
dontDoIt = 0;
}

//RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
//========================================================================================================
// right child.
if(tNode->right!=NULL && tNode->right->strahler>=stopOrder && dontDoIt == 1){
	
new_loc.x = -1.0; new_loc.y = -1.0; new_loc.z = -1.0;
if(tNode->right->assignedOrNot==1){	
	// maybe the right coordinates should be written and not counted as lost? This holds, by definition.

	xu_vectorp.x = xu_vector.x;  xu_vectorp.y = xu_vector.y; xu_vectorp.z = xu_vector.z; // the right solution depends on this xu_vector. So this is not in the if.	
	// call set coordinates right.
	setCoordinatesofTree(tNode->right, rootRCA, rootLA, xu_vectorp, stopOrder, thisOrder);	
	
	
} else if(tNode->right->assignedOrNot!=1){
// set your fraction for order.
// if(tNode->right->strahler>=9) fractionforOrder = 0.10; else fractionforOrder = 1.0;
if(tNode->right->strahler>=9) fractionforOrder = 0.33; else fractionforOrder = 1.0; // make it a bit more linienant.

Ls = tNode->right->length;
nus_vector.x = nus_vector.y = nus_vector.z = 0.0;
supplyVector(xu_vector, &nus_vector, rootRCA, tNode->strahler, tNode->data, Ls);
supplyVector(xu_vector, &nus_vector, rootLA,    tNode->strahler, tNode->data, Ls);
// align it along the space filling, normalise the nus_vector.
	nus_length = sqrt( nus_vector.x * nus_vector.x + nus_vector.y * nus_vector.y + nus_vector.z * nus_vector.z );		
	if(nus_length>0.0){ nusunit_vector.x = nus_vector.x / nus_length;  nusunit_vector.y = nus_vector.y / nus_length; nusunit_vector.z = nus_vector.z / nus_length;  }
	else 			{ printf("The nus_length was 0, exiting. %f\n", nus_length); /* exit(0); */														         }

	// get boundary vector.
	// the vector is the composite direction vector.
	boundaryVector(xu_vector, nusunit_vector, &direction_vector, fractionforOrder, tNode->length, tNode->strahler, tNode->parent_srahler);
// 	Beard only new coordinates.
//	new_loc.x=xu_vector.x+Ls*direction_vector.x; new_loc.y=xu_vector.y+Ls*direction_vector.y; new_loc.z=xu_vector.z+Ls*direction_vector.z;

// Put in Fung theta1 theta2.
// Fung et al. method of Branching plane computation. Phy Med Biol 2011. fung, 56(17): 5651-5663.
// the cs and cb are to be estimated based on where in the 3D space you are.

if(tNode->strahler>=9)	{  cs = 0.1; cb = 0.9;  }
else 					{ cs = 0.5; cb = 0.5;  }

do{
vd.x = cs * nusunit_vector.x + cb * direction_vector.x;    vd.y = cs * nusunit_vector.y + cb * direction_vector.y;   vd.z = cs * nusunit_vector.z + cb * direction_vector.z;
// now get your nb vector with vd and sp.
spdotvd = (vd.x)*(sp.x)+(vd.y)*(sp.y)+(vd.z)*(sp.z); // sp dot vd.
nb.x = vd.x * spdotvd - sp.x; nb.y = vd.y * spdotvd - sp.y;  nb.z = vd.z * spdotvd - sp.z; 
if(fabs(spdotvd-1.0) < 0.001) {cs = genrand_real1(); cb = 1.0 - cs; }
}while(fabs(spdotvd-1.0) < 0.001);

// for the sake of readability, copy it to another vector.
ux = nb.x; uy = nb.y; uz = nb.z;

// left segment.
thetaT = tNode->theta2; cost = cos(thetaT); sint = sin(thetaT); 
if(debug)
printf("theta: %f\n", tNode->theta2);
// now you need the theta values for the following. second variable is row of matrix, 1st variable is column.
rotation[0][0] = cost+ux*ux*(1.0-cost)	    ;  rotation[0][1] = ux*uy*(1.0-cost)-uz*sint; 	rotation[0][2] = ux*uz*(1.0-cost)+uy*sint; 
rotation[1][0]  = uy*ux*(1.0-cost)+uz*sint ;  rotation[1][1]  = cost+uy*uy*(1.0-cost); 		rotation[1][2]  = uy*uz*(1.0-cost)-ux*cost; 
rotation[2][0] = uz*ux*(1.0-cost)-uy*sint ;  rotation[2][1]  = uz*uy*(1.0-cost)+ux*sint; 	rotation[2][2] = cost+uz*uz*(1.0-cost); 

vd1.x = rotation[0][0] * vd.x + rotation[0][1] * vd.y + rotation[0][2] * vd.z;
vd1.y = rotation[1][0] * vd.x + rotation[1][1] * vd.y + rotation[1][2] * vd.z;
vd1.z = rotation[2][0] * vd.x + rotation[2][1] * vd.y + rotation[2][2] * vd.z;

// optimise left, optimise right.
// left.
new_loc.x=xu_vector.x+(tNode->length)*vd1.x; new_loc.y=xu_vector.y+(tNode->length)*vd1.y; new_loc.z=xu_vector.z+(tNode->length)*vd1.z;

stopIter = 0;
if(tNode->right->strahler<7) maxStopIter = 50; else maxStopIter = 200;
if(isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder)!=1)
do{
	// work out closest point on closest surface. work out normal.
	distToSurface_avbdr 			= closestAVBDR		(new_loc, &avbdr_closest, &avbdr_normal, fractionforOrder);
	distToSurface_LVCHAMBER 	= closestLVCHAMBER	(new_loc, &LVCHAMBER_closest, &LVCHAMBER_normal, fractionforOrder);
	distToSurface_RVCHAMBER 	= closestRVCHAMBER	(new_loc, &RVCHAMBER_closest, &RVCHAMBER_normal, fractionforOrder);
	distToSurface_OUTSIDE 		= closestOUTSIDE	(new_loc, &OUTSIDE_closest, &OUTSIDE_normal, fractionforOrder);

	if(distToSurface_avbdr <= distToSurface_LVCHAMBER && distToSurface_avbdr <= distToSurface_RVCHAMBER && distToSurface_avbdr <= distToSurface_OUTSIDE ){
direction_vector.x = direction_vector.x   - epsilon * avbdr_normal.x; direction_vector.y = direction_vector.y  - epsilon * avbdr_normal.y; direction_vector.z = direction_vector.z   - epsilon * avbdr_normal.z;
	}
	else if(distToSurface_LVCHAMBER <=distToSurface_avbdr&&distToSurface_LVCHAMBER <=distToSurface_RVCHAMBER&&distToSurface_LVCHAMBER <=distToSurface_OUTSIDE){
direction_vector.x = direction_vector.x   - epsilon * LVCHAMBER_normal.x; direction_vector.y = direction_vector.y  - epsilon * LVCHAMBER_normal.y;
direction_vector.z = direction_vector.z   - epsilon * LVCHAMBER_normal.z;
	}
	else if(distToSurface_RVCHAMBER<=distToSurface_avbdr&&distToSurface_RVCHAMBER<=distToSurface_LVCHAMBER
	&&distToSurface_RVCHAMBER<=distToSurface_OUTSIDE ){
direction_vector.x = direction_vector.x   - epsilon * RVCHAMBER_normal.x; direction_vector.y = direction_vector.y  - epsilon * RVCHAMBER_normal.y; direction_vector.z = direction_vector.z   - epsilon * RVCHAMBER_normal.z;
	}
	else{
direction_vector.x = direction_vector.x   - epsilon * OUTSIDE_normal.x; direction_vector.y = direction_vector.y  - epsilon * OUTSIDE_normal.y; direction_vector.z = direction_vector.z   - epsilon * OUTSIDE_normal.z;
	}

	nus_length = sqrt(direction_vector.x * direction_vector.x + direction_vector.y * direction_vector.y + direction_vector.z * direction_vector.z );
	// normalised direction vector:
	direction_vector.x = direction_vector.x / nus_length; direction_vector.y = direction_vector.y / nus_length; direction_vector.z = direction_vector.z / nus_length;

	new_loc.x=xu_vector.x+Ls*direction_vector.x; new_loc.y=xu_vector.y+Ls*direction_vector.y; new_loc.z=xu_vector.z+Ls*direction_vector.z;

	stopIter++;
}while(isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder)==0 && stopIter < maxStopIter);

	if( (isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, 0.33)==1) || (isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, 1.0)==1) ){
	printf("Right Optimised, in or out: %d, %d. radius: %f; strahler: %d; which sub-tree? %d, at %f %f %f. data: %d\n", isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder), stopIter, tNode->right->radius, tNode->right->strahler, tNode->right->whichCase, new_loc.x, new_loc.y, new_loc.z, tNode->right->data);
	xu_vectorp.x = xu_vector.x;  xu_vectorp.y = xu_vector.y; xu_vectorp.z = xu_vector.z;	
	tNode->right->x_cood = new_loc.x; tNode->right->y_cood = new_loc.y; tNode->right->z_cood = new_loc.z;  // this is left, not me - you calculated starting at me.
	tNode->right->assignedOrNot = 1; // switch this based on if the coods are in tissue or not.
	setCoordinatesofTree(tNode->right, rootRCA, rootLA, xu_vectorp, stopOrder, thisOrder);
} else {
	printf("Right Optimised, in or out: %d, %d. radius: %f; strahler: %d; which sub-tree? %d, at %f %f %f. data: %d\n", isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder), stopIter, tNode->right->radius, tNode->right->strahler, tNode->right->whichCase, new_loc.x, new_loc.y, new_loc.z, tNode->right->data);
	
		output = fopen("logLoss.data","a+");
		fprintf(output, "Right lost, %d %d %d %f %f %f %d %f %f %f %d %d\n", tNode->right->strahler, tNode->strahler, isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder), new_loc.x, new_loc.y, new_loc.z, tNode->right->whichCase, xu_vector.x, xu_vector.y, xu_vector.z, tNode->right->data, tNode->data);
		fclose(output);
	}

} // end of assigned or not.

} // end of tNode->right calculation.

//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
//**************************************************************************************************************************************************
if(tNode->left!=NULL && tNode->left->strahler>=stopOrder && dontDoIt == 1 ){ 
new_loc.x = -1.0; new_loc.y = -1.0; new_loc.z = -1.0;  // placeholder.
if( tNode->left->assignedOrNot==1 ){
	xu_vectorp.x = xu_vector.x;  xu_vectorp.y = xu_vector.y; xu_vectorp.z = xu_vector.z; // the right solution depends on this xu_vector. So this is not in the if.	
	// call set coordinates left.
	setCoordinatesofTree(tNode->left, rootRCA, rootLA, xu_vectorp, stopOrder, thisOrder);
}
else if(tNode->left->assignedOrNot!=1 ){
// set your fraction for order.
if(tNode->left->strahler>=9) fractionforOrder = 0.01; else fractionforOrder = 1.0;

Ls = tNode->left->length;
nus_vector.x = nus_vector.y = nus_vector.z = 0.0;
supplyVector(xu_vector, &nus_vector, rootRCA, tNode->strahler, tNode->data, Ls);
supplyVector(xu_vector, &nus_vector, rootLA,    tNode->strahler, tNode->data, Ls);
// align it along the space filling, normalise the nus_vector.
	nus_length = sqrt( nus_vector.x * nus_vector.x + nus_vector.y * nus_vector.y + nus_vector.z * nus_vector.z );		
	if(nus_length>0.0){ nusunit_vector.x = nus_vector.x / nus_length;  nusunit_vector.y = nus_vector.y / nus_length; nusunit_vector.z = nus_vector.z / nus_length;  }
	else 			{ printf("The nus_length was 0, exiting. %f\n", nus_length); /* exit(0); */ 															         }

	// get boundary vector.
	// the vector is the composite direction vector.
	boundaryVector(xu_vector, nusunit_vector, &direction_vector, fractionforOrder, tNode->length, tNode->strahler, tNode->parent_srahler);
// 	Beard only new coordinates.
	new_loc.x=xu_vector.x+Ls*direction_vector.x; new_loc.y=xu_vector.y+Ls*direction_vector.y; new_loc.z=xu_vector.z+Ls*direction_vector.z;

// Put in Fung theta1 theta2.
// Fung et al. method of Branching plane computation. Phy Med Biol 2011. fung, 56(17): 5651-5663.
// the cs and cb are to be estimated based on where in the 3D space you are.

if(tNode->strahler>=9)	{    cs = 0.1; cb = 0.9;    }
else 					{    cs = 0.5; cb = 0.5;   }

do{
vd.x = cs * nusunit_vector.x + cb * direction_vector.x;    vd.y = cs * nusunit_vector.y + cb * direction_vector.y;   vd.z = cs * nusunit_vector.z + cb * direction_vector.z;
// now get your nb vector with vd and sp.
spdotvd = (vd.x)*(sp.x)+(vd.y)*(sp.y)+(vd.z)*(sp.z); // sp dot vd.
nb.x = vd.x * spdotvd - sp.x; nb.y = vd.y * spdotvd - sp.y;  nb.z = vd.z * spdotvd - sp.z; 
if(fabs(spdotvd-1.0) < 0.001) {cs = genrand_real1(); cb = 1.0 - cs; }
}while(fabs(spdotvd-1.0) < 0.001);

// for the sake of readability, copy it to another vector.
ux = nb.x; uy = nb.y; uz = nb.z;

// left segment.
thetaT = tNode->theta1; cost = cos(thetaT); sint = sin(thetaT); 
if(debug)
printf("theta: %f\n", tNode->theta1);
// now you need the theta values for the following. second variable is row of matrix, 1st variable is column.
rotation[0][0] = cost+ux*ux*(1.0-cost)	    ;  rotation[0][1] = ux*uy*(1.0-cost)-uz*sint; 	rotation[0][2] = ux*uz*(1.0-cost)+uy*sint; 
rotation[1][0]  = uy*ux*(1.0-cost)+uz*sint ;  rotation[1][1]  = cost+uy*uy*(1.0-cost); 		rotation[1][2]  = uy*uz*(1.0-cost)-ux*cost; 
rotation[2][0] = uz*ux*(1.0-cost)-uy*sint ;  rotation[2][1]  = uz*uy*(1.0-cost)+ux*sint; 	rotation[2][2] = cost+uz*uz*(1.0-cost); 

vd1.x = rotation[0][0] * vd.x + rotation[0][1] * vd.y + rotation[0][2] * vd.z;
vd1.y = rotation[1][0] * vd.x + rotation[1][1] * vd.y + rotation[1][2] * vd.z;
vd1.z = rotation[2][0] * vd.x + rotation[2][1] * vd.y + rotation[2][2] * vd.z;

// optimise left, optimise right.
// left.
new_loc.x=xu_vector.x+(tNode->length)*vd1.x; new_loc.y=xu_vector.y+(tNode->length)*vd1.y; new_loc.z=xu_vector.z+(tNode->length)*vd1.z;

stopIter = 0;

// count how long it takes to do the do loop optimisation.
clock_t begin = clock();
if(tNode->left->strahler<7) maxStopIter = 50; else maxStopIter = 200;
if(isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder)!=1)
do{
	// work out closest point on closest surface. work out normal.
	distToSurface_avbdr 			= closestAVBDR		(new_loc, &avbdr_closest, &avbdr_normal, fractionforOrder);
	distToSurface_LVCHAMBER 	= closestLVCHAMBER	(new_loc, &LVCHAMBER_closest, &LVCHAMBER_normal, fractionforOrder);
	distToSurface_RVCHAMBER 	= closestRVCHAMBER	(new_loc, &RVCHAMBER_closest, &RVCHAMBER_normal, fractionforOrder);
	distToSurface_OUTSIDE 		= closestOUTSIDE	(new_loc, &OUTSIDE_closest, &OUTSIDE_normal, fractionforOrder);

	if(distToSurface_avbdr <= distToSurface_LVCHAMBER && distToSurface_avbdr <= distToSurface_RVCHAMBER && distToSurface_avbdr <= distToSurface_OUTSIDE ){
direction_vector.x = direction_vector.x   - epsilon * avbdr_normal.x; direction_vector.y = direction_vector.y  - epsilon * avbdr_normal.y; direction_vector.z = direction_vector.z   - epsilon * avbdr_normal.z;
	}
	else if(distToSurface_LVCHAMBER <=distToSurface_avbdr&&distToSurface_LVCHAMBER <=distToSurface_RVCHAMBER&&distToSurface_LVCHAMBER <=distToSurface_OUTSIDE){
direction_vector.x = direction_vector.x   - epsilon * LVCHAMBER_normal.x; direction_vector.y = direction_vector.y  - epsilon * LVCHAMBER_normal.y;
direction_vector.z = direction_vector.z   - epsilon * LVCHAMBER_normal.z;
	}
	else if(distToSurface_RVCHAMBER<=distToSurface_avbdr&&distToSurface_RVCHAMBER<=distToSurface_LVCHAMBER
	&&distToSurface_RVCHAMBER<=distToSurface_OUTSIDE ){
direction_vector.x = direction_vector.x   - epsilon * RVCHAMBER_normal.x; direction_vector.y = direction_vector.y  - epsilon * RVCHAMBER_normal.y; direction_vector.z = direction_vector.z   - epsilon * RVCHAMBER_normal.z;
	}
	else{
direction_vector.x = direction_vector.x   - epsilon * OUTSIDE_normal.x; direction_vector.y = direction_vector.y  - epsilon * OUTSIDE_normal.y; direction_vector.z = direction_vector.z   - epsilon * OUTSIDE_normal.z;
	}

	nus_length = sqrt(direction_vector.x * direction_vector.x + direction_vector.y * direction_vector.y + direction_vector.z * direction_vector.z );
	// normalised direction vector:
	direction_vector.x = direction_vector.x / nus_length; direction_vector.y = direction_vector.y / nus_length; direction_vector.z = direction_vector.z / nus_length;

	new_loc.x=xu_vector.x+Ls*direction_vector.x; new_loc.y=xu_vector.y+Ls*direction_vector.y; new_loc.z=xu_vector.z+Ls*direction_vector.z;

	stopIter++;
}while(isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder)==0 && stopIter < maxStopIter);

clock_t end = clock();
double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
// printf("The optimisation took %f seconds.\n", time_spent);
output = fopen("timeTaken.data","a+");
fprintf(output, "%d %f %d %d\n", tNode->left->strahler, time_spent, stopIter, isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder));
fclose(output);

	// recursion along left branch.	
	// now just assign the coordinates on or off.
	if( (isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, 0.33)==1) || (isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, 1.0)==1) ){
	xu_vectorp.x = xu_vector.x;  xu_vectorp.y = xu_vector.y; xu_vectorp.z = xu_vector.z; // the right solution depends on this xu_vector. So this is not in the if.		
	printf("Left Optimised, in or out: %d, %d. radius: %f; strahler: %d; which sub-tree? %d, at %f %f %f. data: %d\n", isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder), stopIter, tNode->left->radius, tNode->left->strahler, tNode->left->whichCase, new_loc.x, new_loc.y, new_loc.z, tNode->left->data);
	
	tNode->left->x_cood = new_loc.x; tNode->left->y_cood = new_loc.y; tNode->left->z_cood = new_loc.z;  // this is left, not me - you calculated starting at me.
	tNode->left->assignedOrNot = 1; // switch this based on if the coods are in tissue or not.	
	setCoordinatesofTree(tNode->left, rootRCA, rootLA, xu_vectorp, stopOrder, thisOrder);
	}
	else{
	printf("Left Optimised, in or out: %d, %d. radius: %f; strahler: %d; which sub-tree? %d, at %f %f %f. data: %d\n", isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder), stopIter, tNode->left->radius, tNode->left->strahler, tNode->left->whichCase, new_loc.x, new_loc.y, new_loc.z, tNode->left->data);
	
		output = fopen("logLoss.data","a+");
		fprintf(output, "Left lost, %d %d %d %f %f %f %d %f %f %f %d %d\n", tNode->left->strahler, tNode->strahler, isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, fractionforOrder), new_loc.x, new_loc.y, new_loc.z, tNode->left->whichCase, xu_vector.x, xu_vector.y, xu_vector.z, tNode->left->data, tNode->data);
		fclose(output);
	}


} // end of assigned or not.
//****************************************************************************************************************************************************************


} // end of tNode->left calculation.

} // end of set coordinates recursion.

//***********************************************************************************************************************************************
//***********************************************************************************************************************************************
//***********************************************************************************************************************************************

void writeTreeSoFar(node *tNode){

FILE *output;

if(tNode->assignedOrNot==1 ){
	if(tNode->whichCase==2 || tNode->whichCase==3 )
output = fopen("rawSegsLA.data","a+");
else
	output = fopen("rawSegsRCA.data","a+");
	
fprintf(output, "%50.50lf %50.50lf %50.50lf %d %d %50.50lf %50.50lf %50.50lf %50.50lf %d %d\n", tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->data, tNode->parent_data, tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->strahler, tNode->parent_srahler ); 
fclose(output);

// go to left and right only if me is already assigned.
if(tNode->left!=NULL) writeTreeSoFar(tNode->left);
if(tNode->right!=NULL) writeTreeSoFar(tNode->right);


}
else
	return;
	
}

void mainNodesNotDefined(node *tNode){

	if(tNode->whichCase==1){
	if(tNode->strahler==11 && tNode->assignedOrNot!=1) printf("Node in RCA undefined: %d %f %f %f\n", tNode->strahler, tNode->x_cood, tNode->y_cood, tNode->z_cood);
	if(tNode->left!=NULL && tNode->left->strahler==11) mainNodesNotDefined(tNode->left);
	if(tNode->right!=NULL && tNode->right->strahler==11) mainNodesNotDefined(tNode->right);
	}

	if(tNode->whichCase==2){
	if(tNode->strahler==11 && tNode->assignedOrNot!=1) printf("Node in LAD undefined: %d %f %f %f\n", tNode->strahler, tNode->x_cood, tNode->y_cood, tNode->z_cood);
	if(tNode->left!=NULL && tNode->left->strahler==11) mainNodesNotDefined(tNode->left);
	if(tNode->right!=NULL && tNode->right->strahler==11) mainNodesNotDefined(tNode->right);
	
	if(tNode->left!=NULL && tNode->left->whichCase==3 && tNode->left->strahler==10) mainNodesNotDefined(tNode->left);
	if(tNode->right!=NULL && tNode->right->whichCase==3 && tNode->right->strahler==10) mainNodesNotDefined(tNode->right);
	}

	if(tNode->whichCase==3){
	if(tNode->strahler==10 && tNode->assignedOrNot!=1) printf("Node in LCX undefined: %d %f %f %f\n", tNode->strahler, tNode->x_cood, tNode->y_cood, tNode->z_cood);
	if(tNode->left!=NULL && tNode->left->strahler==10) mainNodesNotDefined(tNode->left);
	if(tNode->right!=NULL && tNode->right->strahler==10) mainNodesNotDefined(tNode->right);
	
	}

}

void pruneTree(node *tNode){

	FILE *output;
	
	if(tNode!=NULL && tNode->assignedOrNot==1){
/*		This worked at test. */

		if(tNode->whichCase==1){
	output = fopen("PrunedrawSegsRCA.data","a+");
	fprintf(output, "%50.50lf %50.50lf %50.50lf %d %d %50.50lf %50.50lf %50.50lf %50.50lf %d %d\n", tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->data, tNode->parent_data, tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->strahler, tNode->parent_srahler ); // the last one was thisOrder in pre- 22 Sept versions.
	fclose(output);
		} else {
	output = fopen("PrunedrawSegsLA.data","a+");
	fprintf(output, "%50.50lf %50.50lf %50.50lf %d %d %50.50lf %50.50lf %50.50lf %50.50lf %d %d\n", tNode->x_cood, tNode->y_cood, tNode->z_cood, tNode->data, tNode->parent_data, tNode->radius, tNode->resistanceSegment, tNode->pressure, tNode->perfusion, tNode->strahler, tNode->parent_srahler ); // the last one was thisOrder in pre- 22 Sept versions.
	fclose(output);
		}

		if(tNode->left!=NULL)    pruneTree(tNode->left);
		if(tNode->right!=NULL) pruneTree(tNode->right);
	} else {
		deleteTree(tNode);
		tNode = NULL;
	}
}

void pruneTree2(node *tNode){
// only tNode that is 1 goes here.	
if(tNode->left==NULL){ return; }
else
if(tNode->left!=NULL && tNode->left->assignedOrNot!=1)		{ deleteTree(tNode->left); 		tNode->left = NULL; 	} 
else if(tNode->left!=NULL && tNode->left->assignedOrNot==1)  
	pruneTree2(tNode->left)		;

if(tNode->right==NULL){ return; }
else
if(tNode->right!=NULL && tNode->right->assignedOrNot!=1) { deleteTree(tNode->right); 	tNode->right = NULL; 	} 
else if(tNode->right!=NULL && tNode->right->assignedOrNot==1)
	pruneTree2(tNode->right)	;

}



void reviseLengthResistance(node *tNode, double x0, double y0, double z0){
	
	if(tNode==NULL) return;
	
	if(tNode->data>-1 && tNode->assignedOrNot==1){
	tNode->length = sqrt( (tNode->x_cood-x0)*(tNode->x_cood-x0)+(tNode->y_cood-y0)*(tNode->y_cood-y0)+(tNode->z_cood-z0)*(tNode->z_cood-z0) );
	tNode->resistanceSegment 	= 8.0 * mu_0 * (tNode->length)/ ( M_PI * (tNode->radius) * (tNode->radius) * (tNode->radius) * (tNode->radius) );	
	x0 = tNode->x_cood; y0 = tNode->y_cood; z0 = tNode->z_cood;
	if(tNode->left!=NULL && tNode->left->assignedOrNot==1)  reviseLengthResistance(tNode->left, x0, y0, z0);
	if(tNode->right!=NULL && tNode->right->assignedOrNot==1)  reviseLengthResistance(tNode->right, x0, y0, z0);
	}

}


void readNodes(node *tNode){

FILE *output;
int i;

// for the viz.
int    viz_temp_x, viz_srkNodes0, viz_srkNodes1, viz_ThisNodeIsALeaf, viz_srkStrahlerNumberOfNode;
double viz_coodsMe0, viz_coodsMe1, viz_coodsMe2, viz_coodsParent0, viz_coodsParent1, viz_coodsParent2, viz_radiusSegment;
double viz_resistivity, viz_pressure, viz_perfusion, viz_downResistance;
double viz_theta1, viz_theta2, viz_theta, viz_leng, viz_xc, viz_yc, viz_zc;
int viz_assignedOrNot, viz_onOff, viz_whichCase;
double tx, ty, tz; // some temp coordinates.
int intp1, intp2, intme, intparent, intc1, intc2, intme_other, intparent_other;
double mean, sigma, U1, U2, radius, tRad;
int order, temp_order;

intp1 = intp2 = 0;
printf("out: %d %d %d\n", intp1, tNode->data, intp2);

output = fopen("ALL_98/NodeData.data.1","r");
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF){

	if(intp1==tNode->data){
	
		printf("%d %d %d\n", intp1, tNode->data, intp2);
	tNode->x_cood = tx;
	tNode->y_cood = ty;
	tNode->z_cood = tz;
	
//	tNode->data 	= intp1;
	tNode->parent_data 	= intp2;
	
	tNode->radius		= tRad;
	
	tNode->resistance 	= viz_resistivity;	
	tNode->pressure 		= viz_pressure;		
	tNode->perfusion 	= viz_perfusion;
	tNode->strahler 	= order;	
	tNode->parent_srahler = temp_order;
	break;
	}
	
}
fclose(output);

}

void doLengths(node *tNode){

// it always starts at root, which is always 1.
double xp, yp, zp, xc, yc, zc, leng;

xp = tNode->x_cood; yp = tNode->y_cood; zp = tNode->z_cood;

		if(tNode->left!=NULL && tNode->left->assignedOrNot==1){
		xc = tNode->left->x_cood; yc = tNode->left->y_cood; zc = tNode->left->z_cood;
		leng = sqrt( (xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp)   );
		tNode->left->length = leng;
		doLengths(tNode->left);
		}

		if(tNode->right!=NULL && tNode->right->assignedOrNot==1){
		xc = tNode->right->x_cood; yc = tNode->right->y_cood; zc = tNode->right->z_cood;
		leng = sqrt( (xc-xp)*(xc-xp) + (yc-yp)*(yc-yp) + (zc-zp)*(zc-zp)   );
		tNode->right->length = leng;
		doLengths(tNode->right);
		}
}
