/* Declarations. */
FILE *output, *random, *strahler;
/* Program reduced to use of FD with colouring. */
//  Vec            u		    			; /* solution, residual vectors. You need the r for SNES */
//  AppCtx         user		    	; /* user-defined work context                           */
//  DM             da		    		;
//  PetscInt       tissue_type		;
//  PetscInt       mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z, usr_i, usr_j, usr_k, temp_x, temp_y, temp_z, temp_x1, temp_x2;
  int mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z, usr_i, usr_j, usr_k, temp_x, temp_y, temp_z, temp_x1, temp_x2, tissue_type;

  int mysize_x1, mysize_y1, mysize_z1,  mysize_x2, mysize_y2, mysize_z2, mysize_x3, mysize_y3, mysize_z3;
  
// create the tree.
/*
root is the tree object.
rootR is the rearranged tree object.
*/

double temperature;

double x0, y0, z0;

node *root, *rootR, *rootRCA, *rootLA;  root = NULL; rootR = NULL; rootRCA = NULL; rootLA = NULL;
  int useThisValue = 0;
  int theParent, srkHeight, hgt, isDefined, range_x;
  int mex, mey, mez;  
  double mexR, meyR, mezR, normR, sgn;
  double pT;  
  int rad_int;
  double rad_double;
  char 		str[1000]; // so you do not keep creating and destroying this.  
  double slope_m;
  double a1, a2, a3, a4, b1, b2, b3, b4, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4;
  double theta, small_q, small_r; // q is size in x direction of small rectangle in the latin hypercube sampling and r is size in y direction.
double tx, ty, tz; // some temp coordinates.
double dist, dist_old, small_dist;
int hold_leaf_number;

double CF_tissue_supply,           volume,           CF_tissue_volume,           CF_power;
double total_cost, total_cost_new;
double CF_tissue_supply_new, volume_new, CF_tissue_volume_new, CF_power_new;
double mb, Pout;
double q0, radius_V_mm, QN, qq;
double length, lengthTempMax, bCost;
int radius_vasc;
char command[50];
int nodeTemp, strahlerTemp;
int numberofChildren;
int locLeft, locRight, locMy, leftStrahler, rightStrahler, allRadiiDefined, allPressuresDefined, allPressuresDefined2, allLengthsDefined, radTemp, ifLeafTemp;
double randLeft, randRight, randTemp, resLeft, resRight, resMy, perfusionTemp;
int iter_cood, accepted_cood;
int tempMV;
int b_v;
double fL, fR;

int bigLoop;

int o1, o2, o3, o4;
double temp_r, temp_theta;

// rearrage the tree.
// variables in rearrange have a NewR or rearrange in them, most of the time.
// this is not in the loop.
int rearrange_n1, rearrange_n2; // randomly chosen nodes in rearrangement.
int val_n1, val_n2, val_p1, val_p2, inx1, inx2, inx3;
FILE *pTest, *pp1, *pp2, *pp3, *p_levell;
int count_value, count_max;
int *in_inorderRCA, *in_preorderRCA, *in_postorderRCA;
int *in_inorderLA, *in_preorderLA, *in_postorderLA;
int notALeaf;
int intp1, intp2, intme, intparent, intc1, intc2, intme_other, intparent_other;
int subTreeMoved = 0;

	int numberOfDo;

/* Declare the required nodes. Ideally, these need to be destroyed.                                                                                      */
int *ThisNodeIsALeaf, **srkNodes, **srkStrahlerNumberOfNode, **srkLevelOfNode, **coodsParent, **NewcoodsParent, *setUp;

//revised datatype from into to double.
double **coodsMe, **NewcoodsMe;
int *parent; // size of nuNodes.
double **boundingBox, **NewboundingBox, **NewRboundingBox; // nuNodes x 4 array.

double          **fractionNodes,         *radiusSegment; // this is the same size as Nodes but is 1D - it holds the radius of the segment between Me and my Parent.
double *lengthSegment, *NewlengthSegment, *NewRlengthSegment;
double **NewfractionNodes, *NewradiusSegment; // this is the same size as Nodes but is 1D - it holds the radius of the segment between Me and my Parent.
double *resistivity,                  *perfusion,          *pressure; // the resistance and perfusion are properties of segments. pressure belongs to each node.


double *Newresistivity, *Newperfusion, *Newpressure; // the resistance and perfusion are properties of segments. pressure belongs to each node.
double *downResistance, *NewdownResistance, *NewRdownResistance;

int *leaves; // this array has indices of the leaves from the coodMe array.
double ROOT[2]; // this is the root, it is a constant.

int start, end;

// this is for when we start moving sub-trees around.
// these arrays are created once.
// Construct arrays that will let us calculate the new cost.
// these are not all random, but depend on the old arrays.
int **NewRsrkNodes; // this is my list of the new nodes. it can just be new?
int **NewRsrkLevelOfNode;
int **NewRsrkStrahlerNumberOfNode;
int *NewRThisNodeIsALeaf;
int **NewRcoodsMe;
int **NewRcoodsParent;

double *NewRradiusSegment; // this is the same size as Nodes but is 1D - it holds the radius of the segment between Me and my Parent.
double *NewRresistivity, *NewRperfusion, *NewRpressure; // the resistance and perfusion are properties of segments. pressure belongs to each node.

double angle, beta, totalTreeResistance, totalTreeResistanceRCA, totalTreeResistanceLA;
double downrTemp;

// for the viz.
int    viz_temp_x, viz_srkNodes0, viz_srkNodes1, viz_ThisNodeIsALeaf, viz_srkStrahlerNumberOfNode;
double viz_coodsMe0, viz_coodsMe1, viz_coodsMe2, viz_coodsParent0, viz_coodsParent1, viz_coodsParent2, viz_radiusSegment;
double viz_resistivity, viz_pressure, viz_perfusion, viz_downResistance;
double viz_theta1, viz_theta2, viz_theta, viz_leng, viz_xc, viz_yc, viz_zc;
int viz_assignedOrNot, viz_onOff, viz_whichCase;
int viz_data, viz_parent_data, viz_strahler, viz_parent_strahler;
double viz_radius, viz_length, viz_resistance, viz_resistanceSegment;

int temp_me_data;

double p0_x, p0_y, p0_z, p1_x, p1_y; // this is just mpx, mpy, tx, ty
// segment 2, all other segments except for segmentj1 and segmentj2.
double p2_x, p2_y, p3_x, p3_y, p1_z, p2_z, p3_z;
double q1_x, q1_y, q1_z, q2_x, q2_y, q2_z, q3_x, q3_y, q3_z;
int p2idx, p3idx;
int collided;
double lgt, numerat, slopej, slopejpr, segA, segB, x_0, y_0, parTx, parTy, parT;
int intersected;

int cont, cn;
double x_min, x_max, y_min, y_max, z_min, z_max;

int numSegsSoFar, segmentj, terminals;
double mpx, mpy, mpz; // for many segments, this is an array or a struct.
int mecoodsIdx, parentcoodsIdx;
	double LV_x_c, LV_y_c, LV_z_c;
	double xtmp, ytmp, ztmp;

int numberofOrders; // = 11; // use this after it settles into the recursion.
int thisOrder;
double mean, sigma, U1, U2, radius, tRad;
int numberOfSegments;

double amx[MAXORDERNUM];

double segmentElementRatio_meanRCA[MAXORDERNUM]; 
double      segmentElementRatio_SDRCA[MAXORDERNUM]; 
double           elementLengths_meanRCA[MAXORDERNUM]; 
double                elementLengths_SDRCA[MAXORDERNUM]; 
double          segmentLengths_meanRCA[MAXORDERNUM]; 
double               segmentLengths_SDRCA[MAXORDERNUM]; 
double      elementDiameters_meanRCA[MAXORDERNUM];

double segmentElementRatio_meanLAD[MAXORDERNUM]; 
double      segmentElementRatio_SDLAD[MAXORDERNUM]; 
double           elementLengths_meanLAD[MAXORDERNUM]; 
double                elementLengths_SDLAD[MAXORDERNUM]; 
double          segmentLengths_meanLAD[MAXORDERNUM]; 
double               segmentLengths_SDLAD[MAXORDERNUM]; 
double      elementDiameters_meanLAD[MAXORDERNUM];

double segmentElementRatio_meanLCX[MAXORDERNUM]; 
double      segmentElementRatio_SDLCX[MAXORDERNUM]; 
double           elementLengths_meanLCX[MAXORDERNUM]; 
double                elementLengths_SDLCX[MAXORDERNUM]; 
double          segmentLengths_meanLCX[MAXORDERNUM]; 
double               segmentLengths_SDLCX[MAXORDERNUM]; 
double      elementDiameters_meanLCX[MAXORDERNUM];

int numNodesInOrder[MAXORDERNUM];


int order, temp_order;
double **connRCA;
double **connLAD;
double **connLCX;
double **prob_connRCA, 	**prob_connLAD, 	**prob_connLCX;
// double cuprob_connRCA[12][12];
double **cuprob_connRCA;
double **cuprob_connLAD;
double **cuprob_connLCX;
double total_col1, total_col2, total_col3;
int colum, row, row2;
int whichCase;
unsigned int *nodeData, *usednodeData; // dynamically allocate this.

// seed the random number generator with system time.
/* Initialising the random number generator is important. Do this once. */
random = fopen("/dev/urandom","rb");
fread(&somerandomnumber, sizeof(unsigned int), 1, random);
fclose(random);
init_genrand(somerandomnumber);


// initialise PetSc and MPI.
//  PetscInitialize(&argc,&argv,(char*)0,help);
  if(debug)
    printf("Human bi-ventricle model geometry, half ellipsoids, 1 proc job. \n");
    if(debug)
  printf("Rectangle model size: %d %d\n", usr_MX, usr_MY);
/* Construct data structures. */
		/* Declare the geometry memory. This declares memory on each proc. */
/*
			user.geometry	 								= (PetscInt***)	calloc(usr_MY,sizeof(PetscInt **));
			user.vasculature		 							= (PetscReal** )		calloc(usr_MY,sizeof(PetscReal *)) ;
user.vasc_output		 							= (PetscReal** )		calloc(usr_MY,sizeof(PetscReal *)) ;			
		user.spheres    									= (PetscInt  **)  	calloc(usr_MY, sizeof(PetscInt * ));			
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.geometry[usr_j] 							= (PetscInt **) 	calloc(usr_MX,sizeof(PetscInt*))  ;
				user.vasculature[usr_j] 						= (PetscReal * ) 	calloc(usr_MX,sizeof(PetscReal))   ;								
				user.vasc_output[usr_j] 						= (PetscInt * ) 	calloc(usr_MX,sizeof(PetscInt))   ;				
				user.spheres[usr_j] 							= (PetscInt*)  		calloc(usr_MX, sizeof(PetscInt ));				
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.geometry[usr_j][usr_i] 				= (PetscInt*  ) 	calloc(NBS,sizeof(PetscInt))         ;				
				}
			}


									 		ThisNodeIsALeaf  			   = (int  *)     calloc(nuNodes, sizeof(int))	 ;
									 		setUp			  			   = (int  *)     calloc(nuNodes, sizeof(int))	 ;
									 		parent			  			   = (int  *)     calloc(nuNodes, sizeof(int))	 ;
									        
// srkNodes are the segments from me to Parent. Total number of segments is nuNodes = 2 * (number of leaves) - 1.
									 		srkNodes                 			   = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 srkNodes[temp_x] 				= (int*) 	calloc(2, sizeof(int))			;
									 		srkStrahlerNumberOfNode            = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 srkStrahlerNumberOfNode[temp_x] = (int*) 	calloc(2, sizeof(int))			;
									 		srkLevelOfNode                 		   = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 srkLevelOfNode[temp_x] 			= (int*) 	calloc(2, sizeof(int))			; // 0 is me, 1 is my parent.



									 		coodsMe                 			   = (double  **)     calloc(nuNodes, sizeof(double*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 coodsMe[temp_x] 				= (double*) 	calloc(3, sizeof(double))			; // 0 is x, 1 is y, 2 is z.
									 		NewcoodsMe                 		   = (double  **)     calloc(nuNodes, sizeof(double*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewcoodsMe[temp_x] 			= (double*) 	calloc(3, sizeof(double))			; // 0 is x, 1 is y, 2 is z.


									 		boundingBox        	= (double  **)  calloc(nuNodes, sizeof(double*));
for(temp_x = 0; temp_x < nuNodes; temp_x++) boundingBox[temp_x] = (double*) 	calloc(4, sizeof(double))		; // 0 is x, 1 is y, 2 is z.

									 	NewboundingBox        	= (double  **)  calloc(nuNodes, sizeof(double*));
for(temp_x = 0; temp_x < nuNodes; temp_x++) NewboundingBox[temp_x] = (double*) 	calloc(4, sizeof(double))		; // 0 is x, 1 is y, 2 is z.

									 	NewRboundingBox        	= (double  **)  calloc(nuNodes, sizeof(double*));
for(temp_x = 0; temp_x < nuNodes; temp_x++) NewRboundingBox[temp_x] = (double*) 	calloc(4, sizeof(double))		; // 0 is x, 1 is y, 2 is z.


									 		coodsParent                 			   = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 coodsParent[temp_x] 		       = (int*) 	calloc(3, sizeof(int))			; // 0 is x, 1 is y, 2 is z.
									 		NewcoodsParent                 	   = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewcoodsParent[temp_x] 		= (int*) 	calloc(3, sizeof(int))			; // 0 is x, 1 is y, 2 is z.

// this is set at initialisation time.
// fractions may not be needed.
									 		fractionNodes                 		  = (double  **)     calloc(nuNodes, sizeof(double*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 fractionNodes[temp_x] 			= (double*) 	calloc(2, sizeof(double))			; // 0 is left, 1 is right.
// this is calculated after.
									 		radiusSegment                 		   = (double  *)     calloc(nuNodes, sizeof(double))	 ;

									 		NewfractionNodes                 	   = (double  **)     calloc(nuNodes, sizeof(double*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewfractionNodes[temp_x] 		= (double*) 	calloc(2, sizeof(double))			; // 0 is left, 1 is right.
									 		NewradiusSegment                 	   = (double  *)     calloc(nuNodes, sizeof(double))	 ;

									 		resistivity					   = (double  *)     calloc(nuNodes, sizeof(double));
									 		perfusion					   = (double  *)     calloc(nuNodes, sizeof(double));
									 		pressure						   = (double  *)     calloc(nuNodes, sizeof(double));
											downResistance				   = (double  *)     calloc(nuNodes, sizeof(double));
											NewdownResistance		  	  = (double  *)     calloc(nuNodes, sizeof(double));
											NewRdownResistance		  = (double  *)     calloc(nuNodes, sizeof(double));

									 	Newresistivity					   = (double  *)     calloc(nuNodes, sizeof(double));
									 	Newperfusion					   = (double  *)     calloc(nuNodes, sizeof(double));
									 	Newpressure					   = (double  *)     calloc(nuNodes, sizeof(double));

									 	NewRresistivity					   = (double  *)     calloc(nuNodes, sizeof(double));
									 	NewRperfusion					   = (double  *)     calloc(nuNodes, sizeof(double));
									 	NewRpressure					   = (double  *)     calloc(nuNodes, sizeof(double));

									 	lengthSegment                 		   = (double  *)     calloc(nuNodes, sizeof(double));
										NewlengthSegment          		   = (double  *)     calloc(nuNodes, sizeof(double));
										NewRlengthSegment       		   = (double  *)     calloc(nuNodes, sizeof(double));

										leaves 			 = (int *) calloc(nuLeaves, sizeof(int));




			for (usr_j = 0; usr_j < usr_MY; usr_j++)
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.geometry[usr_j][usr_i][0] 				= 1;
				}


// data for when the sub-trees start moving around.
									 		NewRsrkNodes                 = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewRsrkNodes[temp_x] = (int*) 	calloc(2, sizeof(int))			;

									 		NewRsrkLevelOfNode                 = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewRsrkLevelOfNode[temp_x] = (int*) 	calloc(2, sizeof(int))			;

									 		NewRThisNodeIsALeaf  = (int  *)     calloc(nuNodes, sizeof(int))	 ;
									 		NewRsrkStrahlerNumberOfNode                 = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewRsrkStrahlerNumberOfNode[temp_x] = (int*) 	calloc(2, sizeof(int))			;

									 		NewRcoodsMe                 = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewRcoodsMe[temp_x] = (int*) 	calloc(3, sizeof(int))			; // 0 is x, 1 is y, 2 is z.
									 		NewRcoodsParent                 = (int  **)     calloc(nuNodes, sizeof(int*))	 ;
for(temp_x = 0; temp_x < nuNodes; temp_x++) 	 NewRcoodsParent[temp_x] = (int*) 	calloc(3, sizeof(int))			; // 0 is x, 1 is y, 2 is z.



									               NewRradiusSegment                 = (double  *)     calloc(nuNodes, sizeof(double))	 ;

in_inorder  	= (int *) calloc(nuNodes-1, sizeof(int)) ; in_preorder  	= (int *) calloc(nuNodes-1, sizeof(int)) ; in_postorder  = (int *) calloc(nuNodes-1, sizeof(int)) ;
*/

									 cuprob_connRCA 			= (double  **)  calloc(12, sizeof(double*));
									 cuprob_connLAD 			= (double  **)  calloc(12, sizeof(double*));
									 cuprob_connLCX 			= (double  **)  calloc(12, sizeof(double*));
									 
for(temp_x=0; temp_x<12; temp_x++) 		cuprob_connRCA[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;
for(temp_x=0; temp_x<12; temp_x++) 		cuprob_connLAD[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;
for(temp_x=0; temp_x<12; temp_x++) 		cuprob_connLCX[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;

									 prob_connRCA 			= (double  **)  calloc(12, sizeof(double*));
									 prob_connLAD 			= (double  **)  calloc(12, sizeof(double*));
									 prob_connLCX 			= (double  **)  calloc(12, sizeof(double*));
									 
for(temp_x=0; temp_x<12; temp_x++) 		prob_connRCA[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;
for(temp_x=0; temp_x<12; temp_x++) 		prob_connLAD[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;
for(temp_x=0; temp_x<12; temp_x++) 		prob_connLCX[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;

									 connRCA 			= (double  **)  calloc(12, sizeof(double*));
									 connLAD 			= (double  **)  calloc(12, sizeof(double*));
									 connLCX 			= (double  **)  calloc(12, sizeof(double*));
									 
for(temp_x=0; temp_x<12; temp_x++) 		connRCA[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;
for(temp_x=0; temp_x<12; temp_x++) 		connLAD[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;
for(temp_x=0; temp_x<12; temp_x++) 		connLCX[temp_x] 	  = (double*)     calloc(12, sizeof(double)) ;


nodeData 		= (unsigned int *) calloc(billion, sizeof(unsigned int)) ;
usednodeData 	= (unsigned int *) calloc(billion, sizeof(unsigned int)) ;

