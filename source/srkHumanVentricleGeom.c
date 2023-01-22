/*
Sanjay Kharche.
16 January 2017.

Idea is to generate a bi-ventricular geometry with smooth (C0, C1) fibre orientation.
Starting point is Kuhl/Sermesant papers, and their reference to Mercier 1982 paper.

The Z axis is the LV long axis. It is shared with the long axis of the RV.

1 Feb 2017.

Generate a Purkinje network after having defined the 3D ventricles, and their boundary conditions for fibres.
*/

static char help[] = "Synthetic 3D geometry for human ventricles, 1 process job.\n";

#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscsys.h>

/* my standard sundials headers, constants, and macros that will fit Petsc.
   Petsc cannot handle cvodes, at least not my installations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>       /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>            /* prototype for CVDense                */
#include <sundials/sundials_dense.h>      /* definitions DlsMat DENSE_ELEM        */
#include <sundials/sundials_types.h>      /* definition of type realtype          */
#include <sundials/sundials_math.h>

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */
#define RTOL        RCONST(1.0e-12)   	  /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-6)        /* scalar absolute tolerance components */
#define MAXSTEPS    5000
#define ZERO        RCONST(0.0)

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */

/*
This works at 0.1 mm, 0.2 mm, and probably works at 0.4. As far as possible, keep DX an integral multiple of 0.1 mm.
*/
#define DX               0.4                 /* X internode spacing                       */
// #define DY              0.2                 /* Y internode spacing                       */
// #define DZ	          0.2 		  		/* Z internode spacing                       */


/*
Dimensions of the model. These dimenstions are taken from Comput Mech (2010) 45: 227-243, page 240, Figure 7.
A margin of 5 on each side of the model.
*/
#define LV_X_RADIUS 30
#define LV_Y_RADIUS 30
#define LV_Z_RADIUS 70
#define LV_THICKNESS 12

#define RV_X_RADIUS 30
#define RV_Y_RADIUS 52
#define RV_Z_RADIUS 60
#define RV_THICKNESS 6

#define EPI_PHI        1.0
#define ENDO_PHI  -1.0

#define MARGIN 5

/*
EPI is the outer surface of the ventricles, as well as the suface of the LV facing the RV in the septum. Septum is not segmented as of now.
ENDO is the inner surface of the ventricles. The Septum has inner surface on the LV side and outer surface on the RV side.
*/
#define EPI      1000
#define ENDO 2000

/*
This needs working out for each resolution you choose.
5 nodes on each side for margins.
The LV is 70 mm tall. That makes 350 nodes in Z direction.
The LV radius is 30 mm in x-y plane, diameter is 60 mm. That makes it 300 nodes.
*/
// #define usr_MX 310    // max in x direction: 0 to usr_MX
// #define usr_MY 420    // max in x direction: 0 to usr_MY
// #define usr_MZ 360     // max in x direction: 0 to usr_MZ. The cylinder is along the Z axis. This is the size at production runs
#define NBS    7            // 3D arrays have 7 neighbours: itself, and 6 surrounding it in a standard 1st order FD stencil.

/* User-defined data structures and routines           */
/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da; // DM instance. This is required.
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
/* all skPetsc programs use geomtry as: 0 = itself, 1 = y+1, 2 = x+1, 3 = y -1, 4 = x - 1, 5 = z + 1, 6 = z - 1 */
  PetscInt     ****geometry    ; // 3D models have 4D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1.
  PetscInt      ***epiendo       ;

// fibres and sheets arrays.
PetscReal ****n;
PetscReal ****c;
PetscReal ****ncz;
PetscReal ****p;
PetscReal ****sheets;
PetscReal ****fibres;
 
} AppCtx;

// arguments are: the purkinje 4D array, starting and end point x,y,z's, and the array size, radius of the pf cylinder.
void cylinder(int ****pf, int xx1, int yy1, int zz1, int xx2, int yy2, int zz2, int usr_MX, int usr_MY, int usr_MZ, int pfrad){

	int usr_i, usr_j, usr_k;
	double t;
	int x0, y0, z0;
	int diff;

	for(t = 0; t <= 1.0; t = t + 0.005){

	x0 = xx1 + (int)((double)(xx2 - xx1) * t);
	y0 = yy1 + (int)((double)(yy2 - yy1) * t);	
	z0 = zz1 + (int)((double)(zz2 - zz1) * t);	
	for(usr_k = 1; usr_k <= usr_MZ-1; usr_k++)
		for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
			for(usr_i = 1; usr_i < usr_MX-1; usr_i++){
				diff = (usr_i - x0) * (usr_i - x0) + (usr_j - y0) * (usr_j - y0) + (usr_k - z0) * (usr_k - z0) - pfrad * pfrad;
				if( diff < 0 ) pf[usr_k][usr_j][usr_i][0] = 1;	
	}

	} // end of t loop.

} // end of void cylinder.

int main(int argc,char **argv)
{
/* Program reduced to use of FD with colouring. */
  Vec            u		    			; /* solution, residual vectors. You need the r for SNES */
  AppCtx         user		    	; /* user-defined work context                           */
  DM             da		    		;
  PetscInt       tissue_type		;

	// sk variables.
	  PetscInt       mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z, usr_i, usr_j, usr_k;
	  PetscScalar    ***u_localptr;
	  char 		str[1000]; // so you do not keep creating and destroying this.
          PetscInt       temp_x, temp_y, temp_z;
         FILE               *output;

         int ****purkinje; // the Purkinje fibre array. 
	int PF_radius;	// radius of cylinders representing Purkinje fibres in number of FD nodes.
         
         int LV_x_c, LV_y_c, LV_z_c;
         int LV_x_radius,   LV_y_radius,   LV_z_radius;
         int LV_x_radius2, LV_y_radius2, LV_z_radius2;

        int RV_x_c, RV_y_c, RV_z_c;
         int RV_x_radius,   RV_y_radius,   RV_z_radius;
         int RV_x_radius2, RV_y_radius2, RV_z_radius2;
         
         double tempx1, tempy1, tempz1, tempx2, tempy2, tempz2;
         PetscInt usr_MX, usr_MY, usr_MZ, tempRADIUS;
         
         double alpha;
         double A_star, B_star, C_star, N_amplitude;
         
// initialise PetSc and MPI.
  PetscInitialize(&argc,&argv,(char*)0,help);
/*
#define LV_X_RADIUS 30
#define LV_Y_RADIUS 30
#define LV_Z_RADIUS 70

#define RV_X_RADIUS 30
#define RV_Y_RADIUS 52
#define RV_Z_RADIUS 60
*/
  // you need the usr_MX/MY/MZ before creating the DMDA objects.
// the 2 ventricles are in parallel along the X axis.
if(LV_X_RADIUS >= RV_X_RADIUS) tempRADIUS = LV_X_RADIUS; else tempRADIUS = RV_X_RADIUS;
  usr_MX = (PetscInt) ( (double) ( ( 2* tempRADIUS ) ) / (double)(DX) )  + 2 * MARGIN;
  // the 2 ventricles are in series along the X axis. This is not necessarily a radius that needs doubling, see the Figure 7 in Comput Mech 2010 45: 227-243.
 usr_MY = (PetscInt) ( (double) ( ( LV_X_RADIUS + RV_Y_RADIUS ) ) / (double)(DX) ) + 2 * MARGIN;  
 // the 2 ventricles are in parallel along the X axis.
if(LV_Z_RADIUS >= RV_Z_RADIUS) tempRADIUS = LV_Z_RADIUS; else tempRADIUS = RV_Z_RADIUS;
  usr_MZ = (PetscInt) ( (double) ( ( tempRADIUS ) ) / (double)(DX) ) + 2 * MARGIN;

  printf("%d %d %d\n", usr_MX, usr_MY, usr_MZ);
	     /* Initialize user application context */
	      user.da           = NULL;
	     // the -11, -11 becomes usr_MX, usr_MY when I have #define usr_MX and #define usr_MY (potentially #define usr_MZ) put those values in place of -11,-11
   DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,usr_MX,usr_MY,usr_MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,NULL,&da);  /* this decides that u is non-ghosted 1st order star stencil vector. */
	     user.da = da		;
	     DMCreateGlobalVector(da,&u); 
	   DMDAGetCorners(da,&mybase_x,&mybase_y,&mybase_z,&mysize_x,&mysize_y,&mysize_z);  // you need this for running the Sundials loops.

  printf("Human bi-ventricle model geometry, half ellipsoids, 1 proc job. \n");
		/* Declare the geometry memory. This declares memory on each proc. */
		user.geometry = (PetscInt  ****) 	calloc(usr_MZ, sizeof(PetscInt ***));
		user.epiendo    = (PetscInt  ***)  	calloc(usr_MZ, sizeof(PetscInt ** ));
		purkinje    	   = (int  		****)  calloc(usr_MZ, sizeof(int *** ));				
		for(usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.geometry[usr_k] = (PetscInt***)  	calloc(usr_MY,sizeof(PetscInt **));
			user.epiendo[usr_k] = (PetscInt**)  		calloc(usr_MY,sizeof(PetscInt *));
			purkinje[usr_k]    = (int  ***)  			calloc(usr_MY, sizeof(int ** ));							
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.geometry[usr_k][usr_j] = (PetscInt **) 	calloc(usr_MX,sizeof(PetscInt*));
				user.epiendo[usr_k][usr_j] = (PetscInt *) 		calloc(usr_MX,sizeof(PetscInt));				
				purkinje[usr_k][usr_j]    = (int  **)  			calloc(usr_MX, sizeof(int * ));				
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.geometry[usr_k][usr_j][usr_i] = (PetscInt*) 	calloc(NBS,sizeof(PetscInt));				
				purkinje[usr_k][usr_j][usr_i]    = (int  *)  				calloc(NBS, sizeof(int  ));									
				}
			}
		}

		user.sheets     	= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));		
		user.fibres 		= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));
		user.n		 	= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));				
		user.ncz		 	= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));						
		user.p		 	= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));						
		user.c		 	= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));				
		
		for(usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.sheets[usr_k] 	= (PetscReal***)  calloc(usr_MY,sizeof(PetscReal **));
			user.fibres[usr_k] 	= (PetscReal***)  calloc(usr_MY,sizeof(PetscReal **));

		user.n[usr_k]		 	= (PetscReal  ***) calloc(usr_MY, sizeof(PetscReal **));				
		user.ncz[usr_k]	 	= (PetscReal  ***) calloc(usr_MY, sizeof(PetscReal **));						
		user.p[usr_k]		 	= (PetscReal  ***) calloc(usr_MY, sizeof(PetscReal **));						
		user.c[usr_k]		 	= (PetscReal  ***) calloc(usr_MY, sizeof(PetscReal **));				
			
			
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.sheets[usr_k][usr_j] 		= (PetscReal**)  calloc(usr_MX,sizeof(PetscReal *));
				user.fibres[usr_k][usr_j] 	= (PetscReal**)  calloc(usr_MX,sizeof(PetscReal *));

		user.n[usr_k][usr_j]	 	= (PetscReal  **) calloc(usr_MX, sizeof(PetscReal * ));				
		user.ncz[usr_k][usr_j]	 	= (PetscReal  **) calloc(usr_MX, sizeof(PetscReal * ));						
		user.p[usr_k][usr_j]	 	= (PetscReal  **) calloc(usr_MX, sizeof(PetscReal * ));						
		user.c[usr_k][usr_j]	 	= (PetscReal  **) calloc(usr_MX, sizeof(PetscReal * ));				
				
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.sheets[usr_k][usr_j][usr_i] = (PetscReal*) calloc(3, sizeof(PetscReal));
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal*) calloc(3, sizeof(PetscReal));									
		user.n[usr_k][usr_j][usr_i]	 	= (PetscReal  *) calloc(3, sizeof(PetscReal ));				
		user.ncz[usr_k][usr_j][usr_i] 	= (PetscReal  *) calloc(3, sizeof(PetscReal ));						
		user.p[usr_k][usr_j][usr_i]	 	= (PetscReal  *) calloc(3, sizeof(PetscReal ));						
		user.c[usr_k][usr_j][usr_i]	 	= (PetscReal  *) calloc(3, sizeof(PetscReal ));				
				}
			}
		}

/*
The LV is at min X and min Y of the model. The RV is further out.
#define LV_X_RADIUS 30
#define LV_Y_RADIUS 30
#define LV_Z_RADIUS 70

#define RV_X_RADIUS 30
#define RV_Y_RADIUS 52
#define RV_Z_RADIUS 60
*/
	
// centre of the LV ellipsoids.				
LV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX) + MARGIN; 
LV_y_c 	       	= (PetscInt)((double)LV_Y_RADIUS / DX) + MARGIN; 
LV_z_c 	   	= MARGIN;
LV_x_radius	= (PetscInt)((double)LV_X_RADIUS / DX);      // outer ellipsoid
LV_y_radius 	= (PetscInt)((double)LV_Y_RADIUS / DX);      
LV_z_radius 	= (PetscInt)((double)LV_Z_RADIUS / DX);
LV_x_radius2 = (PetscInt)((double)LV_X_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX);    // inner ellipsoid
LV_y_radius2 = (PetscInt)((double)LV_Y_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_z_radius2 = (PetscInt)((double)LV_Z_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX);
printf("LV: %d %d %d %d %d %d %d %d %d\n",LV_x_c, LV_y_c, LV_z_c, LV_x_radius, LV_y_radius, LV_z_radius, LV_x_radius2, LV_y_radius2, LV_z_radius2);
		for(usr_k = MARGIN; usr_k < usr_MZ; usr_k++)
			for (usr_j = 0; usr_j < usr_MY; usr_j++)
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
// epicardial surface
tempx1 = (double) ( ( usr_i - LV_x_c) * ( usr_i - LV_x_c) )/    (double) (LV_x_radius * LV_x_radius);
tempy1 = (double) ( ( usr_j - LV_y_c) * ( usr_j - LV_y_c) )/  (double) (LV_y_radius * LV_y_radius);
tempz1 = (double) ( ( usr_k - LV_z_c) * ( usr_k - LV_z_c) )/ (double) (LV_z_radius * LV_z_radius );
// endocardial surface
tempx2 = (double) ( ( usr_i - LV_x_c) * ( usr_i - LV_x_c) )/    (double) (LV_x_radius2 * LV_x_radius2);
tempy2 = (double) ( ( usr_j - LV_y_c) * ( usr_j - LV_y_c) )/  (double) (LV_y_radius2 * LV_y_radius2);
tempz2 = (double) ( ( usr_k - LV_z_c) * ( usr_k - LV_z_c) )/ (double) (LV_z_radius2 * LV_z_radius2 );

if( (tempx1 + tempy1 + tempz1 <= 1) && (tempx2 + tempy2 + tempz2 >=1)  ) user.geometry[usr_k][usr_j][usr_i][0] = 1; else user.geometry[usr_k][usr_j][usr_i][0] = 0;
				}

// now the RV
/*
RV_x_c = 150;  RV_y_c = 205; RV_z_c = 5;
RV_x_radius   = 150;         RV_y_radius   = 260;         RV_z_radius 	= 300;
RV_x_radius2 = 150 - 30; RV_y_radius2 = 260 - 30; RV_z_radius2 	= 300 - 30;

#define LV_X_RADIUS 30
#define LV_Y_RADIUS 30
#define LV_Z_RADIUS 70
#define LV_THICKNESS 12

#define RV_X_RADIUS 30
#define RV_Y_RADIUS 52
#define RV_Z_RADIUS 60
#define RV_THICKNESS 6

#define MARGIN 5

*/
RV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX)+MARGIN; 
RV_y_c 	       	= (PetscInt)((double)(LV_Y_RADIUS) / DX) + MARGIN;  // the X and Y coordinates of this ellipsoid are the same as for LV.
RV_z_c 	   	= MARGIN;
RV_x_radius	= (PetscInt)((double)RV_X_RADIUS / DX);      
RV_y_radius = (PetscInt)((double)RV_Y_RADIUS / DX);      
RV_z_radius = (PetscInt)((double)RV_Z_RADIUS / DX);
RV_x_radius2 = (PetscInt)((double)RV_X_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX); 
RV_y_radius2 = (PetscInt)((double)RV_Y_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX); 
RV_z_radius2 = (PetscInt)((double)RV_Z_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX);

printf("RV: %d %d %d %d %d %d %d %d %d\n",RV_x_c, RV_y_c, RV_z_c, RV_x_radius, RV_y_radius, RV_z_radius, RV_x_radius2, RV_y_radius2, RV_z_radius2);
		for(usr_k = MARGIN; usr_k < usr_MZ; usr_k++)
			for (usr_j = RV_y_c; usr_j < usr_MY; usr_j++)
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
				// epicardial surface
				tempx1 = (double) ( ( usr_i - RV_x_c) * ( usr_i - RV_x_c) )/    (double) (RV_x_radius * RV_x_radius);
				tempy1 = (double) ( ( usr_j - RV_y_c) * ( usr_j - RV_y_c) )/  (double) (RV_y_radius * RV_y_radius);
				tempz1 = (double) ( ( usr_k - RV_z_c) * ( usr_k - RV_z_c) )/ (double) (RV_z_radius * RV_z_radius );
				// endocardial surface
				tempx2 = (double) ( ( usr_i - RV_x_c) * ( usr_i - RV_x_c) )/    (double) (RV_x_radius2 * RV_x_radius2);
				tempy2 = (double) ( ( usr_j - RV_y_c) * ( usr_j - RV_y_c) )/  (double) (RV_y_radius2 * RV_y_radius2);
				tempz2 = (double) ( ( usr_k - RV_z_c) * ( usr_k - RV_z_c) )/ (double) (RV_z_radius2 * RV_z_radius2 );
				if( (tempx1 + tempy1 + tempz1 <= 1) && (tempx2 + tempy2 + tempz2 >=1) && user.geometry[usr_k][usr_j][usr_i][0] == 0 ){
					user.geometry[usr_k][usr_j][usr_i][0] = 2;
	
				// LV check.
				// epicardial surface
				tempx1 = (double) ( ( usr_i - LV_x_c) * ( usr_i - LV_x_c) )/    (double) (LV_x_radius * LV_x_radius);
				tempy1 = (double) ( ( usr_j - LV_y_c) * ( usr_j - LV_y_c) )/  (double) (LV_y_radius * LV_y_radius);
				tempz1 = (double) ( ( usr_k - LV_z_c) * ( usr_k - LV_z_c) )/ (double) (LV_z_radius * LV_z_radius );
				// endocardial surface
				tempx2 = (double) ( ( usr_i - LV_x_c) * ( usr_i - LV_x_c) )/    (double) (LV_x_radius2 * LV_x_radius2);
				tempy2 = (double) ( ( usr_j - LV_y_c) * ( usr_j - LV_y_c) )/  (double) (LV_y_radius2 * LV_y_radius2);
				tempz2 = (double) ( ( usr_k - LV_z_c) * ( usr_k - LV_z_c) )/ (double) (LV_z_radius2 * LV_z_radius2 );
				if( (tempx2 + tempy2 + tempz2 < 1) && user.geometry[usr_k][usr_j][usr_i][0] == 2) user.geometry[usr_k][usr_j][usr_i][0] = 0;	
				} // end of if tempx1 etc.
			} // end of outer loops.


// septum is intersection of RV and LV inner spaces. Septum as of now is all tissue_type = 1. Do it after.


// neighbours.
// do neighbours here.
// get the neighbours in my numbering scheme.
for(usr_k = 0; usr_k < usr_MZ; usr_k++)
for(usr_j = 0; usr_j < usr_MY; usr_j++)
for(usr_i = 0; usr_i < usr_MX; usr_i++){
// printf("got here %d %d %d\n",usr_i,usr_j,usr_k);
		// top y + 1
		temp_y = usr_j + 1;
		if(temp_y<usr_MY) user.geometry[usr_k][usr_j][usr_i][1] = user.geometry[usr_k][temp_y][usr_i][0];
		else 		  user.geometry[usr_k][usr_j][usr_i][1] = -1;
		// right
		temp_x = usr_i + 1;
		if(temp_x<usr_MX) user.geometry[usr_k][usr_j][usr_i][2] = user.geometry[usr_k][usr_j][temp_x][0];
		else 	          user.geometry[usr_k][usr_j][usr_i][2] = -1;
		// bottom
		temp_y = usr_j - 1;
		if(temp_y>=0) 	    user.geometry[usr_k][usr_j][usr_i][3] = user.geometry[usr_k][temp_y][usr_i][0];
		else 		    user.geometry[usr_k][usr_j][usr_i][3] = -1;
		// left
		temp_x = usr_i - 1;
		if(temp_x>=0) 	    user.geometry[usr_k][usr_j][usr_i][4] = user.geometry[usr_k][usr_j][temp_x][0];
		else 		    user.geometry[usr_k][usr_j][usr_i][4] = -1;
		// +z
		temp_z = usr_k + 1;
		if(temp_z<usr_MZ) user.geometry[usr_k][usr_j][usr_i][5] = user.geometry[temp_z][usr_j][usr_i][0];
		else 		  user.geometry[usr_k][usr_j][usr_i][5] = -1;
		// -z
		temp_z = usr_k - 1;
		if(temp_z>=0) 	    user.geometry[usr_k][usr_j][usr_i][6] = user.geometry[temp_z][usr_j][usr_i][0];
		else 		    user.geometry[usr_k][usr_j][usr_i][6] = -1;
}

// epi and endo array.
// the geometry is defined, using that array.
// none of the FD nodes are close to a hard boundary.
// my previous definition was dense, but not a single star stencil thick. I need it 1 cell/1 star stencil thick so that I can define the Z vector correctly.
for(usr_k = 1; usr_k < usr_MZ-1; usr_k++)
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++)
if( user.geometry[usr_k][usr_j][usr_i][0] == 1 && (user.geometry[usr_k][usr_j][usr_i][1] == 0 || user.geometry[usr_k][usr_j][usr_i][2] == 0 || user.geometry[usr_k][usr_j][usr_i][3] == 0 || user.geometry[usr_k][usr_j][usr_i][4] == 0 || user.geometry[usr_k][usr_j][usr_i][5] == 0 || user.geometry[usr_k][usr_j][usr_i][6] == 0) )  {

LV_x_radius2 = (PetscInt)((double)LV_X_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / ( 2.0 * DX) );    // inner ellipsoid
LV_y_radius2 = (PetscInt)((double)LV_Y_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / ( 2.0 * DX) ); 
LV_z_radius2 = (PetscInt)((double)LV_Z_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / ( 2.0 * DX ));	
	
	// epicardial surface
tempx1 = (double) ( ( usr_i - LV_x_c) * ( usr_i - LV_x_c) )/    (double) (LV_x_radius * LV_x_radius);
tempy1 = (double) ( ( usr_j - LV_y_c) * ( usr_j - LV_y_c) )/  (double) (LV_y_radius * LV_y_radius);
tempz1 = (double) ( ( usr_k - LV_z_c) * ( usr_k - LV_z_c) )/ (double) (LV_z_radius * LV_z_radius );
// endocardial surface
tempx2 = (double) ( ( usr_i - LV_x_c) * ( usr_i - LV_x_c) )/    (double) (LV_x_radius2 * LV_x_radius2);
tempy2 = (double) ( ( usr_j - LV_y_c) * ( usr_j - LV_y_c) )/  (double) (LV_y_radius2 * LV_y_radius2);
tempz2 = (double) ( ( usr_k - LV_z_c) * ( usr_k - LV_z_c) )/ (double) (LV_z_radius2 * LV_z_radius2 );

user.epiendo[usr_k][usr_j][usr_i] = 3000;

// epi LV surface
if(tempx2 + tempy2 + tempz2 > 1.0 ) user.epiendo[usr_k][usr_j][usr_i] = EPI; 

// endo LV surface
if(tempx2 + tempy2+ tempz2 < 1.0 ) user.epiendo[usr_k][usr_j][usr_i] = ENDO; 

if(usr_k==MARGIN) user.epiendo[usr_k][usr_j][usr_i] = EPI;   // this needs to be improved.

}

// RV.
for(usr_k = 1; usr_k < usr_MZ-1; usr_k++)
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++)
if( user.geometry[usr_k][usr_j][usr_i][0] == 2 && (user.geometry[usr_k][usr_j][usr_i][1] == 0 || user.geometry[usr_k][usr_j][usr_i][2] == 0 || user.geometry[usr_k][usr_j][usr_i][3] == 0 || user.geometry[usr_k][usr_j][usr_i][4] == 0 || user.geometry[usr_k][usr_j][usr_i][5] == 0 || user.geometry[usr_k][usr_j][usr_i][6] == 0) )  {

RV_x_radius2 = (PetscInt)((double)RV_X_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS /  (2.0 * DX ) ); 
RV_y_radius2 = (PetscInt)((double)RV_Y_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / (2.0 * DX) ); 
RV_z_radius2 = (PetscInt)((double)RV_Z_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / (2.0 * DX) );
	
	
	// epicardial surface
tempx1 = (double) ( ( usr_i - RV_x_c) * ( usr_i - RV_x_c) )/    (double) (RV_x_radius * RV_x_radius);
tempy1 = (double) ( ( usr_j - RV_y_c) * ( usr_j - RV_y_c) )/  (double) (RV_y_radius * RV_y_radius);
tempz1 = (double) ( ( usr_k - RV_z_c) * ( usr_k - RV_z_c) )/ (double) (RV_z_radius * RV_z_radius );
// endocardial surface
tempx2 = (double) ( ( usr_i - RV_x_c) * ( usr_i - RV_x_c) )/    (double) (RV_x_radius2 * RV_x_radius2);
tempy2 = (double) ( ( usr_j - RV_y_c) * ( usr_j - RV_y_c) )/  (double) (RV_y_radius2 * RV_y_radius2);
tempz2 = (double) ( ( usr_k - RV_z_c) * ( usr_k - RV_z_c) )/ (double) (RV_z_radius2 * RV_z_radius2 );

user.epiendo[usr_k][usr_j][usr_i] = 4000;

// epi RV surface
if(tempx2 + tempy2 + tempz2 > 1.0 ) user.epiendo[usr_k][usr_j][usr_i] = EPI; 

// endo RV surface
if(tempx2 + tempy2+ tempz2 < 1.0 ) user.epiendo[usr_k][usr_j][usr_i] = ENDO; 

if(usr_k==MARGIN) user.epiendo[usr_k][usr_j][usr_i] = EPI; 

}


/*
Info.

LV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX)+MARGIN; 
LV_y_c 	       	= (PetscInt)((double)LV_Y_RADIUS / DX) + MARGIN; 
LV_z_c 	   	= MARGIN;
LV_x_radius	= (PetscInt)((double)LV_X_RADIUS / DX);      
LV_y_radius 	= (PetscInt)((double)LV_Y_RADIUS / DX);      
LV_z_radius 	= (PetscInt)((double)LV_Z_RADIUS / DX);
LV_x_radius2 = (PetscInt)((double)LV_X_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_y_radius2 = (PetscInt)((double)LV_Y_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_z_radius2 = (PetscInt)((double)LV_Z_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX);
*/

// see the local coordinate system first, then the sheets and fibres vectors.
// use the fibres and sheets arrays for now.
// node normal to ellipsoidal surface.
for(usr_k = 1; usr_k < usr_MZ-1; usr_k++)
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++){

LV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX)+MARGIN; 
LV_y_c 	       	= (PetscInt)((double)LV_Y_RADIUS / DX) + MARGIN; 
LV_z_c 	   	= MARGIN;
LV_x_radius	= (PetscInt)((double)LV_X_RADIUS / DX);      
LV_y_radius 	= (PetscInt)((double)LV_Y_RADIUS / DX);      
LV_z_radius 	= (PetscInt)((double)LV_Z_RADIUS / DX);
LV_x_radius2 = (PetscInt)((double)LV_X_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_y_radius2 = (PetscInt)((double)LV_Y_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_z_radius2 = (PetscInt)((double)LV_Z_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX);

RV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX)+MARGIN; 
RV_y_c 	       	= (PetscInt)((double)(LV_Y_RADIUS) / DX) + MARGIN;  // the X and Y coordinates of this ellipsoid are the same as for LV.
RV_z_c 	   	= MARGIN;
RV_x_radius	= (PetscInt)((double)RV_X_RADIUS / DX);      
RV_y_radius = (PetscInt)((double)RV_Y_RADIUS / DX);      
RV_z_radius = (PetscInt)((double)RV_Z_RADIUS / DX);
RV_x_radius2 = (PetscInt)((double)RV_X_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX); 
RV_y_radius2 = (PetscInt)((double)RV_Y_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX); 
RV_z_radius2 = (PetscInt)((double)RV_Z_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX);


	// LV epi.
	if(user.epiendo[usr_k][usr_j][usr_i] == EPI && user.geometry[usr_k][usr_j][usr_i][0] == 1){
		A_star = 1.0/ (double) (LV_x_radius * LV_x_radius);
		B_star = 1.0/ (double) (LV_x_radius * LV_y_radius);
		C_star = 1.0/ (double) (LV_x_radius * LV_z_radius);

		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - LV_x_c)  * (usr_i - LV_x_c)  + B_star * B_star * (usr_j - LV_y_c)  * (usr_j - LV_y_c)  + C_star * C_star * (usr_k - LV_z_c)  * (usr_k - LV_z_c) ) ) ; // this is half of what is written in the notes.
		user.n[usr_k][usr_j][usr_i][0] = A_star / N_amplitude * (usr_i - LV_x_c); 
		user.n[usr_k][usr_j][usr_i][1] =  B_star / N_amplitude * (usr_j - LV_y_c);
		user.n[usr_k][usr_j][usr_i][2] = C_star / N_amplitude * (usr_k - LV_z_c);
	}
// LV endo
	if(user.epiendo[usr_k][usr_j][usr_i] == ENDO && user.geometry[usr_k][usr_j][usr_i][0] == 1){
		A_star = 1.0/ (double) (LV_x_radius2 * LV_x_radius);
		B_star = 1.0/ (double) (LV_x_radius2 * LV_y_radius);
		C_star = 1.0/ (double) (LV_x_radius2 * LV_z_radius);
		
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - LV_x_c)  * (usr_i - LV_x_c)  + B_star * B_star * (usr_j - LV_y_c)  * (usr_j - LV_y_c)  + C_star * C_star * (usr_k - LV_z_c)  * (usr_k - LV_z_c) ) ) ; // this is half of what is written in the notes.
		user.n[usr_k][usr_j][usr_i][0] = A_star / N_amplitude * (usr_i - LV_x_c); 
		user.n[usr_k][usr_j][usr_i][1] =  B_star / N_amplitude * (usr_j - LV_y_c);
		user.n[usr_k][usr_j][usr_i][2] = C_star / N_amplitude * (usr_k - LV_z_c);
	}
// RV epi
	if(user.epiendo[usr_k][usr_j][usr_i] == EPI && user.geometry[usr_k][usr_j][usr_i][0] == 2){
		A_star = 1.0/ (double) (RV_x_radius * RV_x_radius);
		B_star = 1.0/ (double) (RV_x_radius * RV_y_radius);
		C_star = 1.0/ (double) (RV_x_radius * RV_z_radius);
		
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - RV_x_c)  * (usr_i - RV_x_c)  + B_star * B_star * (usr_j - RV_y_c)  * (usr_j - RV_y_c)  + C_star * C_star * (usr_k - RV_z_c)  * (usr_k - RV_z_c) ) ) ; // this is half of what is written in the notes.
		user.n[usr_k][usr_j][usr_i][0] = A_star / N_amplitude * (usr_i - RV_x_c); 
		user.n[usr_k][usr_j][usr_i][1] =  B_star / N_amplitude * (usr_j - RV_y_c);
		user.n[usr_k][usr_j][usr_i][2] = C_star / N_amplitude * (usr_k - RV_z_c);
	}
// RV endo
	if(user.epiendo[usr_k][usr_j][usr_i] == ENDO && user.geometry[usr_k][usr_j][usr_i][0] == 2){
		A_star = 1.0/ (double) (RV_x_radius2 * RV_x_radius2);
		B_star = 1.0/ (double) (RV_x_radius2 * RV_y_radius2);
		C_star = 1.0/ (double) (RV_x_radius2 * RV_z_radius2);
		
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - RV_x_c)  * (usr_i - RV_x_c)  + B_star * B_star * (usr_j - RV_y_c)  * (usr_j - RV_y_c)  + C_star * C_star * (usr_k - RV_z_c)  * (usr_k - RV_z_c) ) ) ; // this is half of what is written in the notes.
		user.n[usr_k][usr_j][usr_i][0] = A_star / N_amplitude * (double)(usr_i - RV_x_c); 
		user.n[usr_k][usr_j][usr_i][1] =  B_star / N_amplitude * (double)(usr_j - RV_y_c);
		user.n[usr_k][usr_j][usr_i][2] = C_star / N_amplitude * (double)(usr_k - RV_z_c);
	}
	
	if(user.epiendo[usr_k][usr_j][usr_i] > 0 && user.geometry[usr_k][usr_j][usr_i][0] > 0 && usr_k==MARGIN ){	
		user.n[usr_k][usr_j][usr_i][0] = 0.0; 
		user.n[usr_k][usr_j][usr_i][1] =  0.0;
		user.n[usr_k][usr_j][usr_i][2] = 1.0;		
	}
	
	
} // end of for loops.

// write
		sprintf(str,"n.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_fibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				fprintf(output,"%10.10F %10.10F %10.10F\n", user.n[usr_k][usr_j][usr_i][0], user.n[usr_k][usr_j][usr_i][1], user.n[usr_k][usr_j][usr_i][2] );
			}
		fclose(output);


for(usr_k = 0; usr_k < usr_MZ; usr_k++)
for(usr_j = 0; usr_j < usr_MY; usr_j++)
for(usr_i = 0; usr_i < usr_MX; usr_i++){
	if(user.epiendo[usr_k][usr_j][usr_i] == EPI){
		user.sheets[usr_k][usr_j][usr_i][0] =	user.n[usr_k][usr_j][usr_i][0];
		user.sheets[usr_k][usr_j][usr_i][1]  = 	user.n[usr_k][usr_j][usr_i][1];
		user.sheets[usr_k][usr_j][usr_i][2] =	user.n[usr_k][usr_j][usr_i][2];
	}

	if(user.epiendo[usr_k][usr_j][usr_i] == ENDO){
		user.sheets[usr_k][usr_j][usr_i][0] =	-user.n[usr_k][usr_j][usr_i][0];
		user.sheets[usr_k][usr_j][usr_i][1]  = 	-user.n[usr_k][usr_j][usr_i][1];
		user.sheets[usr_k][usr_j][usr_i][2] =	-user.n[usr_k][usr_j][usr_i][2];
	}
	
	
}

// write
		sprintf(str,"sheets.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_fibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				fprintf(output,"%10.10F %10.10F %10.10F\n", user.sheets[usr_k][usr_j][usr_i][0], user.sheets[usr_k][usr_j][usr_i][1], user.sheets[usr_k][usr_j][usr_i][2] );
			}
		fclose(output);


// take the Z vector along long axis. My long axis is along global Z axis.
// this is the c vector. ncz vector and c vectors both have the same amplitude. so do them together.
for(usr_k = 1; usr_k <= usr_MZ-1; usr_k++)
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++){
	
	if(user.epiendo[usr_k][usr_j][usr_i] == EPI && user.geometry[usr_k][usr_j][usr_i][0] == 1){ // LV epi
		A_star = 1.0/ (double) (LV_x_radius * LV_x_radius);
		B_star = 1.0/ (double) (LV_x_radius * LV_y_radius);
		C_star = 1.0/ (double) (LV_x_radius * LV_z_radius);
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - LV_x_c)  * (usr_i - LV_x_c)  + B_star * B_star * (usr_j - LV_y_c)  * (usr_j - LV_y_c) ) ) ; // this is half of what is written in the notes.
if(N_amplitude > 0.0){
		user.c[usr_k][usr_j][usr_i][0] = - B_star / N_amplitude * (usr_j - LV_y_c); 
		user.c[usr_k][usr_j][usr_i][1] =  A_star / N_amplitude * (usr_i - LV_x_c); 
		user.c[usr_k][usr_j][usr_i][2] = 0.0;

		user.ncz[usr_k][usr_j][usr_i][0] = A_star / N_amplitude * (usr_i - LV_x_c); 
		user.ncz[usr_k][usr_j][usr_i][1] =  B_star / N_amplitude * (usr_j - LV_y_c); 
		user.ncz[usr_k][usr_j][usr_i][2] = 0.0;
}
	}

		if(user.epiendo[usr_k][usr_j][usr_i] == ENDO && user.geometry[usr_k][usr_j][usr_i][0] == 1){ // LV endo
		A_star = 1.0/ (double) (LV_x_radius2 * LV_x_radius2);
		B_star = 1.0/ (double) (LV_x_radius2 * LV_y_radius2);
		C_star = 1.0/ (double) (LV_x_radius2 * LV_z_radius2);
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - LV_x_c)  * (usr_i - LV_x_c)  + B_star * B_star * (usr_j - LV_y_c)  * (usr_j - LV_y_c) ) ) ; // this is half of what is written in the notes.
if(N_amplitude > 0.0){
		user.c[usr_k][usr_j][usr_i][0] =  B_star / N_amplitude * (usr_j - LV_y_c);  // multiplied by the sign.
		user.c[usr_k][usr_j][usr_i][1] =  -A_star / N_amplitude * (usr_i - LV_x_c); 
		user.c[usr_k][usr_j][usr_i][2] = 0.0;
		
		user.ncz[usr_k][usr_j][usr_i][0] = -A_star / N_amplitude * (usr_i - LV_x_c); 
		user.ncz[usr_k][usr_j][usr_i][1] =  -B_star / N_amplitude * (usr_j - LV_y_c); 
		user.ncz[usr_k][usr_j][usr_i][2] = 0.0;		
}
	}

	if(user.epiendo[usr_k][usr_j][usr_i] == EPI && user.geometry[usr_k][usr_j][usr_i][0] == 2){ // RV epi
		A_star = 1.0/ (double) (RV_x_radius * RV_x_radius);
		B_star = 1.0/ (double) (RV_x_radius * RV_y_radius);
		C_star = 1.0/ (double) (RV_x_radius * RV_z_radius);
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - RV_x_c)  * (usr_i - RV_x_c)  + B_star * B_star * (usr_j - RV_y_c)  * (usr_j - RV_y_c) ) ) ; // this is half of what is written in the notes.
if(N_amplitude > 0.0){
		user.c[usr_k][usr_j][usr_i][0] = - B_star / N_amplitude * (usr_j - RV_y_c); 
		user.c[usr_k][usr_j][usr_i][1] =  A_star / N_amplitude * (usr_i - RV_x_c); 
		user.c[usr_k][usr_j][usr_i][2] = 0.0;
		
		user.ncz[usr_k][usr_j][usr_i][0] = A_star / N_amplitude * (usr_i - RV_x_c); 
		user.ncz[usr_k][usr_j][usr_i][1] =  B_star / N_amplitude * (usr_j - RV_y_c); 
		user.ncz[usr_k][usr_j][usr_i][2] = 0.0;
		
}
	}

		if(user.epiendo[usr_k][usr_j][usr_i] == ENDO && user.geometry[usr_k][usr_j][usr_i][0] == 2){ // RV endo
		A_star = 1.0/ (double) (RV_x_radius2 * RV_x_radius2);
		B_star = 1.0/ (double) (RV_x_radius2 * RV_y_radius2);
		C_star = 1.0/ (double) (RV_x_radius2 * RV_z_radius2);
		N_amplitude = sqrt((double) (A_star * A_star * (usr_i - RV_x_c)  * (usr_i - RV_x_c)  + B_star * B_star * (usr_j - RV_y_c)  * (usr_j - RV_y_c) ) ) ; // this is half of what is written in the notes.
if(N_amplitude > 0.0){
		user.c[usr_k][usr_j][usr_i][0] =  B_star / N_amplitude * (usr_j - RV_y_c);  // multiplied by the sign.
		user.c[usr_k][usr_j][usr_i][1] =  -A_star / N_amplitude * (usr_i - RV_x_c); 
		user.c[usr_k][usr_j][usr_i][2] = 0.0;
		
		user.ncz[usr_k][usr_j][usr_i][0] = -A_star / N_amplitude * (usr_i - RV_x_c); 
		user.ncz[usr_k][usr_j][usr_i][1] =  -B_star / N_amplitude * (usr_j - RV_y_c); 
		user.ncz[usr_k][usr_j][usr_i][2] = 0.0;
}
	}
	
} // end of loops.	

// write
		sprintf(str,"c.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_fibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				fprintf(output,"%10.10F %10.10F %10.10F\n", user.c[usr_k][usr_j][usr_i][0], user.c[usr_k][usr_j][usr_i][1], user.c[usr_k][usr_j][usr_i][2] );
			}
		fclose(output);

		sprintf(str,"ncz.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_fibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				fprintf(output,"%10.10F %10.10F %10.10F\n", user.ncz[usr_k][usr_j][usr_i][0], user.ncz[usr_k][usr_j][usr_i][1], user.ncz[usr_k][usr_j][usr_i][2] );
			}
		fclose(output);
		

	// now do the p vector.
for(usr_k = 0; usr_k < usr_MZ; usr_k++)
for(usr_j = 0; usr_j < usr_MY; usr_j++)
for(usr_i = 0; usr_i < usr_MX; usr_i++){
					user.p[usr_k][usr_j][usr_i][0] = 0.0; 
					user.p[usr_k][usr_j][usr_i][1] =  0.0; 
					user.p[usr_k][usr_j][usr_i][2] = 0.0;	

			if(user.epiendo[usr_k][usr_j][usr_i] == EPI || user.epiendo[usr_k][usr_j][usr_i] == ENDO){
			if(user.epiendo[usr_k][usr_j][usr_i] == EPI)	alpha = -70.0 * M_PI / 180.0; // epi.
			if(user.epiendo[usr_k][usr_j][usr_i] == ENDO)	alpha =  80.0 * M_PI / 180.0; // endo

			tempx1 = user.c[usr_k][usr_j][usr_i][0] * cos(alpha);
			tempy1 = user.c[usr_k][usr_j][usr_i][1] * cos(alpha);
			tempz1 = sin(alpha);
			N_amplitude = sqrt(tempx1 * tempx1 + tempy1 * tempy1 + tempz1 * tempz1);
			if(N_amplitude > 0){
					user.p[usr_k][usr_j][usr_i][0] = tempx1 / N_amplitude; 
					user.p[usr_k][usr_j][usr_i][1] =  tempy1 / N_amplitude; 
					user.p[usr_k][usr_j][usr_i][2] = tempz1 / N_amplitude;
			}else
			{
					user.p[usr_k][usr_j][usr_i][0] = 0.0; 
					user.p[usr_k][usr_j][usr_i][1] =  0.0; 
					user.p[usr_k][usr_j][usr_i][2] = 0.0;
				
			}
			
			}					
}



		sprintf(str,"p.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_fibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
//				if(user.epiendo[usr_k][usr_j][usr_i] > 0)
				fprintf(output,"%1.3F %1.3F %1.3F\n", user.p[usr_k][usr_j][usr_i][0], user.p[usr_k][usr_j][usr_i][1], user.p[usr_k][usr_j][usr_i][2] );
//				else
//				fprintf(output,"%10.10F %10.10F %10.10F\n", 0.0, 0.0, 0.0);					
			}
		fclose(output);


// finally, the fibres.
for(usr_k = 1; usr_k <= usr_MZ-1; usr_k++)
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++){
	user.fibres[usr_k][usr_j][usr_i][0] = 0.0;
	user.fibres[usr_k][usr_j][usr_i][1] = 0.0;	
	user.fibres[usr_k][usr_j][usr_i][2] = 0.0;					

			if(user.epiendo[usr_k][usr_j][usr_i] == EPI || user.epiendo[usr_k][usr_j][usr_i] == ENDO){	
	tempx1 = - ( user.p[usr_k][usr_j][usr_i][0] * user.sheets[usr_k][usr_j][usr_i][0] + user.p[usr_k][usr_j][usr_i][1] * user.sheets[usr_k][usr_j][usr_i][1] + user.p[usr_k][usr_j][usr_i][2] * user.sheets[usr_k][usr_j][usr_i][2]) * user.ncz[usr_k][usr_j][usr_i][0] + (user.ncz[usr_k][usr_j][usr_i][0] * user.sheets[usr_k][usr_j][usr_i][0] + user.ncz[usr_k][usr_j][usr_i][1] * user.sheets[usr_k][usr_j][usr_i][1] + user.ncz[usr_k][usr_j][usr_i][2] * user.sheets[usr_k][usr_j][usr_i][2]) * user.p[usr_k][usr_j][usr_i][0];
	
	tempy1 = - ( user.p[usr_k][usr_j][usr_i][0] * user.sheets[usr_k][usr_j][usr_i][0] + user.p[usr_k][usr_j][usr_i][1] * user.sheets[usr_k][usr_j][usr_i][1] + user.p[usr_k][usr_j][usr_i][2] * user.sheets[usr_k][usr_j][usr_i][2]) * user.ncz[usr_k][usr_j][usr_i][1] + (user.ncz[usr_k][usr_j][usr_i][0] * user.sheets[usr_k][usr_j][usr_i][0] + user.ncz[usr_k][usr_j][usr_i][1] * user.sheets[usr_k][usr_j][usr_i][1] + user.ncz[usr_k][usr_j][usr_i][2] * user.sheets[usr_k][usr_j][usr_i][2]) * user.p[usr_k][usr_j][usr_i][1];
	
	tempz1 = - ( user.p[usr_k][usr_j][usr_i][0] * user.sheets[usr_k][usr_j][usr_i][0] + user.p[usr_k][usr_j][usr_i][1] * user.sheets[usr_k][usr_j][usr_i][1] + user.p[usr_k][usr_j][usr_i][2] * user.sheets[usr_k][usr_j][usr_i][2]) * user.ncz[usr_k][usr_j][usr_i][2] + (user.ncz[usr_k][usr_j][usr_i][0] * user.sheets[usr_k][usr_j][usr_i][0] + user.ncz[usr_k][usr_j][usr_i][1] * user.sheets[usr_k][usr_j][usr_i][1] + user.ncz[usr_k][usr_j][usr_i][2] * user.sheets[usr_k][usr_j][usr_i][2]) * user.p[usr_k][usr_j][usr_i][2];

	N_amplitude = sqrt(tempx1 * tempx1 + tempy1 * tempy1 + tempz1 * tempz1);
	
			if(N_amplitude > 0){	
	user.fibres[usr_k][usr_j][usr_i][0] = tempx1 / N_amplitude;
	user.fibres[usr_k][usr_j][usr_i][1] = tempy1 / N_amplitude;	
	user.fibres[usr_k][usr_j][usr_i][2] = tempz1 / N_amplitude;	
			}
			else
			{
	user.fibres[usr_k][usr_j][usr_i][0] = 0.0;
	user.fibres[usr_k][usr_j][usr_i][1] = 0.0;	
	user.fibres[usr_k][usr_j][usr_i][2] = 0.0;					
			}
			}
}


// test if all fibres are on epi-endo, if all epi-endo is on surface.
for(usr_k = 1; usr_k <= usr_MZ-1; usr_k++)
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++)
if( user.geometry[usr_k][usr_j][usr_i][0] > 0 && (user.geometry[usr_k][usr_j][usr_i][1] == 0 || user.geometry[usr_k][usr_j][usr_i][2] == 0 || user.geometry[usr_k][usr_j][usr_i][3] == 0 || user.geometry[usr_k][usr_j][usr_i][4] == 0 || user.geometry[usr_k][usr_j][usr_i][5] == 0 || user.geometry[usr_k][usr_j][usr_i][6] == 0) )  {
	if(user.epiendo[usr_k][usr_j][usr_i] == 0)
	printf("%d %d %d %d %d\n", usr_i, usr_j, usr_k, 	user.epiendo[usr_k][usr_j][usr_i], user.geometry[usr_k][usr_j][usr_i][0]);
}


/************************************************** Purkinje fibres ****************************************************************************/
// NOTE: AS OF NOW THE PF NETWORK IS A SEPARATE ARRAY: USING USER.GEOMETRY FOR THE PF WILL BOTHER THE FIBRES CALCULATION.
// SORT THIS OUT AT SET UP.

// the first part of this is to generate a HIS system as given in Dux-Santoy 2013, IEEE EMBS Osaka, 3-7 July.
// to do this, examine the geometry using paraview and the epi-endo surfaces, get points coordinates along the HIS, LBB, and RBB.
// Change of plan: 
// I just realised that to keep the HIS bundle simple, I can just do parts of ellipses here.
// just follow the ventricle EPI and ENDO surfaces on the septum.

double theta = 0.0;
int last_i, last_j, last_k;
int last_i2, last_j2, last_k2;
int last_i3, last_j3, last_k3;
double t;

int first_i, first_j, first_k;
int first_i2, first_j2, first_k2;
int first_i3, first_j3, first_k3;

int my_i, my_j, my_k; // my AVNode location I got from paraview by visual inspection.

first_i = first_j = first_k = first_i2 = first_j2 = first_k2 = first_i3 = first_j3 = first_k3 = -1; // impossible values.

// 8 mm depth is where the AV node is.
for(usr_k = (int)(8.0/DX); usr_k <= usr_MZ-1; usr_k++)
for(usr_j = LV_y_c; usr_j < usr_MY-1; usr_j++){	
			// the two ellipses are on the x = LV_x_c location.
			usr_i = LV_x_c;
			theta = atan2( (double)(usr_k - LV_z_c) , (double) (usr_j - LV_y_c) ); // This function computes the arctangent of y/x, but the signs of both arguments are used to determine the quadrant of the result, and x is permitted to be zero. The return value is given in radians and is in the range -pi to pi, inclusive. 

			// part of LV his bundle.
			if( user.epiendo[usr_k][usr_j][usr_i]==ENDO && user.geometry[usr_k][usr_j][usr_i][0]==1 && theta < M_PI * 85.0 / (2.0 * 90.0) ){
				purkinje[usr_k][usr_j][usr_i][0] = 1; // LV HIS bundle is 1.
					purkinje[usr_k][usr_j][usr_i+1][0]=1;
					purkinje[usr_k][usr_j][usr_i-1][0]=1;
					purkinje[usr_k][usr_j+1][usr_i][0]=1;
					purkinje[usr_k][usr_j-1][usr_i][0]=1;
					purkinje[usr_k+1][usr_j][usr_i][0]=1;
					purkinje[usr_k-1][usr_j][usr_i][0]=1;	
					last_i3 = usr_i; last_j3 = usr_j; last_k3 = usr_k;
					if(first_i < 0 && first_j < 0 && first_k < 0){
						first_i = usr_i; first_j = usr_j; first_k = usr_k;
					}
			}
			// part of RV his bundle.
			if( user.epiendo[usr_k][usr_j][usr_i]==EPI && user.geometry[usr_k][usr_j][usr_i][0]==1 &&  theta < M_PI * 55.0 / (2.0 * 90.0) ){
					purkinje[usr_k][usr_j][usr_i][0] = 2; // RV HIS bundle is 2.
					purkinje[usr_k][usr_j][usr_i+1][0]=2;
					purkinje[usr_k][usr_j][usr_i-1][0]=2;
					purkinje[usr_k][usr_j+1][usr_i][0]=2;
					purkinje[usr_k][usr_j-1][usr_i][0]=2;
					purkinje[usr_k+1][usr_j][usr_i][0]=2;
					purkinje[usr_k-1][usr_j][usr_i][0]=2;	
					last_i = usr_i; last_j = usr_j; last_k = usr_k;
					if(first_i2 < 0 && first_j2 < 0 && first_k2 < 0){
						first_i2 = usr_i; first_j2 = usr_j; first_k2 = usr_k;
					}
			}
} // end of HIS bundles part 1.

// draw the line segments from AV node to LBB and RBB using my_i, first_i, and first_i2.
my_i = 80; my_j = 140; my_k = 20; // will only work at DX = 0.4!
for(t = 0; t <= 1.0; t = t + 0.001){

	// LV segment 1.
	usr_i = (int) ( (double)my_i + t * (double)(first_i - my_i)   );
	usr_j = (int) ( (double)my_j + t * (double)(first_j - my_j)   );
	usr_k = (int) ( (double)my_k + t * (double)(first_k - my_k)   );
					purkinje[usr_k][usr_j][usr_i][0] = 1; // LV HIS bundle is 1.
					purkinje[usr_k][usr_j][usr_i+1][0]=1;
					purkinje[usr_k][usr_j][usr_i-1][0]=1;
					purkinje[usr_k][usr_j+1][usr_i][0]=1;
					purkinje[usr_k][usr_j-1][usr_i][0]=1;
					purkinje[usr_k+1][usr_j][usr_i][0]=1;
					purkinje[usr_k-1][usr_j][usr_i][0]=1;	
	// RV segment 2.	
	usr_i = (int) ( (double)my_i + t * (double)(first_i2 - my_i)   );
	usr_j = (int) ( (double)my_j + t * (double)(first_j2 - my_j)   );
	usr_k = (int) ( (double)my_k + t * (double)(first_k2 - my_k)   );
					purkinje[usr_k][usr_j][usr_i][0] = 2; // RV HIS bundle is 2.
					purkinje[usr_k][usr_j][usr_i+1][0]=2;
					purkinje[usr_k][usr_j][usr_i-1][0]=2;
					purkinje[usr_k][usr_j+1][usr_i][0]=2;
					purkinje[usr_k][usr_j-1][usr_i][0]=2;
					purkinje[usr_k+1][usr_j][usr_i][0]=2;
					purkinje[usr_k-1][usr_j][usr_i][0]=2;	
}


// I need 1 point on the RV endo suraface. I got the rough coordinates manually from paraview sphere.
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				if(usr_i > 75 && usr_i < 85 && usr_j > 145 && usr_j < 150 && usr_k > 111 && usr_k < 114 && user.epiendo[usr_k][usr_j][usr_i]==ENDO && user.geometry[usr_k][usr_j][usr_i][0]==2){
					printf("Potential starting point for L-system in RV: %d %d %d\n", usr_i, usr_j, usr_k);
					last_i2 = usr_i; last_j2 = usr_j; last_k2 = usr_k;
				}
			}

			// now draw a line between last and last2.
for(t = 0; t <= 1.0; t = t + 0.001){
usr_i = (int) ( (double)last_i + t * (double)(last_i2 - last_i)   );
usr_j = (int) ( (double)last_j + t * (double)(last_j2 - last_j)   );
usr_k = (int) ( (double)last_k + t * (double)(last_k2 - last_k)   );
				purkinje[usr_k][usr_j][usr_i][0] = 2; // RV HIS bundle is 2.
					purkinje[usr_k][usr_j][usr_i+1][0]=2;
					purkinje[usr_k][usr_j][usr_i-1][0]=2;
					purkinje[usr_k][usr_j+1][usr_i][0]=2;
					purkinje[usr_k][usr_j-1][usr_i][0]=2;
					purkinje[usr_k+1][usr_j][usr_i][0]=2;
					purkinje[usr_k-1][usr_j][usr_i][0]=2;
}

// second part of LV - the semi-circle going to the posterior part of the LV.
for(usr_j = 1; usr_j < usr_MY-1; usr_j++)
for(usr_i = 1; usr_i < usr_MX-1; usr_i++){

	usr_k = (int)(15.0/DX);
	theta = atan2( (double)(usr_j - LV_y_c) , (double) (usr_i - LV_x_c) );
			// part of LV his bundle.
			if( user.epiendo[usr_k][usr_j][usr_i]==ENDO && user.geometry[usr_k][usr_j][usr_i][0]==1 && theta > M_PI * 85.0 / (2.0 * 90.0) ){
				purkinje[usr_k][usr_j][usr_i][0] = 1; // LV HIS bundle is 1.
					purkinje[usr_k][usr_j][usr_i+1][0]=1;
					purkinje[usr_k][usr_j][usr_i-1][0]=1;
					purkinje[usr_k][usr_j+1][usr_i][0]=1;
					purkinje[usr_k][usr_j-1][usr_i][0]=1;
					purkinje[usr_k+1][usr_j][usr_i][0]=1;
					purkinje[usr_k-1][usr_j][usr_i][0]=1;	
					if(first_i3 < 0 && first_j3 < 0 && first_k3 < 0){
						first_i3 = usr_i; first_j3 = usr_j; first_k3 = usr_k;
					}
					
			}
	
}

// now I need these 3 points. Write it on screen and to a text file.
/* at 0.4 = DX, we get this:
LV anterior (i,j,k) format: 80 93 144
LV anterior (i,j,k) format: 36 80 37
RV (i,j,k) format: 84 149 113
*/
printf("LV anterior (i,j,k) format: %d %d %d\n", last_i3, last_j3, last_k3);
printf("LV anterior (i,j,k) format: %d %d %d\n", first_i3, first_j3, first_k3);	
printf("RV (i,j,k) format: %d %d %d\n", last_i2, last_j2, last_k2);

/******************************************************************************************************************************************************************************************/
/******************************************************************************************************************************************************************************************/
/*
LV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX)+MARGIN; 
LV_y_c 	       	= (PetscInt)((double)LV_Y_RADIUS / DX) + MARGIN; 
LV_z_c 	   	= MARGIN;
LV_x_radius	= (PetscInt)((double)LV_X_RADIUS / DX);      
LV_y_radius 	= (PetscInt)((double)LV_Y_RADIUS / DX);      
LV_z_radius 	= (PetscInt)((double)LV_Z_RADIUS / DX);
LV_x_radius2 = (PetscInt)((double)LV_X_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_y_radius2 = (PetscInt)((double)LV_Y_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX); 
LV_z_radius2 = (PetscInt)((double)LV_Z_RADIUS / DX) - (PetscInt)((double)LV_THICKNESS / DX);

RV_x_c 	      	= (PetscInt)((double)LV_X_RADIUS / DX)+MARGIN; 
RV_y_c 	       	= (PetscInt)((double)(LV_Y_RADIUS) / DX) + MARGIN;  // the X and Y coordinates of this ellipsoid are the same as for LV.
RV_z_c 	   	= MARGIN;
RV_x_radius	= (PetscInt)((double)RV_X_RADIUS / DX);      
RV_y_radius = (PetscInt)((double)RV_Y_RADIUS / DX);      
RV_z_radius = (PetscInt)((double)RV_Z_RADIUS / DX);
RV_x_radius2 = (PetscInt)((double)RV_X_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX); 
RV_y_radius2 = (PetscInt)((double)RV_Y_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX); 
RV_z_radius2 = (PetscInt)((double)RV_Z_RADIUS / DX) - (PetscInt)((double)RV_THICKNESS / DX);
*/
// epi and endo surfaces in XYZ format for Meshlab. I write the x, y, z, and the normal components to the half ellipsoids to the xyz file.
		sprintf(str,"purkinjeLV.xyz");
		output = fopen(str,"w");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				if(user.epiendo[usr_k][usr_j][usr_i]==ENDO && user.geometry[usr_k][usr_j][usr_i][0] == 1){ // this is for LV endo.
					// normals.
					tempx1 = -2.0 * ( (float)usr_i - (float) LV_x_c ) / (float)(LV_x_radius2);
					tempy1 = -2.0 * ( (float)usr_j - (float) LV_y_c ) / (float)(LV_y_radius2);
					tempz1 = -2.0 * ( (float)usr_k - (float) LV_z_c ) / (float)(LV_z_radius2);
					fprintf(output, "%f %f %f %f %f %f\n", (float)usr_i, (float)usr_j, (float)usr_k, tempx1, tempy1, tempz1);
				}
			}
		fclose(output);

		sprintf(str,"purkinjeRV.xyz");
		output = fopen(str,"w");
			for(usr_k = MARGIN+1; usr_k < usr_MZ-2; usr_k++)
			for(usr_j = MARGIN+1; usr_j < usr_MY-2; usr_j++)
			for(usr_i = MARGIN+1; usr_i < usr_MX-2; usr_i++){
				if(user.epiendo[usr_k][usr_j][usr_i]==ENDO && user.geometry[usr_k][usr_j][usr_i][0] == 2){ // this is for RV endo - septal side.
					// normals, components of a unit vector.
					tempx1 = -2.0 * ( (float) usr_i - (float)  RV_x_c ) / (float)(RV_x_radius2);
					tempy1 = -2.0 * ( (float)usr_j - (float) RV_y_c ) / (float)(RV_x_radius2);
					tempz1 = -2.0 * ( (float)usr_k - (float) RV_z_c ) / (float)(RV_x_radius2);
					fprintf(output, "%f %f %f %f %f %f\n", (float)usr_i, (float)usr_j, (float)usr_k, tempx1, tempy1, tempz1);
				}		

/*				
				tempx2 = ( (float) usr_i * (float)DX - (float)  RV_x_c * (float) DX ) * ( (float) usr_i * (float)DX - (float)  RV_x_c * (float) DX ) / ( (RV_X_RADIUS-RV_THICKNESS) * (RV_X_RADIUS-RV_THICKNESS) );
				tempy2 = ( (float)usr_j  * (float)DX - (float) RV_y_c * (float) DX ) * ( (float)usr_j  * (float)DX - (float) RV_y_c * (float) DX ) / ((RV_Y_RADIUS-RV_THICKNESS) * (RV_Y_RADIUS-RV_THICKNESS) );
				tempz2 = ( (float)usr_k * (float)DX - (float) RV_z_c * (float) DX ) * ( (float)usr_k * (float)DX - (float) RV_z_c * (float) DX ) / ((RV_Z_RADIUS-RV_THICKNESS) * (RV_Z_RADIUS-RV_THICKNESS) );

				if(user.epiendo[usr_k][usr_j][usr_i]==EPI && user.geometry[usr_k][usr_j][usr_i][0] == 1 && (tempx2 + tempy2 + tempz2 < 1.0) && usr_j > usr_MY/2){ // this is for RV endo.
					// normals, components of a unit vector.
					tempx1 = 2.0 * ( (float) usr_i * (float)DX - (float)  LV_x_c * (float) DX ) / (LV_X_RADIUS-LV_THICKNESS);
					tempy1 = 2.0 * ( (float)usr_j  * (float)DX - (float) LV_y_c * (float) DX ) / (LV_Y_RADIUS-LV_THICKNESS);
					tempz1 = 2.0 * ( (float)usr_k * (float)DX - (float) LV_z_c * (float) DX ) / (LV_Z_RADIUS-LV_THICKNESS);
					fprintf(output, "%f %f %f %f %f %f\n", (float)usr_i * (float)DX, (float)usr_j * (float)DX, (float)usr_k * (float)DX, tempx1, tempy1, tempz1);
				}				
*/				
				
			}

// here, you also want the septum.
		fclose(output);

/*****************************************************************************************************************************************************/			
/*****************************************************************************************************************************************************/			
// see it
for(tissue_type = 1; tissue_type <= 2; tissue_type++){
		sprintf(str,"humanVentricle%d.vtk", tissue_type);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
//		fprintf(output,"SPACING %g %g %g\n", DX, DX, DX);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS huamnVentricle int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				if(tissue_type==user.geometry[usr_k][usr_j][usr_i][0])
				fprintf(output,"%d ", 1 );
				else
				fprintf(output,"%d ", 0 );					
			}
				fprintf(output,"\n");
			}
		fclose(output);
}

for(tissue_type = 1; tissue_type <= 2; tissue_type++){
		sprintf(str,"humanVentriclePF_HIS%d.vtk", tissue_type);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
//		fprintf(output,"SPACING %g %g %g\n", DX, DX, DX);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS huamnVentriclePFHIS int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				if(tissue_type==purkinje[usr_k][usr_j][usr_i][0])
				fprintf(output,"%d ", 1 );
				else
				fprintf(output,"%d ", 0 );					
			}
				fprintf(output,"\n");
			}
		fclose(output);
}


for(tissue_type = 1000; tissue_type <= 4000; tissue_type = tissue_type + 1000){
		sprintf(str,"epiendo%d.vtk", tissue_type);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
//		fprintf(output,"SPACING %g %g %g\n", DX, DX, DX);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_epiendo int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				if(tissue_type==user.epiendo[usr_k][usr_j][usr_i])
				fprintf(output,"%d ", user.epiendo[usr_k][usr_j][usr_i] );
				else
				fprintf(output,"%d ", 0 );					
			}
				fprintf(output,"\n");
			}
		fclose(output);
}


// test the fibres/sheets.
		sprintf(str,"fibres_on_surface.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HV_fibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
				fprintf(output,"%10.10F %10.10F %10.10F\n", user.fibres[usr_k][usr_j][usr_i][0], user.fibres[usr_k][usr_j][usr_i][1], user.fibres[usr_k][usr_j][usr_i][2] );
			}
		fclose(output);

// write the geometry.
/* Write your geometry file, in case you make another program for the simulations. You need the tissue type, neighbours information here. */

		sprintf(str,"truncatedVentricles.geom");
		output = fopen(str,"w");
			for(usr_k = 0; usr_k < usr_MZ; usr_k++)
			for(usr_j = 0; usr_j < usr_MY; usr_j++)
			for(usr_i = 0; usr_i < usr_MX; usr_i++){
fprintf(output,"%d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f\n",usr_k, usr_j, usr_i, user.geometry[usr_k][usr_j][usr_i][0],user.geometry[usr_k][usr_j][usr_i][1], user.geometry[usr_k][usr_j][usr_i][2], user.geometry[usr_k][usr_j][usr_i][3], user.geometry[usr_k][usr_j][usr_i][4], user.geometry[usr_k][usr_j][usr_i][5], user.geometry[usr_k][usr_j][usr_i][6], user.epiendo[usr_k][usr_j][usr_i], user.fibres[usr_k][usr_j][usr_i][0], user.fibres[usr_k][usr_j][usr_i][1], user.fibres[usr_k][usr_j][usr_i][2], user.sheets[usr_k][usr_j][usr_i][0], user.sheets[usr_k][usr_j][usr_i][1], user.sheets[usr_k][usr_j][usr_i][2]);
	}
	fclose(output);



  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

			   VecDestroy(&u); 
			   DMDestroy(&da); 


   PetscFinalize();
   PetscFunctionReturn(0);
}

