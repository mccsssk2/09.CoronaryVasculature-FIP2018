/*
SRK, LHSC, 6 February 2017.

This file is used after the fractal-tree simulation is finished.
It takes the output produced by the python scripts, and puts it into the geom file.

All paths, constants, and other data are hard coded into this program as of now. Get it working, then make it seamless with the other code.

The files given by the python code are:
srkPurkinje-line_xyz.txt			x, y, z coordinates of the nodes. I am assuming that the line number is the node number, numbering starts at 0.
srkPurkinje-line_ien.txt  			starting and ending nodes of each line segment.
srkPurkinje-line_endnodes.txt		The terminal points of the PF, the data is in node numbers.


Inputs to the program:
a) above 3 txt files.
b) the geom file.

Outputs:
a) geom file.
b) VTK files for each tissue type.
*/

static char help[] = "Synthetic 3D geometry for human ventricle Purkinje fibres (truncated, not the whole ellipsoid), 1 process job.\n";

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

#define NCOODS 3 // x, y, z are the 3 coordinates in the 2D nodes array.

/* User-defined data structures and routines           */
/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da; // DM instance. This is required.
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
/* all skPetsc programs use geomtry as: 0 = itself, 1 = y+1, 2 = x+1, 3 = y -1, 4 = x - 1, 5 = z + 1, 6 = z - 1 */
  PetscInt     ****geometry    ; // 3D models have 4D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1.
  PetscInt      ***epiendo       ;

// fibres and sheets arrays.
PetscReal **nodes;
PetscReal **nodesRV;
PetscReal ****sheets;
PetscReal ****fibres;
 
} AppCtx;

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
	  char 		str[10000]; // so you do not keep creating and destroying this.
          PetscInt       temp_x, temp_y, temp_z;
         FILE               *input, *output, *geometry;
	int         intx, inty, intz, intu0,  intu1 , intu2 , intu3 , intu4 , intu5 , intu6 , intu7 ;
        PetscReal f1, f2, f3, s1, s2, s3;         

         int ****purkinje; // the Purkinje fibre array. 
         
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
		for(usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.geometry[usr_k] = (PetscInt***)  	calloc(usr_MY,sizeof(PetscInt **));
			user.epiendo[usr_k] = (PetscInt**)  		calloc(usr_MY,sizeof(PetscInt *));
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.geometry[usr_k][usr_j] = (PetscInt **) 	calloc(usr_MX,sizeof(PetscInt*));
				user.epiendo[usr_k][usr_j] = (PetscInt *) 		calloc(usr_MX,sizeof(PetscInt));				
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.geometry[usr_k][usr_j][usr_i] = (PetscInt*) 	calloc(NBS,sizeof(PetscInt));									
				}
			}
		}

		user.sheets     	= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));		
		user.fibres 		= (PetscReal  ****) calloc(usr_MZ, sizeof(PetscReal ***));
		
		for(usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.sheets[usr_k] 	= (PetscReal***)  calloc(usr_MY,sizeof(PetscReal **));
			user.fibres[usr_k] 	= (PetscReal***)  calloc(usr_MY,sizeof(PetscReal **));
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.sheets[usr_k][usr_j] 		= (PetscReal**)  calloc(usr_MX,sizeof(PetscReal *));
				user.fibres[usr_k][usr_j] 	= (PetscReal**)  calloc(usr_MX,sizeof(PetscReal *));
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.sheets[usr_k][usr_j][usr_i] = (PetscReal*) calloc(3, sizeof(PetscReal));
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal*) calloc(3, sizeof(PetscReal));									
				}
			}
		}


		// read in the current geometry. This is the geometry produced by the srkHumanVentricleGeom.c program. It may not have all the 
// computed fibres. The fibres and sheets are irrelevant to this calculation.
				/* read in your geometry here.  		*/
				/* This reads the ASCII geometry on each proc. 
				The file is opened for reading on each proc.    */
				geometry = fopen("truncatedVentricles.geom","r");
			while(fscanf(geometry,"%d %d %d %d  %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf ",&intz, &inty, &intx, &intu0, &intu1, &intu2, &intu3, &intu4, &intu5, &intu6, &intu7, &f1, &f2, &f3, &s1, &s2, &s3)!=EOF){
				usr_i = (PetscInt)intx - mybase_x; usr_j = (PetscInt)inty - mybase_y; usr_k = (PetscInt)intz - mybase_z;
				if(usr_k>=0&&usr_k<mysize_z&&usr_j>=0&&usr_j<mysize_y&&usr_i>=0&&usr_i<mysize_x){
				user.geometry[usr_k][usr_j][usr_i][0] = (PetscInt)intu0;
				} // end of reading geomtry.
			}
			fclose(geometry);

			int LV_nodes; double t;
/******************************** LV PF NODES *****************************************************************************************************************************************/			
			// now work out how many lines in the txt file for PF node coordinates.
			// LV nodes.
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/LV/srkPurkinje-line_xyz.txt");
			input = fopen(str,"r");			
			LV_nodes = 0;
			while(fscanf(input,"%lf %lf %lf",&tempx1, &tempy1, &tempz1)!=EOF){
				LV_nodes++;
			}
			fclose(input);
			
			printf("The number of LV nodes is %d\n", LV_nodes);
			// declare a 2D array of LV_nodes size.
			user.nodes = (PetscReal  **) calloc(LV_nodes, sizeof(PetscReal *));
			for (usr_i = 0; usr_i < LV_nodes; usr_i++) user.nodes[usr_i] = (PetscReal*) calloc(3, sizeof(PetscReal));
			
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/LV/srkPurkinje-line_xyz.txt");
			input = fopen(str,"r");			
			LV_nodes = 0;
			while(fscanf(input,"%lf %lf %lf",&tempx1, &tempy1, &tempz1)!=EOF){
				user.nodes[LV_nodes][0] = tempx1/DX;  user.nodes[LV_nodes][1] = tempy1/DX; user.nodes[LV_nodes][2] = tempz1/DX; 
//				printf("%d %f %f %f\n", LV_nodes, user.nodes[LV_nodes][0], user.nodes[LV_nodes][1], user.nodes[LV_nodes][2]);
				LV_nodes++;
			}
			fclose(input);
			
			// now look up the end nodes of each line segment,
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/LV/srkPurkinje-line_ien.txt");
			input = fopen(str,"r");
			while(fscanf(input,"%d %d",&intu1, &intu2)!=EOF){
				/*
				printf("P1: %f %f %f. ", LV_nodes, user.nodes[intu1][0], user.nodes[intu1][1], user.nodes[intu1][2]);
				printf("P2: %f %f %f\n", LV_nodes, user.nodes[intu2][0], user.nodes[intu2][1], user.nodes[intu2][2]);								
				*/
				// parameteric equation of the line segment.
				for(t=0;t<=1;t=t+0.01){
					tempx1 = user.nodes[intu1][0] + t * (user.nodes[intu2][0] - user.nodes[intu1][0]);
					tempy1 = user.nodes[intu1][1] + t * (user.nodes[intu2][1] - user.nodes[intu1][1]);					
					tempz1 = user.nodes[intu1][2] + t * (user.nodes[intu2][2] - user.nodes[intu1][2]);					
					// convert this to usr_i, usr_j, usr_k.
					usr_i = (int)(tempx1 );	usr_j = (int)(tempy1 );	usr_k = (int)(tempz1 );
					if(usr_i> 0 && usr_i < usr_MX-1 && usr_j > 0 && usr_j < usr_MY-1 && usr_k > 0 && usr_k < usr_MZ-1){
					user.geometry[usr_k][usr_j][usr_i][0] = 3; // LV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j][usr_i+1][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j][usr_i-1][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j+1][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j-1][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k+1][usr_j][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k-1][usr_j][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					}
				}
			}
			fclose(input);
	// now apply the end nodes. this is tissue type 4.		
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/LV/srkPurkinje-line_endnodes.txt");
			input = fopen(str,"r");
			while(fscanf(input,"%d ",&intu1)!=EOF){
					tempx1 = user.nodes[intu1][0];
					tempy1 = user.nodes[intu1][1];					
					tempz1 = user.nodes[intu1][2];					
					// convert this to usr_i, usr_j, usr_k.
					usr_i = (int)(tempx1 );	usr_j = (int)(tempy1 );	usr_k = (int)(tempz1 );
					if(usr_i> 0 && usr_i < usr_MX-1 && usr_j > 0 && usr_j < usr_MY-1 && usr_k > 0 && usr_k < usr_MZ-1){					
					user.geometry[usr_k][usr_j][usr_i][0] = 4; // LV Purkinje-myocardial junctions are tissue type 4 in LV.
					user.geometry[usr_k][usr_j][usr_i+1][0] = 4; 
					user.geometry[usr_k][usr_j][usr_i-1][0] = 4; 
					user.geometry[usr_k][usr_j+1][usr_i][0] = 4; 
					user.geometry[usr_k][usr_j-1][usr_i][0] = 4; 
					user.geometry[usr_k][usr_j+1][usr_i][0] = 4; 
					user.geometry[usr_k][usr_j-1][usr_i][0] = 4; 
					}

			}

/******************************** RV PF NODES *****************************************************************************************************************************************/			
// RV nodes.
int RV_nodes; 
			// RV nodes.
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/RV/srkPurkinje-line_xyz.txt");
			input = fopen(str,"r");			
			RV_nodes = 0;
			while(fscanf(input,"%lf %lf %lf",&tempx1, &tempy1, &tempz1)!=EOF){
				RV_nodes++;
			}
			fclose(input);
			
			printf("The number of RV nodes is %d\n",RV_nodes);
			// declare a 2D array of RV_nodes size.
			user.nodesRV = (PetscReal  **) calloc(RV_nodes, sizeof(PetscReal *));
			for (usr_i = 0; usr_i < RV_nodes; usr_i++) user.nodesRV[usr_i] = (PetscReal*) calloc(3, sizeof(PetscReal));
			
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/RV/srkPurkinje-line_xyz.txt");
			input = fopen(str,"r");			
			RV_nodes = 0;
			while(fscanf(input,"%lf %lf %lf",&tempx1, &tempy1, &tempz1)!=EOF){
				user.nodesRV[RV_nodes][0] = tempx1/DX;  user.nodesRV[RV_nodes][1] = tempy1/DX; user.nodesRV[RV_nodes][2] = tempz1/DX; 
//				printf("%d %f %f %f\n", RV_nodes, user.nodesRV[RV_nodes][0], user.nodesRV[RV_nodes][1], user.nodesRV[RV_nodes][2]);
				RV_nodes++;
			}
			fclose(input);
			
			// now look up the end nodes of each line segment,
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/RV/srkPurkinje-line_ien.txt");
			input = fopen(str,"r");
			while(fscanf(input,"%d %d",&intu1, &intu2)!=EOF){
				/*
				printf("P1: %f %f %f. ", RV_nodes, user.nodesRV[intu1][0], user.nodesRV[intu1][1], user.nodesRV[intu1][2]);
				printf("P2: %f %f %f\n", RV_nodes, user.nodesRV[intu2][0], user.nodesRV[intu2][1], user.nodesRV[intu2][2]);								
				*/
				// parameteric equation of the line segment.
				for(t=0;t<=1;t=t+0.01){
					tempx1 = user.nodesRV[intu1][0] + t * (user.nodesRV[intu2][0] - user.nodesRV[intu1][0]);
					tempy1 = user.nodesRV[intu1][1] + t * (user.nodesRV[intu2][1] - user.nodesRV[intu1][1]);					
					tempz1 = user.nodesRV[intu1][2] + t * (user.nodesRV[intu2][2] - user.nodesRV[intu1][2]);					
					// convert this to usr_i, usr_j, usr_k.
					usr_i = (int)(tempx1);	usr_j = (int)(tempy1);	usr_k = (int)(tempz1);
					if(usr_i> 0 && usr_i < usr_MX-1 && usr_j > 0 && usr_j < usr_MY-1 && usr_k > 0 && usr_k < usr_MZ-1){					
					user.geometry[usr_k][usr_j][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j][usr_i+1][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j][usr_i-1][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j+1][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j-1][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j+1][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					user.geometry[usr_k][usr_j-1][usr_i][0] = 3; // RV Purkinje fibres are tissue type 3.
					}
				}
			}
			fclose(input);
	// now apply the end nodes. this is tissue type 4.		
			sprintf(str,"/home/kharches/projects/McIntyre/SyntheticHearts/1.EllenKuhl/fractal-tree/RV/srkPurkinje-line_endnodes.txt");
			input = fopen(str,"r");
			while(fscanf(input,"%d ",&intu1)!=EOF){
					tempx1 = user.nodesRV[intu1][0];
					tempy1 = user.nodesRV[intu1][1];					
					tempz1 = user.nodesRV[intu1][2];					
					// convert this to usr_i, usr_j, usr_k.
					usr_i = (int)(tempx1);	usr_j = (int)(tempy1);	usr_k = (int)(tempz1);
					if(usr_i> 0 && usr_i < usr_MX-1 && usr_j > 0 && usr_j < usr_MY-1 && usr_k > 0 && usr_k < usr_MZ-1){					
					user.geometry[usr_k][usr_j][usr_i][0] = 4; // LV Purkinje-myocardial junctions are tissue type 4 in RV.
					user.geometry[usr_k][usr_j][usr_i+1][0] = 4; 
					user.geometry[usr_k][usr_j][usr_i-1][0] = 4; 
					user.geometry[usr_k][usr_j+1][usr_i][0] = 4; 
					user.geometry[usr_k][usr_j-1][usr_i][0] = 4; 
					user.geometry[usr_k][usr_j+1][usr_i][0] = 4; 
					user.geometry[usr_k][usr_j-1][usr_i][0] = 4; 
					}
			}

			
/******************************** end PF NODES *****************************************************************************************************************************************/			

// do the neighbours.
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

			
/*****************************************************************************************************************************************************/			
/*****************************************************************************************************************************************************/			
// see it, now we have 4 tissue types.
for(tissue_type = 1; tissue_type <= 4; tissue_type++){
		sprintf(str,"humanVentricleWithPF%d.vtk", tissue_type);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,usr_MZ);
//		fprintf(output,"SPACING %g %g %g\n", DX, DX, DX);
		fprintf(output,"SPACING %g %g %g\n", 1, 1, 1);
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

// write the geometry.
/* Write your geometry file, in case you make another program for the simulations. You need the tissue type, neighbours information here. */
// this is the geometry with the PF tissue types in.

		sprintf(str,"truncatedVentricles2.geom");
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

