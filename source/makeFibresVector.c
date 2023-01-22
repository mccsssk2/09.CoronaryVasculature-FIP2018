/*
Sanjay Kharche.
SRK codes.
21 Jan 2017.

This is a single proc run for writing the 3 components of the fibre vector.

Run the srkHumanVentricleGeom.c to get a geometry at given resolution.
Run the srk3d.c for the Poisson equation to get stable components of fibre vector (3 runs, one for each component. See Kuhl 2014 and other papers.)
Then read the last bin from each of these 3 directories, and write it to a vtk and geom file. Revise the geom as well.

Inputs:
truncatedVentricles.geom
fx/my205.bin (last file after the diffusion is steady).
fy/my205.bin
fy/my205.bin

Outputs:

revised geom file.
revised fibres VTK file.
*/
static char help[] = "Mouse 2015 3D SAN measurements, post-processing.\n";

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

#define NEQ         3                     /* number of reaction variables. This is FK3 */

#define DX           0.4                 /* X internode spacing                       */
#define DY           0.4                 /* Y internode spacing                       */
#define DZ	     0.4

// pass this as data in simulations to make it dynamic. On the other hand, the resolution is set in the geometry creation.
#define usr_MX 160
#define usr_MY 215
#define usr_MZ 185

/* The stencil is u0 = x,y,z; , u1 = y+1, u2 = x+1, u3 = y-1, u4 = x-1, u5 = z+1, u6 = z - 1 */
#define NBS    7    // 3D arrays have 7 units. Itself, and 6 surrounding it in a standard 1st order FD stencil.

/* User-defined data structures and routines           */
/* AppCtx: used by FormIFunction() and FormIJacobian() */
typedef struct {
  DM        da; // DM instance
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
  PetscInt  ****geometry; // 3D models have 4D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1, 5 is z+1, 6 is z-1
  PetscInt      ***epiendo       ;
  // fibres and sheets arrays.  
PetscReal ****sheets;
PetscReal ****fibres;   // srk3d.c Posson solver has 3D arrays. I need a 4D array here.
} AppCtx;

		//! Byte swap int
		int32_t swap_int32( int32_t val )
		{
		    val = ((val << 8) & 0xFF00FF00) | ((val >> 8) & 0xFF00FF ); 
		    return (val << 16) | ((val >> 16) & 0xFFFF);
		}

		float ReverseFloat( const float inFloat )
		{
		   float retVal;
		   char *floatToConvert = ( char* ) & inFloat;
		   char *returnFloat = ( char* ) & retVal;

		   // swap the bytes into a temporary buffer
		   returnFloat[0] = floatToConvert[3];
		   returnFloat[1] = floatToConvert[2];
		   returnFloat[2] = floatToConvert[1];
		   returnFloat[3] = floatToConvert[0];

		   return retVal;
		}


int main(int argc,char **argv)
{

  PetscInt   file_Counter = 0;
  FILE      *geometry, *output;
  char       str[1000]; // so you do not keep creating and destroying this.
int         intx, inty, intz, intu0,  intu1 , intu2 , intu3 , intu4 , intu5 , intu6 , intu7 ;

/* Program reduced to use of FD with colouring. */
  Vec            u;              /* solution, residual vectors. You need the r for SNES */
  AppCtx         user;                  /* user-defined work context  */
  DM             da;
  PetscInt       usr_i, usr_j, usr_k, time_int;
  PetscScalar    ***u_localptr;
  PetscInt start, end;
  PetscReal f1, f2, f3, s1, s2, s3;
  PetscInt       mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z;         
         

  PetscInitialize(&argc,&argv,(char*)0,help);

// fileCounter loop start
// dir 4 needs:
start = 40;
end   = 41;
for(file_Counter = start; file_Counter < end; file_Counter++){ // this is where I saw 1 activation of the SAN model.

  /* Initialize user application context */
  user.da           = NULL;

   DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_BOX,usr_MX,usr_MY,usr_MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,NULL,&da);  /* this decides that u is non-ghosted 1st order star/box stencil vector. */
   user.da = da;
   DMDASetUniformCoordinates(da,0,usr_MX, 0, usr_MY, 0, usr_MZ); // this is to set up the dimensions so that we get a good vts file.
  DMCreateGlobalVector(da,&u);      // at this time.
	   DMDAGetCorners(da,&mybase_x,&mybase_y,&mybase_z,&mysize_x,&mysize_y,&mysize_z);  // For uniformity.

        /* Declare the geometry memory. This declares memory on each proc. */
		user.geometry = (PetscInt  ****) calloc(usr_MZ, sizeof(PetscInt ***));
		user.epiendo    = (PetscInt  ***)  calloc(usr_MZ, sizeof(PetscInt ** ));	
		user.fibres	   = (PetscReal  ****)  calloc(usr_MZ, sizeof(PetscReal *** ));			
		for(usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.geometry[usr_k] = (PetscInt***)  calloc(usr_MY,sizeof(PetscInt **));
			user.epiendo[usr_k] = (PetscInt**)  calloc(usr_MY,sizeof(PetscInt *));
			user.fibres[usr_k] 	    = (PetscReal***)  calloc(usr_MY,sizeof(PetscReal **));			
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.geometry[usr_k][usr_j] = (PetscInt **) calloc(usr_MX,sizeof(PetscInt*));
				user.epiendo[usr_k][usr_j] = (PetscInt *) calloc(usr_MX,sizeof(PetscInt));
				user.fibres[usr_k][usr_j] = (PetscReal **) calloc(usr_MX,sizeof(PetscReal*));	
				for (usr_i = 0; usr_i < usr_MX; usr_i++){
					user.geometry[usr_k][usr_j][usr_i] = (PetscInt*) calloc(NBS,sizeof(PetscInt));				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal *) calloc(3,sizeof(PetscReal*));	
				}
			}
		}

				/* read in your geometry here.  		*/
				/* This reads the ASCII geometry on each proc. 
				The file is opened for reading on each proc.    */
				geometry = fopen("truncatedVentricles.geom","r");
			while(fscanf(geometry,"%d %d %d %d  %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf ",&intz, &inty, &intx, &intu0, &intu1, &intu2, &intu3, &intu4, &intu5, &intu6, &intu7, &f1, &f2, &f3, &s1, &s2, &s3)!=EOF){
				usr_i = (PetscInt)intx - mybase_x; usr_j = (PetscInt)inty - mybase_y; usr_k = (PetscInt)intz - mybase_z;
				if(usr_k>=0&&usr_k<mysize_z&&usr_j>=0&&usr_j<mysize_y&&usr_i>=0&&usr_i<mysize_x){
				user.geometry[usr_k][usr_j][usr_i][0] = (PetscInt)intu0;
				user.geometry[usr_k][usr_j][usr_i][1] = (PetscInt)intu1;
				user.geometry[usr_k][usr_j][usr_i][2] = (PetscInt)intu2;
				user.geometry[usr_k][usr_j][usr_i][3] = (PetscInt)intu3;
				user.geometry[usr_k][usr_j][usr_i][4] = (PetscInt)intu4;
				user.geometry[usr_k][usr_j][usr_i][5] = (PetscInt)intu5;
				user.geometry[usr_k][usr_j][usr_i][6] = (PetscInt)intu6;
				user.epiendo[usr_k][usr_j][usr_i]         = (PetscInt)intu7;		
				} // end of reading geomtry.
			}
			fclose(geometry);


   // inputs are Petsc binary files.
			PetscViewer viewer_in_fx;
			sprintf(str,"fx/my_3d%d_1.bin",file_Counter);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_READ,&viewer_in_fx);
			VecLoad(u,viewer_in_fx);
			PetscViewerDestroy(&viewer_in_fx);
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	 user.fibres[usr_k][usr_j][usr_i][0] = (float)u_localptr[usr_k][usr_j][usr_i]; /* LHS is my 3D data array from appContext, RHS is PetSc's internal representation. */
    DMDAVecRestoreArray(da,u,&u_localptr);

			PetscViewer viewer_in_fy;
			sprintf(str,"fy/my_3d%d_2.bin",file_Counter);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_READ,&viewer_in_fy);
			VecLoad(u,viewer_in_fy);
			PetscViewerDestroy(&viewer_in_fy);
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	 user.fibres[usr_k][usr_j][usr_i][1] = (float)u_localptr[usr_k][usr_j][usr_i]; /* LHS is my 3D data array from appContext, RHS is PetSc's internal representation. */
    DMDAVecRestoreArray(da,u,&u_localptr);

			PetscViewer viewer_in_fz;
			sprintf(str,"fz/my_3d%d_3.bin",file_Counter);
			PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_READ,&viewer_in_fz);
			VecLoad(u,viewer_in_fz);
			PetscViewerDestroy(&viewer_in_fz);
    DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	 user.fibres[usr_k][usr_j][usr_i][2] = (float)u_localptr[usr_k][usr_j][usr_i]; /* LHS is my 3D data array from appContext, RHS is PetSc's internal representation. */
    DMDAVecRestoreArray(da,u,&u_localptr);
    
    // writing the fibres throughout the model.

		sprintf(str,"HVF%d.vtk", file_Counter);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX, usr_MY, usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HVFibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	{
		if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i]  > 0)
		fprintf(output,"%1.3f %1.3f %1.3f\n", user.fibres[usr_k][usr_j][usr_i][0], user.fibres[usr_k][usr_j][usr_i][1], user.fibres[usr_k][usr_j][usr_i][2]);
		else
		fprintf(output,"%1.3f %1.3f %1.3f\n", 0.0, 0.0, 0.0);			
	}
		fclose(output);


		sprintf(str,"HVF2_%d.vtk", file_Counter);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX, usr_MY, usr_MZ);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*usr_MZ);
		fprintf(output,"SCALARS HVFibres float 3\n");
		fprintf(output,"LOOKUP_TABLE default\n");
	for(usr_k = 0; usr_k < usr_MZ; usr_k++)
	for(usr_j = 0; usr_j < usr_MY; usr_j++)
	for(usr_i = 0; usr_i < usr_MX; usr_i++)
	{
		if(user.geometry[usr_k][usr_j][usr_i][0] > 0)
		fprintf(output,"%1.3f %1.3f %1.3f\n", user.fibres[usr_k][usr_j][usr_i][0], user.fibres[usr_k][usr_j][usr_i][1], user.fibres[usr_k][usr_j][usr_i][2]);
		else
		fprintf(output,"%1.3f %1.3f %1.3f\n", 0.0, 0.0, 0.0);			
	}
		fclose(output);
		
		
/**********************************************************************************************/

// this is in your time loop for now.
   VecDestroy(&u); 
   DMDestroy(&da); 

 } // fileCounter loop end

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Free work space.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

   PetscFinalize();
   PetscFunctionReturn(0);
}

