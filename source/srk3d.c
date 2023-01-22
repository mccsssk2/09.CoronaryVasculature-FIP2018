/*

This needs rewriting. Use ex14.c in snes/

Sanjay Kharche.
Diffusion of theta angle with Dirichlet boundary conditions.
19 January 2017.

Issues:
1) You have no boundary condition at usr_k = 5, i.e. the region close to the base. So you will lose information.
2) The Dirichlet problem is not being solved accurately. Use KSP or other solver for that.

20 January 2017. Do the fibres and sheets according to Kuhl 2014 and her suggestions.
*/

static char help[] = "srk3D.c, skCode2016 (LHSC). A step in working out the synthetic boundary conditions.\n";

// each of these headers have a whole bunch of headers themselves. See the manual/documentation of the header files.
#include <petscdm.h>
#include <petscdmda.h>
#include <petscts.h>
#include <petscts.h>
#include <petscvec.h>
#include <petscsys.h>

#define MAXSTEPS 5000
/* my standard sundials headers, constants, and macros that will fit Petsc. Petsc cannot handle cvodes, at least not installations. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include <cvode/cvode.h>                  	/* prototypes for CVODE fcts., consts.  	*/
#include <nvector/nvector_serial.h>       	/* serial N_Vector types, fcts., macros 	*/
#include <cvode/cvode_dense.h>            	/* prototype for CVDense                		*/
#include <sundials/sundials_dense.h>      /* definitions DlsMat DENSE_ELEM        	*/
#include <sundials/sundials_types.h>      /* definition of type realtype          		*/
#include <sundials/sundials_math.h>
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ        */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ        */
#define RTOL        RCONST(1.0e-12)   	  /* scalar relative tolerance            */
#define ATOL        RCONST(1.0e-6)        /* scalar absolute tolerance components */

// this is in microM
#define DX            0.400                 /* X internode spacing                       */
#define DY            0.400                 /* Y internode spacing                       */
#define DZ	        0.400 		  /* Z internode spacing                       */

// pass this as data in simulations to make it dynamic. On the other hand, the resolution is set in the geometry creation.
#define usr_MX 160
#define usr_MY 215
#define usr_MZ 185

#define DELTAT        0.25		// ms, the time interval when I must have output.
#define diffusion       0.10             // diffusion constant. mm2/ms if your length is in mm and time is in ms. I dont know how much the diffusion is.

#define NEQ        			3
#define NBS        			7    // This is for the FD neighbourshood scheme. 3D arrays have 7 neighbours, itself, and 6 surrounding it in a standard 1st order FD stencil.
#define DATA     			7    // This model has 3 bits of data.
#define NUMPARAMS    	13                   /* ODE parameters. Often or all the time I want to specify these as input */

#define MARGIN 5

/* User-defined data structures and routines. The Petsc application context defined by skCode. */
typedef struct {
/* Application context. */
  DM        da;
/* All model data such as geomtery, fibers, cell type data, extracellular distributions must go through the context. */
  PetscInt  ****geometry; // 3D models have 3D array geometry. 0 = itself, 1 is y+1, 2 is x+1, 3 is y-1, 4 is x-1, 5 is z+1 and 6 is z-1
  PetscInt      ***epiendo       ;
  
  // fibres and sheets arrays.
PetscReal ***sheets;
PetscReal ***fibres;
} AppCtx;

extern PetscErrorCode FormIFunction(TS,PetscReal,Vec,Vec,Vec,void*);

int main(int argc,char **argv)
{
	  TS             	ts;                   			/* nonlinear solver 			*/
	  Vec            	u, r;                  		/* solution, residual vectors 	*/
	  Mat            	J	, Jmf = NULL  ;   	/* Jacobian matrices 			*/
//	  PetscInt       maxsteps = 1000;      	/* iterations for convergence 	*/
	  DM             	da;
	  AppCtx         user;              			/* user-defined work context 	*/
	  SNES           snes;
	// sk variables.
	  PetscInt       mybase_x, mybase_y, mybase_z, mysize_x, mysize_y, mysize_z, usr_i, usr_j, usr_k, time_int, total_time_int;
	  PetscScalar    ***u_localptr;
	  double 	my_x, my_y, my_z;
	  char 		str[1000]; // so you do not keep creating and destroying this.
	  double   	***usr_u0_loc, ***usr_u1_loc, ***usr_u2_loc;
	  realtype       usr_time, usr_t, usr_tout;
          PetscInt       file_Counter = 0, intx, inty, intz;
         FILE       *geometry;
         PetscInt intu1, intu2, intu3, intu4, intu5, intu6, intu7, intu0;
         PetscReal f1, f2, f3, s1, s2, s3;

	  PetscInitialize(&argc,&argv,(char*)0,help);
	
	time_int = 0; total_time_int = 50000;
	for(time_int = 0; time_int < total_time_int; time_int++){ // I decided to do 2 seconds, which takes 1.6 hours on 4 procs.

	     /* Initialize user application context */
	      user.da           = NULL;
	     // the -11, -11 becomes usr_MX, usr_MY when I have #define usr_MX and #define usr_MY (potentially #define usr_MZ) put those values in place of -11,-11
   DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,usr_MX,usr_MY,usr_MZ, PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL,NULL,&da);  /* this decides that u is non-ghosted 1st order star stencil vector. */
	     user.da = da		;
	     DMCreateGlobalVector(da,&u); 
	     VecDuplicate(u,&r)		; // this is the residual vector, which has the same properties as the state vector u.
	   DMDAGetCorners(da,&mybase_x,&mybase_y,&mybase_z,&mysize_x,&mysize_y,&mysize_z);  // you need this for running the Sundials loops.

	     if(time_int<1){
		     
		usr_u0_loc = (double ***) calloc(mysize_z, sizeof(double**));
	for (usr_k = 0; usr_k < mysize_z; usr_k++){ 
				usr_u0_loc[usr_k] = (double **) calloc(mysize_y, sizeof(double*));
				/* the usr_mysize_ and usr_mybase_ arrays are not the same on all processors. */
					for (usr_j = 0; usr_j < mysize_y; usr_j++){ 
						usr_u0_loc[usr_k][usr_j] = (double*) calloc(mysize_x,sizeof(double));
					}
	}
        /* Declare the geometry memory. This declares memory on each proc. */
		user.geometry = (PetscInt  ****) calloc(usr_MZ, sizeof(PetscInt ***));
		user.epiendo    = (PetscInt  ***)  calloc(usr_MZ, sizeof(PetscInt ** ));	
		user.fibres	   = (PetscReal  ***)  calloc(usr_MZ, sizeof(PetscReal ** ));			
		for(usr_k = 0; usr_k < usr_MZ; usr_k++){
			user.geometry[usr_k] = (PetscInt***)  calloc(usr_MY,sizeof(PetscInt **));
			user.epiendo[usr_k] = (PetscInt**)  calloc(usr_MY,sizeof(PetscInt *));
			user.fibres[usr_k] 	    = (PetscReal**)  calloc(usr_MY,sizeof(PetscReal *));			
			for (usr_j = 0; usr_j < usr_MY; usr_j++){
				user.geometry[usr_k][usr_j] = (PetscInt **) calloc(usr_MX,sizeof(PetscInt*));
				user.epiendo[usr_k][usr_j] = (PetscInt *) calloc(usr_MX,sizeof(PetscInt));
				user.fibres[usr_k][usr_j] = (PetscReal *) calloc(usr_MX,sizeof(PetscReal));	
				for (usr_i = 0; usr_i < usr_MX; usr_i++)
					user.geometry[usr_k][usr_j][usr_i] = (PetscInt*) calloc(NBS,sizeof(PetscInt));				
			}
		}

	// read in geometry here. geomtry is in a specific srk format with extension .geom.
				geometry = fopen("truncatedVentricles.geom","r");
			while(fscanf(geometry,"%d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf ",&intz, &inty, &intx, &intu0, &intu1, &intu2, &intu3, &intu4, &intu5, &intu6, &intu7, &f1, &f2, &f3, &s1, &s2, &s3)!=EOF){
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
if(atoi(argv[1])==1){				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal)f1; // these are the initial conditions for this Poisson solver.
if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] > 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)f1;
}

if(atoi(argv[1])==2){				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal)f2; // these are the initial conditions for this Poisson solver.
if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] > 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)f2;
}

if(atoi(argv[1])==3){				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal)f3; // these are the initial conditions for this Poisson solver.
if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] > 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)f3;
}

if(atoi(argv[1])==4){				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal)s1; // these are the initial conditions for this Poisson solver.
if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] > 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)s1;
}

if(atoi(argv[1])==5){				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal)s2; // these are the initial conditions for this Poisson solver.
if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] > 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)s2;
}

if(atoi(argv[1])==6){				
					user.fibres[usr_k][usr_j][usr_i] = (PetscReal)s3; // these are the initial conditions for this Poisson solver.
if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] > 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)s3;
}

if(atoi(argv[1]) < 1 ||atoi(argv[1]) > 6){ printf("Argv1 is from 1 to 6, for fx, fy, fz, sx, sy, sz. Exiting. Try again.\n"); exit(1); }

if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i] == 0) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)0.0;
if(user.geometry[usr_k][usr_j][usr_i][0] == 0 											) usr_u0_loc[usr_k][usr_j][usr_i] = (PetscReal)(-10.0);	
				} // end of reading geomtry.
			}
			fclose(geometry);
			printf("%d %d %d %d %d %d\n", mysize_x, mysize_y, mysize_z, usr_MX, usr_MY, usr_MZ);

	} // end of time_int < 1 ****************************************************************************************************************************************

// exit(0); // just testing reading the geometry.

	/* boundary conditions                                                             */
		for(usr_k = 0; usr_k < mysize_z; usr_k++)
		for(usr_j  = 0; usr_j  < mysize_y; usr_j++)
		for(usr_i  = 0; usr_i  < mysize_x; usr_i++) {
		if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i]==1000)   usr_u0_loc[usr_k][usr_j][usr_i] = user.fibres[usr_k][usr_j][usr_i];
		if(user.geometry[usr_k][usr_j][usr_i][0] > 0 && user.epiendo[usr_k][usr_j][usr_i]==2000)  usr_u0_loc[usr_k][usr_j][usr_i] = user.fibres[usr_k][usr_j][usr_i];	
		}

			// Petsc part, solve for 1 time step from t = 0 to t = DELTAT.
			//     Create timestepping solver context
		   	TSCreate(PETSC_COMM_WORLD,&ts)		; 
		   	TSSetProblemType(ts,TS_NONLINEAR)		; 
		   	TSSetType(ts,TSBEULER)			; 
		   	TSSetDM(ts,da)				; 
		   	TSSetIFunction(ts,r,FormIFunction,&user)	; 
		   	TSSetDuration(ts,MAXSTEPS,DELTAT)		; // this needs to be the time step/output step of my application.
			/*
			Input Parameter
				ts 	- the time-step context
				eftopt 	- exact final time option
			 TS_EXACTFINALTIME_STEPOVER    - Don't do anything if final time is exceeded
			 TS_EXACTFINALTIME_INTERPOLATE - Interpolate back to final time
			 TS_EXACTFINALTIME_MATCHSTEP   - Adapt final time step to match the final time
			*/
			   TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); // get the solution at the exact time of end. This is new to LHSC SRK*.C codes.
			//     Set initial conditions ********************************************
    	DMDAVecGetArray(da,u,&u_localptr);
	for(usr_k = 0; usr_k < mysize_z; usr_k++) for(usr_j = 0; usr_j < mysize_y; usr_j++) for(usr_i = 0; usr_i < mysize_x; usr_i++)
	u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i] = usr_u0_loc[usr_k][usr_j][usr_i];
    	DMDAVecRestoreArray(da,u,&u_localptr);
                        //  ******************************************************************************
			   TSSetSolution(ts,u); 
			   TSSetInitialTimeStep(ts,0.0,DELTAT); // if convergence fails, then you can reduce DELTAT by several fractions.
				//   Set Jacobian evaluation routine using colouring.
			   DMSetMatType(da,MATAIJ); 
			   DMCreateMatrix(da,&J); 
			   TSGetSNES(ts,&snes); 
			   MatCreateSNESMF(snes,&Jmf); 
			   /* Use coloring to compute  finite difference J efficiently */
			   SNESSetJacobian(snes,Jmf,J,SNESComputeJacobianDefaultColor,0); 
				//   Sets various TS parameters from user options
			   TSSetFromOptions(ts); 
				// get a solution from 0 to DELTAT
			   TSSolve(ts,u); 
				// get the u into my non PetSc array.
			    DMDAVecGetArray(da,u,&u_localptr);
			   	 /* LHS is my 2D array, RHS is PetSc's */
    				DMDAVecGetArray(da,u,&u_localptr);
    				for(usr_k = 0; usr_k < mysize_z; usr_k++) for(usr_j = 0; usr_j < mysize_y; usr_j++) for(usr_i = 0; usr_i < mysize_x; usr_i++) 
    				usr_u0_loc[usr_k][usr_j][usr_i] = u_localptr[mybase_z + usr_k][mybase_y + usr_j][mybase_x + usr_i];
    				DMDAVecRestoreArray(da,u,&u_localptr);

				// do a vts view.
				if(time_int % 100 == 0){
							sprintf(str,"my_3d%d.vts", file_Counter);
							PetscViewer viewer;  PetscViewerCreate(PETSC_COMM_WORLD, &viewer); PetscViewerSetType(viewer, PETSCVIEWERVTK);
							PetscViewerFileSetName(viewer, str); VecView(u, viewer); PetscViewerDestroy(&viewer);
				// binary file output.
		                sprintf(str,"my_3d%d_%d.bin",file_Counter, atoi(argv[1]));
				PetscViewer viewer2;
				PetscViewerBinaryOpen(PETSC_COMM_WORLD,str,FILE_MODE_WRITE,&viewer2);
				VecView(u,viewer2);
				PetscViewerDestroy(&viewer2);					
				file_Counter++;
				}

			   MatDestroy(&J); 
			   MatDestroy(&Jmf); 
			   VecDestroy(&u); 
			   VecDestroy(&r); 
			   TSDestroy(&ts); 
			   DMDestroy(&da); 

			// track time
			usr_time = DELTAT * time_int;

} // end of time loop.

// necessary.
   PetscFinalize();
  PetscFunctionReturn(0);
} // end of main function.

/* --------------------------------------------------------------------- */
/*
  FormIFunction = Udot - RHSFunction, RHS minus RHSFunction.
*/
PetscErrorCode FormIFunction(TS ts,PetscReal t,Vec U,Vec Udot,Vec F,void *ctx)
{
 AppCtx         *user=(AppCtx*)ctx;
  DM             da   = (DM)user->da;
  PetscInt       i,j,k,xs,ys,zs,xm,ym,zm,gj, gi, gk;
  PetscReal      hx,hy,hz,sx,sy,sz;
  PetscReal      U0, U1, U2, U3, U4, U5, U6;
  PetscScalar    uxx,uyy,uzz, ***uarray, ***f, ***udot;
  PetscScalar    ux, uy, uz;
  Vec            localU;

PetscReal dd;

  PetscFunctionBeginUser;
   DMGetLocalVector(da,&localU);

  /*
     Scatter ghost points to local vector,using the 2-step process
        DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
     By placing code between these two statements, computations can be
     done while messages are in transition.
  */
   DMGlobalToLocalBegin(da,U,INSERT_VALUES,localU); 
   DMGlobalToLocalEnd(da,U,INSERT_VALUES,localU); 
  /* Get pointers to vector data */
   DMDAVecGetArray(da,localU,&uarray); 
   DMDAVecGetArray(da,F,&f); 
   DMDAVecGetArray(da,Udot,&udot); 
  /* Get local grid boundaries */
   DMDAGetCorners(da,&xs,&ys,&zs,&xm,&ym,&zm); 
/*
This calculation of the RHS, or udot - RHS as they call it is
wrong. I started with ex15.c in the ts examples. So fix this.
Numerically, it is stable. Still some work to be done.
*/
  /* Compute function over the locally owned part of the grid */
  for (k=zs; k<zs+zm; k++)
  for (j=ys; j<ys+ym; j++)
    for (i=xs; i<xs+xm; i++)
	f[k][j][i] = 0.0; // so that i do not have to deal with empty space.

    for (k=zs; k<zs+zm; k++)
    for (j=ys; j<ys+ym; j++) 
    for (i=xs; i<xs+xm; i++) {

// the g's' go from start of box to end of box in each direction.
gk = k - zs;
gj = j - ys;
gi = i - xs;

  hx = DX; sx = diffusion/(hx*hx); // atrium and adipose tissue.
  hy = DY; sy = diffusion/(hy*hy);
  hz = DZ; sz = diffusion/(hz*hz);
  
		/* This is still hacky. I am not happy with the way boundaries are broken down into separate components. */
//		if(user->geometry[gk][gj][gi][0]>0){
      if (i == 0 || j == 0 || i == usr_MX-1 || j == usr_MY-1 || k == 0 || k == usr_MZ-1 ) {
/* Neumann BC. Necessary evil at the walls of the model. I am missing something. */
// surfaces
          	        if (i == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k][j][i+1];  // y-z plane at x = 0
		else if (i == usr_MX-1)  		f[k][j][i] = uarray[k][j][i] - uarray[k][j][i-1];  // y-z plane at x = usr_MX-1
		else if (j ==  0) 			f[k][j][i] = uarray[k][j][i] - uarray[k][j+1][i]; 
		else if (j == usr_MY-1) 		f[k][j][i] = uarray[k][j][i] - uarray[k][j-1][i]; 
          	else if (k == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k+1][j][i]; 
          	else if (k == usr_MZ-1)  	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j][i]; 

// edges
          	  else if (i == 0                 &&   j == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k][j+1][i+1]; 
          	  else if (i == usr_MX-1   &&   j == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k][j+1][i-1];
          	  else if (i == 0                 &&   j == usr_MY-1) 	f[k][j][i] = uarray[k][j][i] - uarray[k][j-1][i+1]; 
           	  else if (i == usr_MX-1   &&   j == usr_MY-1)  	f[k][j][i] = uarray[k][j][i] - uarray[k][j-1][i-1]; 
          	  else if (i == 0                 &&   k == 0)    		f[k][j][i] = uarray[k][j][i] - uarray[k+1][j][i+1]; 
          	  else if (i == usr_MX-1   &&   k == 0)    		f[k][j][i] = uarray[k][j][i] - uarray[k+1][j][i-1]; 
          	  else if (i == 0                 &&   k == usr_MZ-1) 	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j][i+1]; 
           	  else if (i == usr_MX-1   &&   k == usr_MZ-1)  	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j][i-1]; 
          	  else if (k == 0                &&   j == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k+1][j][i+1]; 
          	  else if (k == usr_MZ-1  &&   j == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k+1][j][i-1];
          	  else if (k == 0                &&   j == usr_MY-1) 	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j][i+1]; 
           	  else if (k == usr_MZ-1  &&   j == usr_MY-1)  	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j][i-1]; 
// points
          	else if (i == 0               && j == 0               && k == 0)    			f[k][j][i] = uarray[k][j][i] - uarray[k+1][j+1][i+1]; 
		else if (i == usr_MX-1 && j == 0               && k == 0)      		f[k][j][i] = uarray[k][j][i] - uarray[k+1][j+1][i-1]; 
		else if (i == 0               && j == usr_MY-1 && k == 0) 			f[k][j][i] = uarray[k][j][i] - uarray[k+1][j-1][i+1]; 
		else if (i == usr_MX-1 && j == usr_MY-1 && k == 0) 			f[k][j][i] = uarray[k][j][i] - uarray[k+1][j-1][i-1]; 
          	else if (i == 0               && j == 0               && k == usr_MY-1)    	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j+1][i+1]; 
		else if (i == usr_MX-1 && j == 0               && k == usr_MY-1)       f[k][j][i] = uarray[k][j][i] - uarray[k-1][j+1][i-1]; 
		else if (i == 0               && j == usr_MY-1 && k == usr_MY-1) 	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j-1][i+1]; 
		else if (i == usr_MX-1 && j == usr_MY-1 && k == usr_MY-1) 	f[k][j][i] = uarray[k][j][i] - uarray[k-1][j-1][i-1]; 

      } else if(i>0&&j>0&&k>0&&i<usr_MX-1&&j<usr_MY-1&&k<usr_MZ-1){

if(user->geometry[gk][gj][gi][0]  > 0 && user->epiendo[gk][gj][gi]  == 0){    	      
			                         U1 = uarray[k][j+1][i]; 
		     U4 = uarray[k][j][i-1]; U0 = uarray[k][j  ][i]; U2 = uarray[k][j][i+1]; 
		                             U3 = uarray[k][j-1][i];

								     U5 = uarray[k+1][j][i]; 
								     U6 = uarray[k-1][j][i]; 
}
else
	if(user->geometry[gk][gj][gi][0]  > 0 && user->epiendo[gk][gj][gi]  == 1000){    	      
	U0 = U1 = U2 = U3 = U4 = U5 = U6 = 0.1;		
}
else
		if(user->geometry[gk][gj][gi][0]  > 0 && user->epiendo[gk][gj][gi]  == 2000){    	      
	U0 = U1 = U2 = U3 = U4 = U5 = U6 = 0.9;		
}
else
	U0 = U1 = U2 = U3 = U4 = U5 = U6 = -0.1;

			uxx = (-2.0*U0 + U4 + U2);
			uyy = (-2.0*U0 + U3 + U1);
			uzz = (-2.0*U0 + U5 + U6);
		f[k][j][i] = udot[k][j][i] - (uxx*sx + uyy*sy + uzz*sz ); // the udot is not zero

            } // end of if statement to keep it inside the box.
    } // end of for loops.

  /* Restore vectors */
   DMDAVecRestoreArray(da,localU,&uarray); 

   DMDAVecRestoreArray(da,F,&f);
   DMDAVecRestoreArray(da,Udot,&udot); 
   DMRestoreLocalVector(da,&localU); 
  PetscFunctionReturn(0);}



