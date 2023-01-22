/*
Oct 4, 2017.
Distance calculation. Measurement of vasculature data.
Inputs are inorder and preorder, and associated Nodes files.
*/

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

#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.

int main(int argc,char **argv){
	
#include "declarations.c"	

double resolution1, resolution2, resolution3;
XYZ 	     new_loc;
double 	distToSurface_avbdr, distToSurface_LVCHAMBER, distToSurface_RVCHAMBER, distToSurface_OUTSIDE, dx;
XYZ 	     avbdr_closest, avbdr_normal, LVCHAMBER_closest, LVCHAMBER_normal, RVCHAMBER_closest, RVCHAMBER_normal, OUTSIDE_closest, OUTSIDE_normal;

// read in the rootRCA and rootLCA here.

// declare three arrays that will have total perfusion.
	
resolution1 	= atof(argv[1]);  // mm

mysize_x1=(int)((double)TOTAL_X/(double)resolution1);  mysize_y1=(int)((double)TOTAL_Y/(double)resolution1);   mysize_z1=(int)((double)TOTAL_Z/(double)resolution1);

arr ***res1;

// declare your structs.
		res1												= (arr***)	calloc(mysize_z1,sizeof(arr **));
		for(usr_k = 0; usr_k < mysize_z1; usr_k++){
			res1[usr_k] 									= (arr**)  calloc(mysize_y1,sizeof(arr* ));
			for (usr_j = 0; usr_j < mysize_y1; usr_j++){
				res1[usr_k][usr_j] 							= (arr *) calloc(mysize_x1,sizeof(arr));
			}
		}

// define epi endo. In fact, define how far the box is from the endo surface.
// this is a time consuming step. Calculate this distance, and save it forever and ever.

	printf("The xmax, ymax, zmax are: %d %d %d\n", mysize_x1, mysize_y1, mysize_z1);

		dx = resolution1;
		for(usr_k = 0; usr_k < mysize_z1; usr_k++)
			for (usr_j = 0; usr_j < mysize_y1; usr_j++)
				for(usr_i = 0; usr_i < mysize_x1; usr_i++){
				new_loc.x = dx * (double)usr_i+dx/2.0; new_loc.y = dx * (double)usr_j+dx/2.0; new_loc.z = dx * (double)usr_k+dx/2.0;
				res1[usr_k][usr_j][usr_i].distanceToEndo = 0.0;
					if(isNodeInTissue(new_loc.x, new_loc.y, new_loc.z, 1.0)==1){
						distToSurface_LVCHAMBER 	= closestLVCHAMBER	(new_loc, &LVCHAMBER_closest, &LVCHAMBER_normal, 1.);
						distToSurface_RVCHAMBER 	= closestRVCHAMBER	(new_loc, &RVCHAMBER_closest, &RVCHAMBER_normal, 1.);

						if(distToSurface_LVCHAMBER < distToSurface_RVCHAMBER)
							res1[usr_k][usr_j][usr_i].distanceToEndo = distToSurface_LVCHAMBER;
						else
							res1[usr_k][usr_j][usr_i].distanceToEndo = distToSurface_RVCHAMBER;				
							
//							if(res1[usr_k][usr_j][usr_i].distanceToEndo > 14.0)
							printf("%f %d %d %d\n", res1[usr_k][usr_j][usr_i].distanceToEndo, usr_i, usr_j, usr_k);
					}
				}

		printf("did the distances, now writing.\n");				
		// write a data file as well.
		sprintf(str,"DISTANCE%d_%d_%d.DATUM", mysize_x1, mysize_y1, mysize_z1);
		output = fopen(str, "w");	
		for(usr_k = 0; usr_k < mysize_z1; usr_k++)
			for (usr_j = 0; usr_j < mysize_y1; usr_j++)
				for(usr_i = 0; usr_i < mysize_x1; usr_i++)
					fprintf(output, "%d %d %d %2.9f\n", usr_i, usr_j, usr_k, res1[usr_k][usr_j][usr_i].distanceToEndo);
		fclose(output);
		
		dx = resolution1;
		sprintf(str,"distance%d_%d_%d.vtk", mysize_x1, mysize_y1, mysize_z1);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",mysize_x1, mysize_y1, mysize_z1);
		fprintf(output,"SPACING %f %f %f\n", dx, dx, dx);
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",mysize_x1*mysize_y1*mysize_z1);
		fprintf(output,"SCALARS distance float 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
		for(usr_k = 0; usr_k < mysize_z1; usr_k++)
			for (usr_j = 0; usr_j < mysize_y1; usr_j++){
				for(usr_i = 0; usr_i < mysize_x1; usr_i++)
		fprintf(output,"%2.9f ", res1[usr_k][usr_j][usr_i].distanceToEndo);
		fprintf(output,"\n");
		}
		fclose(output);


return 0;

} // end of main function.


