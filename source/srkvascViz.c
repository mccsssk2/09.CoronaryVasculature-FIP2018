/*
13 May 2017.
The srkVasculature2D.c/3D.c programs write .data files. Read those data files and do visualisation, statistics, in this program.
The visual output is for VTK legacy format.
All variables in srkVasculature2D/3D are available in this code due to the declarations.c.
*/

#include "srkVasculature.h"

int main(int argc,char **argv){
#include "declarations.c"

double perfMin;

perfMin = 100000000000000000.0;

// this loop is embarassingly parallel.
for(accepted_cood = 0; accepted_cood < atoi(argv[1]); accepted_cood=accepted_cood+4){

// read in the data.
sprintf(str,"hVasc2D%06d.data", accepted_cood);
output = fopen(str,"r");
// what goes into the data file: node_number, srkNodes_me, srkNodes_parent, this_node_is_a_leaf, strahler_me,  me_coods, parent_coods, radius, resistivity, pressure, perfusion.

while( fscanf(output, "%d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &viz_temp_x, &viz_srkNodes0, &viz_srkNodes1, &viz_ThisNodeIsALeaf, &viz_srkStrahlerNumberOfNode, &viz_coodsMe0, &viz_coodsMe1, &viz_coodsMe2, &viz_coodsParent0, &viz_coodsParent1, &viz_coodsParent2, &viz_radiusSegment,
&viz_resistivity, &viz_pressure, &viz_perfusion, &viz_downResistance)!=EOF){

 temp_x  							= viz_temp_x;  
 srkNodes[temp_x][0] 				= viz_srkNodes0; 
 srkNodes[temp_x][1]  				= viz_srkNodes1; 
 ThisNodeIsALeaf[temp_x]  			= viz_ThisNodeIsALeaf; 
 srkStrahlerNumberOfNode[temp_x][0]   = viz_srkStrahlerNumberOfNode; 
 coodsMe[temp_x][0]  				= viz_coodsMe0;  
 coodsMe[temp_x][1] 					= viz_coodsMe1;
 coodsMe[temp_x][2]  				= viz_coodsMe2; 
 coodsParent[temp_x][0] 				= viz_coodsParent0;  
 coodsParent[temp_x][1] 				= viz_coodsParent1;  
 coodsParent[temp_x][2]				= viz_coodsParent2;  
radiusSegment[temp_x] 				= viz_radiusSegment;
 resistivity[temp_x] 					= viz_resistivity;  
 pressure[temp_x]					= viz_pressure;  
 perfusion[temp_x] 					= viz_perfusion;
if(perfMin > viz_perfusion) perfMin = viz_perfusion;
 downResistance[temp_x]				= viz_downResistance;
} // end of while loop.
fclose(output);



// do some calculations here.
for(usr_j = 0; usr_j < usr_MY; usr_j++) for(usr_i = 0; usr_i < usr_MX; usr_i++){ user.vasculature[usr_j][usr_i] = user.vasc_output[usr_j][usr_i] = 0.0; }

for(temp_x = 0; temp_x < nuNodes-1; temp_x++){	
	for(pT = 0; pT <= 1.0; pT = pT + 0.01){
		mex = (int)( (coodsMe[temp_x][0] + pT * (coodsParent[temp_x][0] - coodsMe[temp_x][0]))/DX ); 
		mey = (int)( (coodsMe[temp_x][1] + pT * (coodsParent[temp_x][1] - coodsMe[temp_x][1]))/DX );
		mez = 0; 
if(mex>1&&mex<usr_MX-2&&mey>1&&mey<usr_MY-2){
/*
		user.vasculature[mey][mex] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey+1][mex] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey-1][mex] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey][mex+1] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey][mex-1] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey+1][mex+1] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey+1][mex-1] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey-1][mex+1] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
		user.vasculature[mey-1][mex-1] = 1.0 * radiusSegment[temp_x]/radiusSegment[0];  // centre line
*/


/*
		user.vasculature[mey][mex] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey+1][mex] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey-1][mex] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey][mex+1] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey][mex-1] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey+1][mex+1] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey+1][mex-1] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey-1][mex+1] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
		user.vasculature[mey-1][mex-1] = 1.0 * (perfusion[temp_x]-perfMin);  // centre line
*/

/*
		user.vasculature[mey][mex] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey+1][mex] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey-1][mex] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey][mex+1] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey][mex-1] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey+1][mex+1] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey+1][mex-1] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey-1][mex+1] = 1.0 * (resistivity[temp_x] );  // centre line
		user.vasculature[mey-1][mex-1] = 1.0 * (resistivity[temp_x] );  // centre line
*/
		user.vasculature[mey][mex] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey+1][mex] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey-1][mex] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey][mex+1] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey][mex-1] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey+1][mex+1] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey+1][mex-1] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey-1][mex+1] = 1.0 * (pressure[temp_x] );  // centre line
		user.vasculature[mey-1][mex-1] = 1.0 * (pressure[temp_x] );  // centre line

printf("%d %f \n", temp_x, pressure[temp_x]);

}

//		user.vasculature[mey][mex] = 1.0 * pressure[temp_x];  // centre line, for perfusion, resistance, ...

	} 
}


// write VTK files.
// init geometry painted using strahler numbering.
/*
		sprintf(str,"hVasc2D_geom.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,1);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*1);
		fprintf(output,"SCALARS hVessels int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++)
				fprintf(output,"%d ",  (int)user.vasculature[usr_j][usr_i] );
				fprintf(output,"\n");
			}
		fclose(output);
*/


		sprintf(str,"hVasc2D_rad%06d.vtk", accepted_cood);
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,1);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*1);
		fprintf(output,"SCALARS radius float 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++)
				fprintf(output,"%f ", (double)user.vasculature[usr_j][usr_i] );
				fprintf(output,"\n");
			}
		fclose(output);
		
		/*
		sprintf(str,"hVasc2D_Init_spheres.vtk");
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 3.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",usr_MX,usr_MY,1);
		fprintf(output,"SPACING 1 1 1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",usr_MX*usr_MY*1);
		fprintf(output,"SCALARS hSpheres int 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_j = 0; usr_j < usr_MY; usr_j++){
			for(usr_i = 0; usr_i < usr_MX; usr_i++) fprintf(output,"%d ", user.spheres[usr_j][usr_i]);
				fprintf(output,"\n");
			}
		fclose(output);
		*/

/*______________________________________________________________________________________________*/
/*______________________________________________________________________________________________*/
/*______________________________________________________________________________________________*/

} // end of accepted_cood big loop.

   PetscFinalize();
   PetscFunctionReturn(0);
} // end of main.


