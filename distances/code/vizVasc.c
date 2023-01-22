#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.


int main(int argc,char **argv){
	
#include "declarations.c"	

nodeInfo	   *plotThisSegment;

int temp_me_data;

//================================================================================================================
//================================================================================================================
//================================================================================================================
// assign temp_z.
temp_z = 0;
output = fopen("rawSegsRCA.data","r");
if(output==NULL){
	printf("RCA raw data not there yet.\n");
}
else
{
	
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF) temp_z++;
fclose(output);

printf("Number of nodes in RCA tree: %d\n", temp_z);

plotThisSegment = (nodeInfo *) calloc(temp_z, sizeof(nodeInfo) );
output = fopen("rawSegsRCA.data","r");
for(temp_x=0; temp_x < temp_z; temp_x++){
fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order);

	plotThisSegment[temp_x].x_cod = tx;
	plotThisSegment[temp_x].y_cod = ty;
	plotThisSegment[temp_x].z_cod = tz;
	
	plotThisSegment[temp_x].me_data 	= intp1;
	plotThisSegment[temp_x].parent_data 	= intp2;
	
	plotThisSegment[temp_x].radius		= tRad;
	
	plotThisSegment[temp_x].resistance 	= viz_resistivity;	
	plotThisSegment[temp_x].pressure 		= viz_pressure;		
	plotThisSegment[temp_x].perfusion 	= viz_perfusion;
	plotThisSegment[temp_x].me_order 	= order;	
}
fclose(output);


for(temp_x=0; temp_x < temp_z; temp_x++){
		temp_me_data = plotThisSegment[temp_x].me_data;	
		plotThisSegment[temp_x].me_data = temp_x;

	for(temp_y=0; temp_y<temp_z; temp_y++){
		if(temp_me_data==plotThisSegment[temp_y].parent_data){
			plotThisSegment[temp_y].parent_data = temp_x;
		}
	}
		printf("%d %d %d\n", temp_x, plotThisSegment[temp_x].me_data, plotThisSegment[temp_x].parent_data);	
	
}

// make a radius for the tubes here.
	sprintf(str,"segmentsRadius3DRCA.vtk");
	output = fopen(str,"w");
	fprintf(output,"# vtk DataFile Version 3.0\n");
	fprintf(output,"Segments\n");
	fprintf(output,"ASCII\n");
	fprintf(output,"\n");
	fprintf(output,"DATASET POLYDATA\n");
	fprintf(output,"POINTS %d float\n", temp_z);
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%f %f %f\n",plotThisSegment[temp_x].x_cod, plotThisSegment[temp_x].y_cod, plotThisSegment[temp_x].z_cod);
	fprintf(output,"\n");
	fprintf(output,"LINES %d %d\n", temp_z-1, 3*(temp_z-1) ); // you have 2 roots now.
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	if(plotThisSegment[temp_x].parent_data >= 0)
	fprintf(output,"2 %d %d\n",plotThisSegment[temp_x].me_data, plotThisSegment[temp_x].parent_data );
	fprintf(output,"\n");
	fprintf(output, "POINT_DATA %d\n",temp_z);
	fprintf(output, "SCALARS radius double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%20.20f \n", plotThisSegment[temp_x].radius);
	fprintf(output,"\n");
	fprintf(output, "SCALARS resistance double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%20.20f \n", plotThisSegment[temp_x].resistance);
	fprintf(output,"\n");
	fprintf(output, "SCALARS pressure double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%20.20f \n", plotThisSegment[temp_x].pressure);
	fprintf(output,"\n");
	fprintf(output, "SCALARS perfusion double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%lf \n", plotThisSegment[temp_x].perfusion);
	fclose(output);
	free(plotThisSegment);

} // this is if the raw data file exists.
//================================================================================================================	
//================================================================================================================
//================================================================================================================	

// assign temp_z.
temp_z = 0;
output = fopen("rawSegsLA.data","r");
if(output==NULL){
	printf("raw data for LA not there yet.\n");
}
else
{
	
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF) temp_z++;
fclose(output);

printf("Number of nodes in LA tree: %d\n", temp_z);

plotThisSegment = (nodeInfo *) calloc(temp_z, sizeof(nodeInfo) );
output = fopen("rawSegsLA.data","r");
for(temp_x=0; temp_x < temp_z; temp_x++){
fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order);

	plotThisSegment[temp_x].x_cod = tx;
	plotThisSegment[temp_x].y_cod = ty;
	plotThisSegment[temp_x].z_cod = tz;
	
	plotThisSegment[temp_x].me_data 	= intp1;
	plotThisSegment[temp_x].parent_data 	= intp2;
	
	plotThisSegment[temp_x].radius		= tRad;
	
	plotThisSegment[temp_x].resistance 	= viz_resistivity;	
	plotThisSegment[temp_x].pressure 		= viz_pressure;		
	plotThisSegment[temp_x].perfusion 	= viz_perfusion;
	plotThisSegment[temp_x].me_order 	= order;	
}
fclose(output);

for(temp_x=0; temp_x < temp_z; temp_x++){
		temp_me_data = plotThisSegment[temp_x].me_data;	
		plotThisSegment[temp_x].me_data = temp_x;

	for(temp_y=0; temp_y<temp_z; temp_y++){
		if(temp_me_data==plotThisSegment[temp_y].parent_data){
			plotThisSegment[temp_y].parent_data = temp_x;
		}
	}
		printf("%d %d %d\n", temp_x, plotThisSegment[temp_x].me_data, plotThisSegment[temp_x].parent_data);	
	
}

// make a radius for the tubes here.
	sprintf(str,"segmentsRadius3DLA.vtk");
	output = fopen(str,"w");
	fprintf(output,"# vtk DataFile Version 3.0\n");
	fprintf(output,"Segments\n");
	fprintf(output,"ASCII\n");
	fprintf(output,"\n");
	fprintf(output,"DATASET POLYDATA\n");
	fprintf(output,"POINTS %d float\n", temp_z);
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%f %f %f\n",plotThisSegment[temp_x].x_cod, plotThisSegment[temp_x].y_cod, plotThisSegment[temp_x].z_cod);
	fprintf(output,"\n");
	fprintf(output,"LINES %d %d\n", temp_z-1, 3*(temp_z-1) ); // you have 2 roots now.
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	if(plotThisSegment[temp_x].parent_data >= 0)
	fprintf(output,"2 %d %d\n",plotThisSegment[temp_x].me_data, plotThisSegment[temp_x].parent_data );
	fprintf(output,"\n");
	fprintf(output, "POINT_DATA %d\n",temp_z);
	fprintf(output, "SCALARS radius double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%20.20f \n", plotThisSegment[temp_x].radius);
	fprintf(output,"\n");
	fprintf(output, "SCALARS resistance double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%20.20f \n", plotThisSegment[temp_x].resistance);
	fprintf(output,"\n");
	fprintf(output, "SCALARS pressure double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%20.20f \n", plotThisSegment[temp_x].pressure);
	fprintf(output,"\n");
	fprintf(output, "SCALARS perfusion double\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%lf \n", plotThisSegment[temp_x].perfusion);
	fclose(output);
	free(plotThisSegment);
	
}

//================================================================================================================
//================================================================================================================

	return 0;
}
