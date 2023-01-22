// fileNum now may be set to a constant 0. The expt binary runs in ALL_* and I am probably using xargs or parallel or xxargs.
 
// RCA part. =================================================================================================================================
// assign temp_z.
sprintf(str, "rawSegsRCA%d.data", fileNum);
temp_z = 0;
output = fopen(str,"r");

if(output==NULL){
printf("no rawsegs RCA file. moving on.\n");
}
else
{
	
	
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF) temp_z++;
fclose(output);

printf("Number of nodes in RCA sub-tree: %d\n", temp_z);

plotThisSegment = (nodeInfo *) calloc(temp_z, sizeof(nodeInfo) );
output = fopen(str,"r");
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
}

// make a radius for the tubes here.
	sprintf(str,"segmentsRadius3DRCA%d.vtk", fileNum);
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
	fprintf(output,"LINES %d %d\n", temp_z-1, 3*(temp_z-1) );
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
	fprintf(output,"\n");
	fprintf(output, "SCALARS order int\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%d \n", plotThisSegment[temp_x].me_order);
	fclose(output);

free(plotThisSegment)	;

} // end of fopen if.

 // LA ====================================================================================================
// assign temp_z.
sprintf(str,"rawSegsLA%d.data", fileNum);
temp_z = 0;
output = fopen(str,"r");
if(output==NULL){
printf("no raw LA file here, moving on.\n");
}
else
{
	
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF) temp_z++;
fclose(output);

printf("Number of nodes in LA sub-tree: %d\n", temp_z);

plotThisSegment = (nodeInfo *) calloc(temp_z, sizeof(nodeInfo) );
output = fopen(str,"r");
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
}

// make a radius for the tubes here.
	sprintf(str,"segmentsRadius3DLA%d.vtk", fileNum);
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
	fprintf(output,"LINES %d %d\n", temp_z-1, 3*(temp_z-1) );
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
	fprintf(output,"\n");
	fprintf(output, "SCALARS order int\n");
	fprintf(output, "LOOKUP_TABLE default\n");
	for(temp_x = 0; temp_x < temp_z; temp_x++)
	fprintf(output,"%d \n", plotThisSegment[temp_x].me_order);	
	fclose(output);

free(plotThisSegment)	;
	

} // end of LA fopen.

// ========================================================================================================
// ========================================================================================================

