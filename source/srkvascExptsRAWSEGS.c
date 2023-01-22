/*
Oct 26 2017.
First, just calculate FD for control case using 100 instances.
Make a grid at res1 and res2.
Add perfusion in each box.
raw histogram
norm to 100 or however many you have.
normal to max
RD.
FD.


Reading the tree is expensive, so just read it and do the experiments on it: that is part 1 of post-processing.
Part 2 is the FD calculation.
*/

/*
Oct 26 2017.
First, just calculate FD for control case using 100 instances.
Make a grid at res1 and res2.
Add perfusion in each box.
raw histogram
norm to 100 or however many you have.
normal to max
RD.
FD.


Reading the tree is expensive, so just read it and do the experiments on it: that is part 1 of post-processing.
Part 2 is the FD calculation.
*/

#define histNum 100

#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.

int main(int argc,char **argv){
	
#include "declarations.c"	
int fileNum;

nodeInfo	   *plotThisSegment;
double dx, dx2, dx3;

double res1PerfusionCountLV[12], res1PerfusionCountRV[6]; // each LV and RV are divided into layers of 1 mm each.

double ***res2, ***res3;
double ***fres1, ***fres2, ***fres3; // res1 is 1 mm, res2 is 0.1 mm, res3 is 0.2 mm.

double minF2, maxF2, minF3, maxF3;
int countRes2, countRes3;

double meanRes2, meanRes3, totalRes2, totalRes3, total_flow, meanres2, meanres3;

double deltaD_res2, deltaD_res3;

double hist[2][histNum+1]; // this histogram is for the masses.
int myfi_counter;
int temptemp;
double total_counts, hist_mean, hist_total, var_total, sd1, sd6;

// I want the perfusion at a given resolution, resolution1.

// higher resolution.
dx2 = 0.1;
mysize_x2=(int)((double)TOTAL_X/(double)dx2);  mysize_y2=(int)((double)TOTAL_Y/(double)dx2); mysize_z2=(int)((double)TOTAL_Z/(double)dx2);
printf("a dx2, The xmax, ymax, zmax are: %d %d %d\n", mysize_x2, mysize_y2, mysize_z2);
// declare your structs.
		res2												= (double***)	calloc(mysize_z2,sizeof(double **));
		fres2												= (double***)	calloc(mysize_z2,sizeof(double **));		
		for(usr_k = 0; usr_k < mysize_z2; usr_k++){
			res2[usr_k] 									= (double**)  calloc(mysize_y2,sizeof(double* ));
			fres2[usr_k] 									= (double**)  calloc(mysize_y2,sizeof(double* ));			
			for (usr_j = 0; usr_j < mysize_y2; usr_j++){
				res2[usr_k][usr_j] 							= (double *) calloc(mysize_x2,sizeof(double));
				fres2[usr_k][usr_j] 							= (double *) calloc(mysize_x2,sizeof(double));				
			}
		}

// in the middle resolution.		
dx3 = 0.5; // mm		
mysize_x3=(int)((double)TOTAL_X/(double)dx3);  mysize_y3=(int)((double)TOTAL_Y/(double)dx3); mysize_z3=(int)((double)TOTAL_Z/(double)dx3);
printf("at dx3, The xmax, ymax, zmax are: %d %d %d\n", mysize_x3, mysize_y3, mysize_z3);
// declare your structs.
		res3												= (double***)	calloc(mysize_z3, sizeof(double **));
		fres3												= (double***)	calloc(mysize_z3, sizeof(double **));		
		for(usr_k = 0; usr_k < mysize_z3; usr_k++){
			res3[usr_k] 									= (double**)  calloc(mysize_y3, sizeof(double* ));
			fres3[usr_k] 									= (double**)  calloc(mysize_y3, sizeof(double* ));			
			for (usr_j = 0; usr_j < mysize_y3; usr_j++){
				res3[usr_k][usr_j] 							= (double *) calloc(mysize_x3, sizeof(double));
				fres3[usr_k][usr_j] 							= (double *) calloc(mysize_x3, sizeof(double));				
			}
		}

for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++){ res2[usr_k][usr_j][usr_i] = 0; }
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++){ res3[usr_k][usr_j][usr_i]= 0;  }

// ================== declarations/allocations finished. ===================================================================================================
// ================== declarations/allocations finished. ===================================================================================================

countRes2 = countRes3 = 0;

for(fileNum = 1; fileNum <= 120; fileNum++){

sprintf(str, "ALL_%d/rawSegsRCA.data", fileNum );
output = fopen(str,"r");
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF){
	if(order==6 /* && viz_perfusion < 0.01 && viz_perfusion > 0.0001 */ ){
usr_i = (int)floor( (tx / dx2) ); usr_j = (int)floor( (ty / dx2) ); usr_k = (int)floor( (tz / dx2) );
res2[usr_k][usr_j][usr_i] = res2[usr_k][usr_j][usr_i] + viz_perfusion;

usr_i = (int)floor( (tx / dx3) ); usr_j = (int)floor( (ty / dx3) ); usr_k = (int)floor( (tz / dx3) );
res3[usr_k][usr_j][usr_i] = res3[usr_k][usr_j][usr_i] + viz_perfusion;
	}
}
fclose(output);

sprintf(str, "ALL_%d/rawSegsLA.data", fileNum );
output = fopen(str,"r");
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF){
		if(order==6 /* && viz_perfusion < 0.01 && viz_perfusion > 0.0001 */ ){
usr_i = (int)floor( (tx / dx2) ); usr_j = (int)floor( (ty / dx2) ); usr_k = (int)floor( (tz / dx2) );
res2[usr_k][usr_j][usr_i] = res2[usr_k][usr_j][usr_i] + viz_perfusion;

usr_i = (int)floor( (tx / dx3) ); usr_j = (int)floor( (ty / dx3) ); usr_k = (int)floor( (tz / dx3) );
res3[usr_k][usr_j][usr_i] = res3[usr_k][usr_j][usr_i] + viz_perfusion;
		}
}
fclose(output);

} // end of fileNum loop.

// ==========================================================================================================================
// ==========================================================================================================================
// ==========================================================================================================================

// I have not hacked the FD/RD in my code yet: so I am going to set it up for MATLAB.
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
res2[usr_k][usr_j][usr_i] = res2[usr_k][usr_j][usr_i]/120.0;
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
res3[usr_k][usr_j][usr_i] = res3[usr_k][usr_j][usr_i]/120.0;

output = fopen("order6.data","w");
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
	if(res2[usr_k][usr_j][usr_i] > 0.0)
		fprintf(output, "%d %d %d %20.20f\n", usr_i, usr_j, usr_k, res2[usr_k][usr_j][usr_i]);
fclose(output);

exit(0);

	
// min and max of the res2 and res3 - I need them different.
maxF2 = -100000.0; minF2 = 10000000.0;
totalRes2 = 0.0; countRes2 = 0;
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
if(res2[usr_k][usr_j][usr_i]>0.0){
if(minF2 > res2[usr_k][usr_j][usr_i] )  minF2 = res2[usr_k][usr_j][usr_i];
if(maxF2 < res2[usr_k][usr_j][usr_i] )  maxF2 = res2[usr_k][usr_j][usr_i];
countRes2++;
totalRes2 = totalRes2 + res2[usr_k][usr_j][usr_i];
}

meanRes2 = totalRes2 / (double)countRes2;

maxF3 = -100000.0; minF3 = 10000000.0;
totalRes3 = 0.0; countRes3 = 0;
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
if(res3[usr_k][usr_j][usr_i]>0.0){
if(minF3 > res3[usr_k][usr_j][usr_i] )  minF3 = res3[usr_k][usr_j][usr_i];
if(maxF3 < res3[usr_k][usr_j][usr_i] )  maxF3 = res3[usr_k][usr_j][usr_i];
countRes3++;
totalRes3 = totalRes3 + res3[usr_k][usr_j][usr_i];
}

meanRes3 = totalRes3 / (double)countRes3;

printf("at res2, min and max: %20.20f %20.20f %d. total: %f, mean: %f\n", minF2, maxF2, countRes2, totalRes2, meanRes2);
printf("at res3, min and max: %20.20f %20.20f %d. total: %f, mean: %f\n", minF3, maxF3, countRes3, totalRes3, meanRes3);


for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
fres2[usr_k][usr_j][usr_i] = res2[usr_k][usr_j][usr_i]  / meanRes2;

for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
fres3[usr_k][usr_j][usr_i] = res3[usr_k][usr_j][usr_i]  / meanRes3;

maxF2 = -10000000;
minF2 = 100000000000;
total_flow = 0; 
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
	if(fres2[usr_k][usr_j][usr_i] > 0.0){
			if(fres2[usr_k][usr_j][usr_i]>maxF2) maxF2 = fres2[usr_k][usr_j][usr_i];
			if(fres2[usr_k][usr_j][usr_i]<minF2) minF2 = fres2[usr_k][usr_j][usr_i];			
			total_flow = total_flow + fres2[usr_k][usr_j][usr_i];
		}

		mean = total_flow/(double)countRes2; // check that this is 1.

		if(maxF2 > 0.0 && minF2 > 0.0){
		deltaD_res2 = (maxF2 - minF2)/(double)(histNum+1);
		}
		else{
		printf("the range has a zero in it, exiting.\n"); exit(0);
		}

		printf("dx = 0.1. min of fres2: %f max: %f, mean: %f (the mean must be 1). %f\n", minF2, maxF2, mean, deltaD_res2);
		
maxF3 = -10000000;
minF3 = 100000000000;
total_flow = 0; 
		
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
	if(fres3[usr_k][usr_j][usr_i] > 0.0){
			if(fres3[usr_k][usr_j][usr_i]>maxF3) maxF3 = fres3[usr_k][usr_j][usr_i];
			if(fres3[usr_k][usr_j][usr_i]<minF3) minF3 = fres3[usr_k][usr_j][usr_i];			
			total_flow = total_flow + fres3[usr_k][usr_j][usr_i];
		}

		mean = total_flow/(double)countRes3; // check that this is 1.

		if(maxF3 > 0.0 && minF3 > 0.0){
		deltaD_res3 = (maxF3 - minF3)/(double)(histNum+1);
		}
		else{
		printf("the range has a zero in it, exiting.\n"); exit(0);
		}
		
		printf("dx = 0.5. min of fres2: %f max: %f, mean: %f (the mean must be 1). %f\n", minF3, maxF3, mean, deltaD_res3);		
		
		exit(0);
		
// now you have the deltaD, so you can do the resolution 1 histogram.		
for(temptemp=0; temptemp<histNum;temptemp++){
			hist[0][temptemp] = 0.0;			
			hist[1][temptemp] = 0.0;	
}

for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++){	
	hist[0][myfi_counter] = minF2 + deltaD_res2 * ((double)(myfi_counter));	
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
		if(fres2[usr_k][usr_j][usr_i] > 0.0){
//		if( fres2[usr_k][usr_j][usr_i] > (hist[0][myfi_counter] - deltaD_res2/2.0) && fres2[usr_k][usr_j][usr_i] < (hist[0][myfi_counter] + deltaD_res2/2.0) )
		if( (fres2[usr_k][usr_j][usr_i] <= hist[0][myfi_counter]) && (fres2[usr_k][usr_j][usr_i] > hist[0][myfi_counter]-deltaD_res2 * ((double)(myfi_counter))) )
			hist[1][myfi_counter] = hist[1][myfi_counter] + 1.0;
		}
} // end of myfi_counter.

	for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++) hist[1][myfi_counter] = hist[1][myfi_counter]/(deltaD_res2 * (double)countRes2);

//	hist[1][0] = hist[1][histNum] = 0.0;

// write pdf at resolution 1.
		output = fopen("hist_resolution2.data","w");
		for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++)
		fprintf(output,"%f %f \n", hist[0][myfi_counter], hist[1][myfi_counter] );
		fclose(output);

// do RD here.		
total_counts = 0.0;
hist_mean = hist_total = 0.0;
for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++){
	total_counts  = total_counts + hist[1][myfi_counter];
	hist_total = hist_total + hist[0][myfi_counter] * hist[1][myfi_counter];	
}
hist_mean = hist_total / total_counts;	
var_total = 0.0;

for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++) var_total = var_total + hist[1][myfi_counter] * (hist[0][myfi_counter] - hist_mean) * (hist[0][myfi_counter] - hist_mean);
sd1 = sqrt( var_total / total_counts );

printf("total counts, f: %f %f %f SD at dx = 0.1: %f\n", total_counts, hist_total, hist_mean, sd1);

		for(myfi_counter=0;myfi_counter<= histNum; myfi_counter++){
				hist[0][myfi_counter] = 0;
				hist[1][myfi_counter] = 0;		
		}		
		
// ============================================================================================================
// now you have the deltaD, so you can do the resolution 2 histogram.
for(myfi_counter=1;myfi_counter<=histNum; myfi_counter++){	
	hist[0][myfi_counter] = minF3 + deltaD_res3 * ((double)(myfi_counter));	
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
		if(fres3[usr_k][usr_j][usr_i] > 0.0){
//		if( fres3[usr_k][usr_j][usr_i] > (hist[0][myfi_counter] - deltaD_res3/2.0) && fres3[usr_k][usr_j][usr_i] < (hist[0][myfi_counter] + deltaD_res3/2.0) )
		if( (fres3[usr_k][usr_j][usr_i] <= hist[0][myfi_counter]) && (fres3[usr_k][usr_j][usr_i] > hist[0][myfi_counter]-((double)(myfi_counter)))  )
			hist[1][myfi_counter] = hist[1][myfi_counter] + 1.0;
		}
} // end of myfi_counter.

	for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++) hist[1][myfi_counter] = hist[1][myfi_counter]/(deltaD_res3 * (double)countRes3);

//	hist[1][0] = hist[1][histNum] = 0.0;

// write pdf at resolution 2
		output = fopen("hist_resolution3.data","w");
		for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++)
		fprintf(output,"%f %f \n", hist[0][myfi_counter], hist[1][myfi_counter] );
		fclose(output);

// do RD here.		
total_counts = 0.0;
hist_mean = hist_total = 0.0;
for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++){
	total_counts  = total_counts + hist[1][myfi_counter];
	hist_total = hist_total + hist[0][myfi_counter] * hist[1][myfi_counter];	
}
hist_mean = hist_total / total_counts;	

var_total = 0.0;
for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++) var_total = var_total + hist[1][myfi_counter] * (hist[0][myfi_counter] - hist_mean) * (hist[0][myfi_counter] - hist_mean);
sd6 = sqrt( var_total / total_counts );

printf("total counts for fres3, f: %f %f %f SD at dx = 0.5: %f\n", total_counts, hist_total, hist_mean, sd6);

		for(myfi_counter=0;myfi_counter<=histNum; myfi_counter++){
				hist[0][myfi_counter] = 0;
				hist[1][myfi_counter] = 0;		
		}		
		

return 0;

} // end of main.
