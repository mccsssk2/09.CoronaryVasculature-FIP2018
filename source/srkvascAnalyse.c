/*
16 Nov 2017.
I worked out how to run the expt jobs in serial. The jobs by definition make output for each of the 541 directory.
This program takes the order20 from each directory, and make the composite.
*/

#define divs 100

#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.

int main(int argc,char **argv){
	
#include "declarations.c"	
int fileNum;

double 	***res1, ***res2, ***res3;
int 		***ires1, ***ires2, ***ires3;

double ***fres1, ***fres2, ***fres3; // res1 is 1 mm, res2 is 0.1 mm, res3 is 0.2 mm.

double resolution1, resolution2, resolution3, dx, dx2, dx3;
double total_flow, variance, std_dev, meanT, fisquared , fisquared1, numberofVoxelDistant, distto, fisquaredn, tempsumn;
double deltaD;
int total_pixels;
double maxF, minF;
double hist[2][divs+1]; // this histogram is for the masses.
double myfi, temp_f;
int myfi_counter;
double total_counts, hist_mean, hist_total, var_total, sd1, sd6;
double norm;
total_counts = 0.0;
hist_mean = hist_total = 0.0;
double area;

// double gammaVal[12];
int gammaCount[12];
double gammaSTD[12];
double gammaMean[12];

double q00, qL, qR, q_bigger, q_smaller;
double gammaL, gammaR, gamma_big, gamma_small;
int SNp, SNL, SNR;
double prev_mean;


if(argc<3){
printf("I want 3 arguments, start and finish, as well as what case you want FD for.\n");
exit(0);
}


for(usr_i=0;usr_i<12; usr_i++){
//	gammaVal[usr_i] = 0.0;
	gammaCount[usr_i] = 0;
	gammaSTD[usr_i] = 0.0;
	gammaMean[usr_i] = 0.0;	
}

// higher resolution.
dx2 = 1.0;
mysize_x2=(int)((double)TOTAL_X/(double)dx2);  mysize_y2=(int)((double)TOTAL_Y/(double)dx2); mysize_z2=(int)((double)TOTAL_Z/(double)dx2);
printf("at dx2, The xmax, ymax, zmax are: %d %d %d\n", mysize_x2, mysize_y2, mysize_z2);
// declare your structs.
		res2												= (double***)	calloc(mysize_z2,sizeof(double **));
		fres2												= (double***)	calloc(mysize_z2,sizeof(double **));		
		ires2												= (int***)	calloc(mysize_z2,sizeof(int **));
		for(usr_k = 0; usr_k < mysize_z2; usr_k++){
			res2[usr_k] 										= (double**)  calloc(mysize_y2,sizeof(double* ));
			fres2[usr_k] 										= (double**)  calloc(mysize_y2,sizeof(double* ));
			ires2[usr_k] 										= (int**)  calloc(mysize_y2,sizeof(int* ));
			for (usr_j = 0; usr_j < mysize_y2; usr_j++){
				res2[usr_k][usr_j] 							= (double *) calloc(mysize_x2,sizeof(double));
				fres2[usr_k][usr_j] 							= (double *) calloc(mysize_x2,sizeof(double));
				ires2[usr_k][usr_j] 							= (int *) calloc(mysize_x2,sizeof(int));				
			}
		}
// initialise, just make sure.
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++){ res2[usr_k][usr_j][usr_i] = fres2[usr_k][usr_j][usr_i] =  0; }

// dx3 arrays.
// in the middle resolution.		
dx3 = 2.0; // mm		
mysize_x3=(int)((double)TOTAL_X/(double)dx3);  mysize_y3=(int)((double)TOTAL_Y/(double)dx3); // mysize_z3=(int)((double)TOTAL_Z/(double)dx3);
mysize_z3 = mysize_z2;
printf("at dx3, The xmax, ymax, zmax are: %d %d %d\n", mysize_x3, mysize_y3, mysize_z3);
// declare your structs.
		res3												= (double***)	calloc(mysize_z3, sizeof(double **));
		fres3												= (double***)	calloc(mysize_z3, sizeof(double **));		
		ires3												= (int***)	calloc(mysize_z3, sizeof(int **));		
		for(usr_k = 0; usr_k < mysize_z3; usr_k++){
			res3[usr_k] 									= (double**)  calloc(mysize_y3, sizeof(double* ));
			fres3[usr_k] 									= (double**)  calloc(mysize_y3, sizeof(double* ));
			ires3[usr_k] 									= (int**)  calloc(mysize_y3, sizeof(int* ));
			for (usr_j = 0; usr_j < mysize_y3; usr_j++){
				res3[usr_k][usr_j] 							= (double *) calloc(mysize_x3, sizeof(double));
				fres3[usr_k][usr_j] 							= (double *) calloc(mysize_x3, sizeof(double));				
				ires3[usr_k][usr_j] 							= (int *) calloc(mysize_x3, sizeof(int));
			}
		}
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++){ res3[usr_k][usr_j][usr_i] = fres3[usr_k][usr_j][usr_i] = 0;  }

// now just read it in.

start = atoi(argv[1]);
end = atoi(argv[2]);

norm = (double)(end - start);

for(fileNum=start; fileNum<=end; fileNum++){
//	printf("%d  .  ", fileNum);

sprintf(str,"ALL_%d/order6res2%d.data", fileNum, (int)atoi(argv[3]) );
output = fopen(str,"r");
if(output==NULL){
	printf("order6res in %d did not open/doesnt exist. exiting.\n", fileNum);
	exit(0);
}
while(fscanf(output,"%d %d %d %lf", &usr_i, &usr_j, &usr_k, &viz_perfusion)!=EOF){
if(debug)	
printf("%d %d %d %lf\n", usr_i, usr_j, usr_k, viz_perfusion);
	// assign if sensible.
	if( usr_i >= 0 && usr_i < mysize_x2 && usr_j >= 0 && usr_j < mysize_y2 && usr_k >= 0 && usr_k < mysize_z2 ){
		res2[usr_k][usr_j][usr_i] = res2[usr_k][usr_j][usr_i] + viz_perfusion;
	}
}
fclose(output);

sprintf(str,"ALL_%d/Iorder6res2%d.data", fileNum, (int)atoi(argv[3]) );
output = fopen(str,"r");
if(output==NULL){
	printf("order6res in %d did not open/doesnt exist. exiting.\n", fileNum);
	exit(0);
}
while(fscanf(output,"%d %d %d %d", &usr_i, &usr_j, &usr_k, &numberOfDo)!=EOF){
if(debug)	
printf("%d %d %d %lf\n", usr_i, usr_j, usr_k, viz_perfusion);
	// assign if sensible.
	if( usr_i >= 0 && usr_i < mysize_x2 && usr_j >= 0 && usr_j < mysize_y2 && usr_k >= 0 && usr_k < mysize_z2 ){
		ires2[usr_k][usr_j][usr_i] = ires2[usr_k][usr_j][usr_i] + numberOfDo;
	}
}
fclose(output);

sprintf(str,"ALL_%d/anisotropy%d.data", fileNum, (int)atoi(argv[3]) );
output = fopen(str,"r");
if(output==NULL){
	printf("anisotropy data in %d did not open/doesnt exist. exiting.\n", fileNum);
	exit(0);
}
while(  fscanf(output,"%d %d %d %lf %lf %lf %lf %lf %lf %lf\n", &SNp, &SNL, &SNR, &q0, &qL, &qR, &gammaL, &gammaR, &gamma_big, &gamma_small) !=EOF  ){	
	if(SNp > 5 && SNp < 12){
		if(SNp > 0 && SNL > 0 && SNR > 0 && gammaL > 0.0){
		gammaCount[SNp]++;
		prev_mean = gammaMean[SNp];
		if(gammaCount[SNp] > 0)
		gammaMean[SNp] = gammaMean[SNp] + (gammaL - gammaMean[SNp]) / (double)gammaCount[SNp];
		gammaSTD[SNp] = gammaSTD[SNp] + (gammaL - gammaMean[SNp]) * (gammaL - prev_mean);
		}
	}
}
fclose(output);


/*
sprintf(str,"ALL_%d/order6res30.data", fileNum);
output = fopen(str,"r");
if(output==NULL){
	printf("order6res in %d did not open/doesnt exist. exiting.\n", fileNum);
	exit(0);
}
while(fscanf(output,"%d %d %d %lf", &usr_i, &usr_j, &usr_k, &viz_perfusion)!=EOF){
if(debug)	
printf("%d %d %d %lf\n", usr_i, usr_j, usr_k, viz_perfusion);
	// assign if sensible.
	if( usr_i >= 0 && usr_i < mysize_x3 && usr_j >= 0 && usr_j < mysize_y3 && usr_k >= 0 && usr_k < mysize_z3 ){
		res3[usr_k][usr_j][usr_i] = res3[usr_k][usr_j][usr_i] + viz_perfusion;
	}
}
fclose(output);
*/

} // end of fileNum loop.


for(usr_k = 0; usr_k < mysize_z2; usr_k++)
	for(usr_j = 0; usr_j < mysize_y2; usr_j++)
		for(usr_i = 0; usr_i < mysize_x2; usr_i++){
			if(ires2[usr_k][usr_j][usr_i]> 0 )
				res2[usr_k][usr_j][usr_i]  = res2[usr_k][usr_j][usr_i] / ( (double)ires2[usr_k][usr_j][usr_i] );
		}


// once read in, do the histogram and do the FD for res3 whatever res3 you decide to be.
total_flow = 0.0;
total_pixels = 0;
maxF = -10000000;
minF = 100000000000;
for(usr_k = 0; usr_k < mysize_z2; usr_k++)
	for(usr_j = 0; usr_j < mysize_y2; usr_j++)
		for(usr_i = 0; usr_i < mysize_x2; usr_i++){
			if(res2[usr_k][usr_j][usr_i] > 0){
				total_pixels++;
				total_flow = total_flow + res2[usr_k][usr_j][usr_i];
				if(res2[usr_k][usr_j][usr_i]>maxF) maxF = res2[usr_k][usr_j][usr_i];
				if(res2[usr_k][usr_j][usr_i]<minF) minF = res2[usr_k][usr_j][usr_i];							
			}
		} // end of for loops.

deltaD = (maxF - minF)/(double)divs;

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = hist[1][myfi_counter] = 0.0;
}


for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = minF + deltaD * ((double)(myfi_counter));	
for(usr_k = 0; usr_k < mysize_z2; usr_k++)
	for(usr_j = 0; usr_j < mysize_y2; usr_j++)
		for(usr_i = 0; usr_i < mysize_x2; usr_i++)
		if(res2[usr_k][usr_j][usr_i] > 0.0){
		if(res2[usr_k][usr_j][usr_i] > (hist[0][myfi_counter]) && res2[usr_k][usr_j][usr_i] < (hist[0][myfi_counter] + deltaD ) )
			hist[1][myfi_counter] = hist[1][myfi_counter] + 1.0;
		}
} // end of myfi_counter.

// write it.
hist[1][0] = hist[1][divs] = 0.0;
// write raw histogram,
		sprintf(str,"hist_resolutionRAW%d.data", (int)atoi(argv[3]) );
		output = fopen(str, "w");
		for(myfi_counter=0;myfi_counter<=divs; myfi_counter++)
		fprintf(output,"%f %f \n", hist[0][myfi_counter], hist[1][myfi_counter] );
		fclose(output);

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = hist[1][myfi_counter] = 0.0;
}

meanT = total_flow / (double)total_pixels;

printf("\n");
printf("raw flow. total flow: %f; mean: %f; live pixels: %d; total pixels at this resolution: %d\n", total_flow, meanT, total_pixels, mysize_x2 * mysize_y2 * mysize_z2);

total_flow = 0.0;
maxF = -10000000;
minF = 100000000000;

for(usr_k = 0; usr_k < mysize_z2; usr_k++)
	for(usr_j = 0; usr_j < mysize_y2; usr_j++)
		for(usr_i = 0; usr_i < mysize_x2; usr_i++){
			if(res2[usr_k][usr_j][usr_i] > 0){
				 fres2[usr_k][usr_j][usr_i] = res2[usr_k][usr_j][usr_i]/meanT;
				 total_flow = total_flow + fres2[usr_k][usr_j][usr_i]; 
				if(fres2[usr_k][usr_j][usr_i]>maxF) maxF = fres2[usr_k][usr_j][usr_i];
				if(fres2[usr_k][usr_j][usr_i]<minF) minF = fres2[usr_k][usr_j][usr_i];			
			}
		} // end of for loops.
meanT = total_flow / (double)total_pixels;
deltaD = (maxF - minF)/(double)divs;
printf("f flow. total flow: %f. mean (has to be 1): %f. min F: %20.20f. max F: %f. deltaD: %f.\n", total_flow, meanT, minF, maxF, deltaD); // this is okay.

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = hist[1][myfi_counter] = 0.0;
}

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = minF + deltaD * ((double)(myfi_counter));	
for(usr_k = 0; usr_k < mysize_z2; usr_k++)
	for(usr_j = 0; usr_j < mysize_y2; usr_j++)
		for(usr_i = 0; usr_i < mysize_x2; usr_i++)
		if(fres2[usr_k][usr_j][usr_i] > 0.0){
		if( fres2[usr_k][usr_j][usr_i] > (hist[0][myfi_counter]) && fres2[usr_k][usr_j][usr_i] < (hist[0][myfi_counter] + deltaD ) )
			hist[1][myfi_counter] = hist[1][myfi_counter] + 1.0;
		}
} // end of myfi_counter.

printf("\n");

// write it.
hist[1][0] = hist[1][divs] = 0.0;

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
			hist[1][myfi_counter] = hist[1][myfi_counter]/(deltaD * (double)total_pixels);	
}

// write pdf at resolution 1.
		sprintf(str,"hist_resolution1%d.data", (int)atoi(argv[3]) );
		output = fopen(str,"w");
		for(myfi_counter=0;myfi_counter<=divs; myfi_counter++)
		fprintf(output,"%f %f \n", hist[0][myfi_counter], hist[1][myfi_counter] );
		fclose(output);

total_counts = 0.0;
hist_mean = hist_total = 0.0;
for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){
	total_counts  = total_counts + hist[1][myfi_counter];
	hist_total = hist_total + hist[0][myfi_counter] * hist[1][myfi_counter];	// this multipled by deltaD has to work out 1.
}

hist_mean = hist_total / total_counts;	

var_total = 0.0;
for(myfi_counter=0;myfi_counter<=divs; myfi_counter++) var_total = var_total + hist[1][myfi_counter] * (hist[0][myfi_counter] - hist_mean) * (hist[0][myfi_counter] - hist_mean);
sd1 = sqrt( var_total / total_counts );

printf("total counts, f. total counts: %f. area approx. of histogram (must be 1): %f. mean (had to be 1): %f. sd1: %f\n", total_counts, deltaD * hist_total, hist_mean, sd1);
	
// you need to use hist again for the next resolution, so reset to 0.		
		for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){
				hist[0][myfi_counter] = 0;
				hist[1][myfi_counter] = 0;		
		}		
//=================================================================================================================
//=================================dx3================================================================================

 for(usr_k = 0; usr_k < mysize_z2-2; usr_k+=2)
	for(usr_j = 0; usr_j < mysize_y2-2; usr_j+=2)
	for(usr_i = 0; usr_i < mysize_x2-2; usr_i+=2){
 if(
 res2[usr_k][usr_j][usr_i]>0.0 &&res2[usr_k][usr_j][usr_i+1]>0.0&&res2[usr_k][usr_j+1][usr_i]>0.0&&res2[usr_k][usr_j+1][usr_i+1]>0.0 &&res2[usr_k+1][usr_j] [usr_i]>0.0&&res2[usr_k+1][usr_j][usr_i+1]>0.0&&res2[usr_k+1][usr_j+1][usr_i]>0.0&&res2[usr_k+1][usr_j+1][usr_i+1]>0.0  ){
			temp_f = (res2[usr_k][usr_j][usr_i] + res2[usr_k][usr_j][usr_i+1] + res2[usr_k][usr_j+1][usr_i] + res2[usr_k][usr_j+1][usr_i+1]  + 
					res2[usr_k+1][usr_j][usr_i]  + res2[usr_k+1][usr_j][usr_i+1] + res2[usr_k+1][usr_j+1][usr_i] + res2[usr_k+1][usr_j+1][usr_i+1])/8.0;
					if(temp_f > 0.0)
					res3[usr_k][usr_j/2][usr_i/2] = temp_f;
			}
		}

/*
for(usr_k = 0; usr_k < mysize_z3; usr_k++)
	for(usr_j = 0; usr_j < mysize_y3; usr_j++)
		for(usr_i = 0; usr_i < mysize_x3; usr_i++){
			if(res3[usr_k][usr_j][usr_i] > 10.0  || res3[usr_k][usr_j][usr_i] < 0.1 ){
				res3[usr_k][usr_j][usr_i] = 0.0;
			}
		}
*/

		
total_pixels = 0; total_flow = 0.0;
for(usr_k = 0; usr_k < mysize_z3; usr_k++)
	for(usr_j = 0; usr_j < mysize_y3; usr_j++)
		for(usr_i = 0; usr_i < mysize_x3; usr_i++){
			if(res3[usr_k][usr_j][usr_i] > 0.0){
				total_pixels++;
				total_flow = total_flow + res3[usr_k][usr_j][usr_i];
			}
		}

meanT = total_flow / (double)total_pixels;
printf("The mean of the raw flow at res3: %f mean: %f %d %d\n", total_flow, meanT, total_pixels, mysize_x3 * mysize_y3 * mysize_z3);
total_flow = 0.0;

/* maxF = -10000000; minF = 100000000000; */

for(usr_k = 0; usr_k < mysize_z3; usr_k++)
	for(usr_j = 0; usr_j < mysize_y3; usr_j++)
		for(usr_i = 0; usr_i < mysize_x3; usr_i++){
				fres3[usr_k][usr_j][usr_i] = res3[usr_k][usr_j][usr_i]/(meanT);
				 total_flow = total_flow + fres3[usr_k][usr_j][usr_i]; 
				 if(fres3[usr_k][usr_j][usr_i] > 0.0){
				if(fres3[usr_k][usr_j][usr_i]>maxF) maxF = fres3[usr_k][usr_j][usr_i];
				if(fres3[usr_k][usr_j][usr_i]<minF) minF = fres3[usr_k][usr_j][usr_i];							
				 }
			}

meanT = total_flow /( (double)total_pixels );
deltaD = (maxF - minF)/(double)divs;
printf("The mean of the fres3 flow: %f mean (has to be 1): %f %d %d %f %f %f\n", total_flow, meanT, total_pixels, mysize_x3 * mysize_y3 * mysize_z3, minF, maxF, deltaD);

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = hist[1][myfi_counter] = 0.0;
}

for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
	hist[0][myfi_counter] = minF + deltaD * ((double)(myfi_counter));	
for(usr_k = 0; usr_k < mysize_z3; usr_k++)
	for(usr_j = 0; usr_j < mysize_y3; usr_j++)
		for(usr_i = 0; usr_i < mysize_x3; usr_i++)
		if(fres3[usr_k][usr_j][usr_i] > 0.0){
		if( fres3[usr_k][usr_j][usr_i] > (hist[0][myfi_counter]) && fres3[usr_k][usr_j][usr_i] < (hist[0][myfi_counter] + deltaD ) )
			hist[1][myfi_counter] = hist[1][myfi_counter] + 1.0;
		}
} // end of myfi_counter.

// write it.
// hst[1][0] = hist[1][divs] = 0.0;
		
for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){	
			hist[1][myfi_counter] = hist[1][myfi_counter]/(deltaD * (double)total_pixels);	
}

// write pdf at resolution 2.
		sprintf(str,"hist_resolution3%d.data", (int)atoi(argv[3]) );
		output = fopen(str,"w");
		for(myfi_counter=0;myfi_counter<=divs; myfi_counter++)
		fprintf(output,"%f %f \n", hist[0][myfi_counter], hist[1][myfi_counter] );
		fclose(output);

total_counts = 0.0;
hist_mean = hist_total = 0.0;
for(myfi_counter=0;myfi_counter<=divs; myfi_counter++){
	total_counts  = total_counts + hist[1][myfi_counter];
	hist_total = hist_total + hist[0][myfi_counter] * hist[1][myfi_counter];	
}

hist_mean = hist_total / total_counts;	

var_total = 0.0;
for(myfi_counter=0;myfi_counter<=divs; myfi_counter++) var_total = var_total + hist[1][myfi_counter] * (hist[0][myfi_counter] - hist_mean) * (hist[0][myfi_counter] - hist_mean);
sd6 = sqrt( var_total / total_counts );

printf("total counts, f: %f area approx. (has to be 1): %f mean approx. (is it 1?): %f sd3: %f\n", total_counts, deltaD * hist_total, hist_mean, sd6);

	double fd;
	
	fd  = 1.0-log(sd6/sd1)/log(8.0); // 1.0231

	printf("my fractal dimension is %f\n", fd);
	sprintf(str, "fractalRD%d.data", (int)atoi(argv[3]) );
	output = fopen(str,"w");
	fprintf(output, "%d %f %f %f\n", (int)atoi(argv[3]),  fd, sd1, sd6);
	fclose(output);

// my gamma statistics.
sprintf(str, "gammaStats%d.data", (int)atoi(argv[3]) );
output = fopen(str,"w");
for(usr_i=6;usr_i<12; usr_i++){
	fprintf(output, "%d %20.20f %20.20f\n", usr_i, gammaMean[usr_i] , sqrt( gammaSTD[usr_i] / (double) gammaCount[usr_i] ) );
}
fclose(output);
	
	

// see: see the low res, never mind the high res.
// if(debug){
		sprintf(str,"perfRes2_%d.vtk", (int)atoi(argv[3]));
		output = fopen(str,"w");
		fprintf(output,"# vtk DataFile Version 1.0\n");
		fprintf(output,"vtk input\n");
		fprintf(output,"ASCII\n");
		fprintf(output,"DATASET STRUCTURED_POINTS\n");
		fprintf(output,"DIMENSIONS %d %d %d\n",mysize_x2, mysize_y2, mysize_z2);
		fprintf(output,"SPACING 0.1 0.1 0.1\n");
		fprintf(output,"ORIGIN 0 0 0\n");
		fprintf(output,"POINT_DATA %d\n",mysize_x2 * mysize_y2 * mysize_z2);
		fprintf(output,"SCALARS Q float 1\n");
		fprintf(output,"LOOKUP_TABLE default\n");
			for(usr_k = 0; usr_k < mysize_z2; usr_k++)
			for(usr_j = 0; usr_j < mysize_y2; usr_j++){
			for(usr_i = 0; usr_i < mysize_x2; usr_i++)
				if(res2[usr_k][usr_j][usr_i] > 0.0)
				fprintf(output,"%f ",  res2[usr_k][usr_j][usr_i] );
				else
				fprintf(output,"%f ",  -10.0 );					
				fprintf(output,"\n");
			}
		fclose(output);
// }

return 0;
}

