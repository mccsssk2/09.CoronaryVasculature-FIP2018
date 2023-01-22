/*
Oct 26 2017.
First, just calculate FD for control case using 100 instances, or all available instances.
Make a grid at res1 and res2 (3D arrays)
Total up perfusion at terminals in each box, take perfusion only from O(6) capilliaries.
make raw histogram
norm to 100 or however many you have.
normal to max
RD.
FD.


Reading the tree is expensive, so just read it and do the experiments on it: that is part 1 of post-processing.
Part 2 is the FD calculation.

16 Nov. 2017.
fileNum now may be a constant as I may wish to run expt (binary of this C program) in the directory.
*/

#define histNum 1000

#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.

int main(int argc,char **argv){
	
#include "declarations.c"	
int fileNum;
arr ***res1;

nodeInfo	   *plotThisSegment;
int histogram[50]; // the iterator is the perfusion value, the array is the frequency. I am only counting from 1 to 50.
double rawHist[2][101];
double resolution1, resolution2, resolution3, dx, dx2, dx3;

double res1PerfusionCountLV[12], res1PerfusionCountRV[6]; // each LV and RV are divided into layers of 1 mm each.

double ***res2, ***res3;
double ***fres1, ***fres2, ***fres3; // res1 is 1 mm, res2 is 0.1 mm, res3 is 0.2 mm.
int ***ires2, ***ires3;

double histres1[2][histNum+1], histres2[2][histNum+1], histres3[2][histNum+1]; // RD histogram.



double meanRes1, meanRes2, meanRes3, totalRes1, totalRes2, totalRes3, total_flow, meanres1, meanres2, meanres3;
int countRes1, countRes2, countRes3;
double deltaD_res1, deltaD_res2, deltaD_res3;

double maxF, minF;
double myfi, temp_f;
int temptemp;
double maxrawHist;
double total_counts, hist_mean, hist_total, var_total, sd1, sd6;
double Rlb, Q_perf_mod, fraction;

printf("In loops/parameter sweeps, the argv1 may or may not be used.\n");

if(argc<2){
printf("I need an argument, between 0 and 9, all different experiments need different arguments.\n");
 exit(0);
}

// I want the perfusion at a given resolution, resolution1.
dx = 1.0; // mm
mysize_x1=(int)((double)TOTAL_X/(double)dx);  mysize_y1=(int)((double)TOTAL_Y/(double)dx);   mysize_z1=(int)((double)TOTAL_Z/(double)dx);
printf("at dx1, The xmax, ymax, zmax are: %d %d %d\n", mysize_x1, mysize_y1, mysize_z1);
// declare your structs.
		res1												= (arr***)	calloc(mysize_z1,sizeof(arr **));
		for(usr_k = 0; usr_k < mysize_z1; usr_k++){
			res1[usr_k] 									= (arr**)  calloc(mysize_y1,sizeof(arr* ));
			for (usr_j = 0; usr_j < mysize_y1; usr_j++){
				res1[usr_k][usr_j] 							= (arr *) calloc(mysize_x1,sizeof(arr));
			}
		}

/*		
		fres1												= (double***)	calloc(mysize_z1,sizeof(double **));
		for(usr_k = 0; usr_k < mysize_z1; usr_k++){
			fres1[usr_k] 									= (double**)  calloc(mysize_y1,sizeof(double* ));
			for (usr_j = 0; usr_j < mysize_y1; usr_j++){
				fres1[usr_k][usr_j] 							= (double *) calloc(mysize_x1,sizeof(double));
			}
		}
*/				

// higher resolution.
dx2 = 1.0;
mysize_x2=(int)((double)TOTAL_X/(double)dx2);  mysize_y2=(int)((double)TOTAL_Y/(double)dx2); mysize_z2=(int)((double)TOTAL_Z/(double)dx2);
printf("a dx2, The xmax, ymax, zmax are: %d %d %d\n", mysize_x2, mysize_y2, mysize_z2);
// declare your structs.
		res2												= (double***)	calloc(mysize_z2,sizeof(double **));
//		fres2												= (double***)	calloc(mysize_z2,sizeof(double **));		
		ires2												= (int***)	calloc(mysize_z2,sizeof(int **));
		for(usr_k = 0; usr_k < mysize_z2; usr_k++){
			res2[usr_k] 									= (double**)  calloc(mysize_y2,sizeof(double* ));
//			fres2[usr_k] 									= (double**)  calloc(mysize_y2,sizeof(double* ));			
			ires2[usr_k] 										= (int**)  calloc(mysize_y2,sizeof(int* ));
			for (usr_j = 0; usr_j < mysize_y2; usr_j++){
				res2[usr_k][usr_j] 							= (double *) calloc(mysize_x2,sizeof(double));
//				fres2[usr_k][usr_j] 							= (double *) calloc(mysize_x2,sizeof(double));				
				ires2[usr_k][usr_j] 							= (int *) calloc(mysize_x2,sizeof(int));
			}
		}

// in the middle resolution.		
// dx3 = 0.5; 
dx3 = 2.0; // mm		
mysize_x3=(int)((double)TOTAL_X/(double)dx3);  mysize_y3=(int)((double)TOTAL_Y/(double)dx3); mysize_z3=(int)((double)TOTAL_Z/(double)dx3);
printf("at dx3, The xmax, ymax, zmax are: %d %d %d\n", mysize_x3, mysize_y3, mysize_z3);
// declare your structs.
		res3												= (double***)	calloc(mysize_z3, sizeof(double **));
//		fres3												= (double***)	calloc(mysize_z3, sizeof(double **));		
		ires3												= (int***)	calloc(mysize_z3, sizeof(int **));		
		for(usr_k = 0; usr_k < mysize_z3; usr_k++){
			res3[usr_k] 									= (double**)  calloc(mysize_y3, sizeof(double* ));
//			fres3[usr_k] 									= (double**)  calloc(mysize_y3, sizeof(double* ));			
			ires3[usr_k] 									= (int**)  calloc(mysize_y3, sizeof(int* ));
			for (usr_j = 0; usr_j < mysize_y3; usr_j++){
				res3[usr_k][usr_j] 							= (double *) calloc(mysize_x3, sizeof(double));
//				fres3[usr_k][usr_j] 							= (double *) calloc(mysize_x3, sizeof(double));				
				ires3[usr_k][usr_j] 							= (int *) calloc(mysize_x3, sizeof(int));
			}
		}

// ================== declarations/allocations finished. ===================================================================================================
// ================== declarations/allocations finished. ===================================================================================================

// read in the 1 mm distances from sub-dir: distances/DISTANCES64_85_74.DATUM
// res1 is only used for the epi endo heterogeniety diagram.
// use this distance only to set up the epiendo perfusion gradient.
output = fopen("/home/kharches/projects/McIntyre/SyntheticHearts/1.srkCodes/distances/DISTANCE64_85_74.DATUM", "r"); 
if(output==NULL){
	// graham.
	output = fopen("/home/kharches/scratch/distances/DISTANCE64_85_74.DATUM", "r"); 
		if(output==NULL){
		// orca.
		output = fopen("/work/kharches/distances/DISTANCE64_85_74.DATUM", "r"); 
		}
}

if(output==NULL){ printf("distances files!!\n"); exit(2); }

while( fscanf(output, "%d %d %d %lf\n", &o1, &o2, &o3, &dist)!=EOF) res1[o3][o2][o1].distanceToEndo = dist; fclose(output);

printf("finished reading in distances file.\n");

// =============================================================================================================
// Initialise arrays.
// =============================================================================================================
for(usr_k = 0; usr_k < 12; usr_k++) res1PerfusionCountLV[usr_k] = 0.0;
for(usr_k = 0; usr_k < 6; usr_k++)  res1PerfusionCountRV[usr_k] = 0.0;
for(usr_k = 0; usr_k < mysize_z1; usr_k++) for (usr_j = 0; usr_j < mysize_y1; usr_j++) for (usr_i = 0; usr_i < mysize_x1; usr_i++){ res1[usr_k][usr_j][usr_i].count = 0; res1[usr_k][usr_j][usr_i].total_perfusion = 0;		}
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++){ res2[usr_k][usr_j][usr_i] = 0; ires2[usr_k][usr_j][usr_i] = 0; }
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++){ res3[usr_k][usr_j][usr_i]= 0;  ires3[usr_k][usr_j][usr_i]= 0;  }
// =============================================================================================================
// =============================================================================================================
// =============================================================================================================

if(rootRCA!=NULL){ deleteTree(rootRCA); rootRCA = NULL; }   if(rootLA!=NULL){ deleteTree(rootLA); rootLA = NULL; }
// for(usr_k = 0; usr_k < 12; usr_k++) res1PerfusionCountLV[usr_k] = 0.0;   for(usr_k = 0; usr_k < 6; usr_k++)  res1PerfusionCountRV[usr_k] = 0.0;

// read in the inorder/preorder. do RCA first.
// sprintf(str,"ALL_%d/inorderRCA.txt", fileNum);
sprintf(str,"inorderRCA.txt");
output = fopen(str, "r"); 
if(output==NULL){ 
	fclose(output); exit(0); 
} 
else
{
	// count	
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF) count_value++; fclose(output);
	in_inorderRCA = (int *) calloc(count_value, sizeof(int)) ;
	// populate.
	output = fopen(str, "r"); 
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF){ in_inorderRCA[count_value] = intp1; count_value++; } fclose(output); count_max = count_value; 
}
if(debug) printf("In order says RCA has %d nodes.\n", count_max);

// sprintf(str,"ALL_%d/preorderRCA.txt", fileNum);
sprintf(str,"preorderRCA.txt");
output = fopen(str, "r");
if(output==NULL){ 
	fclose(output); exit(0); 
} 
else
{
	// count	
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF) count_value++; fclose(output);
	in_preorderRCA = (int *) calloc(count_value, sizeof(int)) ;
	// populate.
	output = fopen(str, "r");
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF){ in_preorderRCA[count_value] = intp1; count_value++; } fclose(output);
}
printf("Pre order says RCA has %d nodes, and in order said: %d (two integers must be the same)\n", count_value, count_max);

// construct the RCA tree.
preIndex = 0; rootRCA = buildTreeR(in_inorderRCA, in_preorderRCA, 0, count_value-1, 1);
constructTree(rootRCA, 1); 


 printf("I created these many nodes in RCA: %d\n", binarytree_count(rootRCA) );
// ===============================================================================================================
// ===============================================================================================================
// now read in the LA tree.
// read in the inorder/preorder. do LA first.
// sprintf(str,"ALL_%d/inorderLA.txt", fileNum);
sprintf(str,"inorderLA.txt");
output = fopen(str, "r");
if(output==NULL){ 
	fclose(output); exit(0); 
} 
else
{
	// count	
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF) count_value++; fclose(output);
	in_inorderLA = (int *) calloc(count_value, sizeof(int)) ;
	// populate.
	output = fopen(str, "r");
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF){ in_inorderLA[count_value] = intp1; count_value++; } fclose(output); count_max = count_value;
}
if(debug)
printf("In order says LA has %d nodes.\n", count_max);

// sprintf(str,"ALL_%d/preorderLA.txt", fileNum);
sprintf(str,"preorderLA.txt");
output = fopen(str, "r");
if(output==NULL){ 
	fclose(output); exit(0); 
} 
else
{
	// count	
	count_value = 0; while(fscanf(output, "%d ", &intp2)!=EOF) count_value++; fclose(output);
	in_preorderLA = (int *) calloc(count_value, sizeof(int)) ;
	// populate.
	output = fopen(str, "r");
	count_value = 0; while(fscanf(output, "%d ", &intp1)!=EOF){ in_preorderLA[count_value] = intp1; count_value++; } fclose(output);
}
printf("Pre order says LA has %d nodes, and in order said: %d (two integers must be the same)\n", count_value, count_max);

// construct the LA tree.
preIndex = 0;  rootLA = buildTreeR(in_inorderLA, in_preorderLA, 0, count_value-1, 2);
constructTree(rootLA, 2);  
 
 printf("I created these many nodes in LA: %d\n", binarytree_count(rootLA) ); 

 // ===========================================================================================================================================
 // ===========================================================================================================================================

if(debug){
printf("The number of leaves in read RCA are: %d\n", getLeafCount(rootRCA));
printf("The number of leaves in read LA are: %d\n", getLeafCount(rootLA));
 
printf("The total perfusion in read RCA leaves is: %f\n",  totalPerfusionLeaves(rootRCA) );
printf("The total perfusion in read LA leaves is: %f\n",  totalPerfusionLeaves(rootLA) );
}
  
 // ===========================================================================================================================================
 // ===========================================================================================================================================
/* revise the tree:
1. All terminals are now to be O(6).
2. Recalculate perfusion, resistance, pressure, everything.
3. After 2, then do your experiments.
*/

// experiment loops start here.
// for(fileNum = 0; fileNum <= 10; fileNum++){
for(fileNum = atoi(argv[1]) ; fileNum <= atoi(argv[1]); fileNum++){

// =============================================================================================================
// Initialise arrays.
// =============================================================================================================
for(usr_k = 0; usr_k < 12; usr_k++) res1PerfusionCountLV[usr_k] = 0.0;
for(usr_k = 0; usr_k < 6; usr_k++)  res1PerfusionCountRV[usr_k] = 0.0;
for(usr_k = 0; usr_k < mysize_z1; usr_k++) for (usr_j = 0; usr_j < mysize_y1; usr_j++) for (usr_i = 0; usr_i < mysize_x1; usr_i++){ res1[usr_k][usr_j][usr_i].count = 0; res1[usr_k][usr_j][usr_i].total_perfusion = 0;		}
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++){ res2[usr_k][usr_j][usr_i] = 0; ires2[usr_k][usr_j][usr_i] = 0; }
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++){ res3[usr_k][usr_j][usr_i]= 0;  ires3[usr_k][usr_j][usr_i]= 0;  }
// =============================================================================================================
	
double Pinin;
	
// make all terminals O(6) and same radius, length.
reviseTerminals(rootRCA);
reviseTerminals(rootLA);

Pinin = Pin; // control, 100 mmHg.

// expt 1.
//	Pinin = Pin *(2.0 - (double)fileNum / 10.0); Pressure affects perfusion not hetrogeniety run.

// expt. 2.
// if(fileNum>0) stenoseThisArtery(rootLA, fileNum); 

// expt 3. order by order stenosis. The first segment of given order is stenosed to radius 0.04 mm, which is less than most.
// fileNum is the order which goes from 5 to 11. 5 is control.
// stenoseAllArteriesInThisOrder(rootLA, fileNum);    stenoseAllArteriesInThisOrder(rootRCA, fileNum); 

// expt 4: anisotropy ratio. This didnt work yet. I just get very small BF at the ends.
// IncreaseAnisotropyRatio(rootLA, fileNum); IncreaseAnisotropyRatio(rootRCA, fileNum); 

// expt 5: range of pressure: 30, 100, 200 mmHg AND range of stenosis at O(6): 0.1, 0.3, 0.5.
//	if( atoi(argv[2])==1 ) Pinin = 30.0 * 133.33;
// else if( atoi(argv[2])==2 ) Pinin = 100.0 * 133.33;
// else if( atoi(argv[2])==3 ) Pinin = 200.0 * 133.33;
// else if( atoi(argv[2])<1 || atoi(argv[2])>3 ){ printf("no such case in expt 5, exiting.\n"); exit(0); }

// expt 6.
// Pinin = 30.0 * 133.33; // run1.
// Pinin = 100.0 * 133.33; // run2
// Pinin = 200.0 * 133.33; // run3
// stenoseAllArteriesInThisOrder2(rootLA, fileNum);    stenoseAllArteriesInThisOrder2(rootRCA, fileNum);  // fileNum is atoi(  argv[1]  ) this one pinches all arteries at SN 6.

// expt 7.
// Pinin = 30.0 * 133.33; // run1.
// Pinin = 100.0 * 133.33; // run2 FINISHED
Pinin = 200.0 * 133.33; // run3
stenoseAllArteriesInThisOrder3(rootLA, fileNum);    stenoseAllArteriesInThisOrder3(rootRCA, fileNum); 

	totalTreeResistanceRCA  = binaryDownResistance3D(rootRCA); totalTreeResistanceLA     = binaryDownResistance3D(rootLA);
	rootRCA->pressure 	= Pinin; rootLA->pressure = Pinin; // Units: kg/(m-s2)

		if(totalTreeResistanceRCA + totalTreeResistanceLA > 0.0)
		totalTreeResistance = totalTreeResistanceRCA * totalTreeResistanceLA / (totalTreeResistanceRCA + totalTreeResistanceLA);
		else
		printf("your LA or RCA or both resistances are nada, check!");

// I have to calculate at least one: q_perf, Pout,or totalTreeResistance. I choose Pout now that all terminals are O(6).		
		Pout = 20.0 * 133.33; // make sure that your terminal pressure is never negative.
		Q_perf_mod = (Pinin - Pout)/totalTreeResistance;
rootRCA->perfusion = (Pinin - Pout)/ totalTreeResistanceRCA;  rootLA->perfusion = (Pinin - Pout) / totalTreeResistanceLA ;
rootRCA->theta 	   = 0.0; 												 rootLA->theta  	 = 0.0;
printf("Revised tree RCA and LA: %f %f %f %f %f %f %f\n", Q_perf, Q_perf_mod, rootRCA->perfusion, rootLA->perfusion,  rootRCA->perfusion + rootLA->perfusion, Pinin, Pout );
printf("In main, resistances: %f %f %f\n", rootRCA->downResistance, rootLA->downResistance , totalTreeResistance );

perfusionPressure3D(rootRCA, Pout);   									perfusionPressure3D(rootLA, Pout);

 // verified that the rca and LA are read in properly by 1) making raw files; 2) comparing the VTK outputs.
 // ===========================================================================================================================================

// populate the res1 and res2 3D arrays, as well as the epiendo array.
writeAtLeaf(rootRCA, fileNum, res1, dx, res2, dx2, res3, dx3, res1PerfusionCountRV, res1PerfusionCountLV, ires2, ires3); 
writeAtLeaf(rootLA,    fileNum, res1, dx, res2, dx2, res3, dx3, res1PerfusionCountRV, res1PerfusionCountLV, ires2, ires3);

// write LV and RV epi-endo.
sprintf(str,"epiendoRV%d.data", fileNum ); output = fopen(str, "w"); for(o1=0; o1<6;o1++)  fprintf(output, "%f %f\n", (double)o1, res1PerfusionCountRV[o1] ); fclose(output);
sprintf(str,"epiendoLV%d.data", fileNum ); output = fopen(str, "w"); for(o1=0; o1<12;o1++) fprintf(output, "%f %f\n", (double)o1, res1PerfusionCountLV[o1] ); fclose(output);


// these are already sums of all the instances that you are running over.
// write res2. 
sprintf(str, "order6res2%d.data", fileNum );
output = fopen(str,"w");
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
if(res2[usr_k][usr_j][usr_i]>0.0)
fprintf(output, "%d %d %d %20.20f\n", usr_i, usr_j, usr_k, res2[usr_k][usr_j][usr_i]);
fclose(output);

sprintf(str, "Iorder6res2%d.data", fileNum );
output = fopen(str,"w");
for(usr_k = 0; usr_k < mysize_z2; usr_k++) for (usr_j = 0; usr_j < mysize_y2; usr_j++) for (usr_i = 0; usr_i < mysize_x2; usr_i++)
if(ires2[usr_k][usr_j][usr_i]>0)
fprintf(output, "%d %d %d %d\n", usr_i, usr_j, usr_k, ires2[usr_k][usr_j][usr_i]);
fclose(output);

// write res3.
sprintf(str, "order6res3%d.data" , fileNum );
output = fopen(str,"w");
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
//	if(res3[usr_k][usr_j][usr_i]>0.001)
if(res3[usr_k][usr_j][usr_i]>0.0)	
fprintf(output, "%d %d %d %20.20f\n", usr_i, usr_j, usr_k, res3[usr_k][usr_j][usr_i]);
fclose(output);

sprintf(str, "Iorder6res3%d.data" , fileNum );
output = fopen(str,"w");
for(usr_k = 0; usr_k < mysize_z3; usr_k++) for (usr_j = 0; usr_j < mysize_y3; usr_j++) for (usr_i = 0; usr_i < mysize_x3; usr_i++)
if(ires3[usr_k][usr_j][usr_i]>0.0)	
fprintf(output, "%d %d %d %d\n", usr_i, usr_j, usr_k, ires3[usr_k][usr_j][usr_i]);
fclose(output);

//============================================================================================================================================================

writeTreeSoFar2(rootRCA, fileNum);  writeTreeSoFar2(rootLA, fileNum); // this fileNum may be redandant now.
#include "writeVTKS1.c"

// measure anisotropy ratio.
writeAnisotropicRatios(rootRCA, fileNum);   writeAnisotropicRatios(rootLA, fileNum);

// restore data. used in simulations for order by order stenosis, and anisotropy.
// restoreTreeData(rootRCA); restoreTreeData(rootLA);

} // end of fileNum loop.

return 0;

} // end of main.
