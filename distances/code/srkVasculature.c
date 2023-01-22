/*
20 July. Complete revision of the 3D code. Fit the vasculature in the 3D volume given.
As much as possible, do not use any grid, rather just the equations of the volumes, e.g. ellipsoids.


Sanjay Kharche.
13 February 2017.

This program works out the vasculature from arteries till arteriole level in the 3D LHSC heart.
This code is based on Keelan 2016 Phil Trans paper.

1 proc job.

26 Februay 2017.
_______________

What the program does so far:
1) generate the bi-ventricular geometry.
2) generate a binary tree with given number of nodes.
3) calculate tissue supply and volume cost functions.
4) Iterate over minimising these cost functions.
5) The iterations are by one part of the annealing: move 1 node randomly to see if the cost goes down.

The immediate next step would be:
Second part of the annealing, which is transplant sub-trees randomly.

Construct binary tree from arrays:
http://www.geeksforgeeks.org/construct-tree-from-given-inorder-and-preorder-traversal/
http://www.geeksforgeeks.org/construct-a-binary-tree-from-postorder-and-inorder/

http://stackoverflow.com/questions/37481526/preorder-tree-traversal-function-which-returns-an-array-of-integers-in-c
http://www.geeksforgeeks.org/inorder-tree-traversal-without-recursion/
http://www.geeksforgeeks.org/construct-tree-inorder-level-order-traversals/
http://www.ritambhara.in/storing-binary-tree-in-a-file/

2 March 2017.
Both, moving coordinates and rearranging segments (parts of the annealing) are now programmed.
The rearragment has some problems, see the comments below.

We are set for:
1) Power cost to pump blood, sec 2.2 of Keelan.
2) Pressure drop, sec 2.5
3) Exclusion of large blood vessals from intererior.


16 August 2017.

To do:
1) calculate resistance, pressure-perfusion (in that order).
2) big vessels at the surface.
3) revise y, y_hat, nu, nus
4) Main subtree from aorta to RCA and LMA nodes.

MY UNITS ARE:
GRAMS, SECONDS, MM.
CONVERT FROM mmHg TO Pa using 133.3. 1 mmHg = 133.3 Pa
*/

/*_________________________________________________________________________________
____________________________________________________________________________________*/


#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.


int main(int argc,char **argv){
	
#include "declarations.c"	

int counted, counted_atOrder, stopOrder, startOrder;

// assign coordinates and all things along the way.
int 			numNodes2, xi_loc, xs_loc, xs_parent, xs_parents_parent;
int 			tempme, tempparent, temporder;
double 		templength, di_length, Ls, temp_termLsdi, nus_length, tempdistance, distToSurface_avbdr, 
			  distToSurface_LVCHAMBER, distToSurface_RVCHAMBER, distToSurface_OUTSIDE, cs, cb, spdotvd;
nodeInfo 	   *srkNodes2; // array struct.
nodeInfo	   *plotThisSegment;
XYZ 		xu_vector, xu_vectorp, di_vector, yihat_vector, nus_vector, nusunit_vector, avbdr_closest, avbdr_normal, LVCHAMBER_closest, LVCHAMBER_normal, 				RVCHAMBER_closest, RVCHAMBER_normal, OUTSIDE_closest, OUTSIDE_normal, wns_vector, direction_vector, new_loc, vd, nb, sp, vd1, vd2, new_loc_theta1, new_loc_theta2;

double 		rotation[3][3], ux, uy, uz, cost, sint, thetaT, fractionforOrder;

double maxLength_epi, maxLength_transmural;

int segLabel;

segLabel = -1;

// initialise PetSc and MPI. This is in the declarations.h
//  PetscInitialize(&argc,&argv,(char*)0,help);

if(debug)
  printf("Human bi-ventricle model geometry, half ellipsoids, 1 proc job. \n");

  /* Your vasculature is embedded into the 3D ventricle.   */
// writeVentricleInCoordinates();
 
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ------------------------------------------------------------------------*/
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ------------------------------------------------------------------------*/

/* Make your segments topologies for LAD, LCX, RCA. */
/* Morphometric data is taken from Kassab 1993. */

/* Several arrays are defined in this file: Do not change any of these arrays outside of this file for tractability. */
#include "connectivityKassab.c"
#include "lengthsRatiosKassab.c"

/* The whole tree is just 3 separate tree's, or LCX is an offshoot of LAD at the first opportunity and we have RCA and LA. */

if(debug){ printf("%d %d %d \n", usr_MX, usr_MY, usr_MZ); printf("%d \n", usr_MX * usr_MY * usr_MZ); }

 strcpy( command, "rm -f meandParent.txt" ); system(command);

// doing the segments iteratively is a pain: do it recursively by generating the topology first. root is already declared, but not allocated.

whichCase = atoi(argv[1]);
stopOrder = atoi(argv[2]);

if(whichCase!=4){ printf("The algorithm is best for case 4, the whole case. first arg needs to be 4 otherwise you are on your own.\n"); exit(0); }
// stop order dictates how topology is generated.
if(stopOrder<5){ printf("Do you know how long stopOrder < 5 will take? F off!\n"); exit(0); }

double epi_len_min, epi_len_max, trans_len_min, trans_len_max;
int min_num_nodesRCA, max_num_nodesRCA, min_num_nodesLAD, max_num_nodesLAD, min_num_nodesLCX, max_num_nodesLCX;
epi_len_min = 100000000000000000000.0;
epi_len_max = -1.0;
trans_len_min = 100000000000.0;
trans_len_max = -1.0;

/*
RCAcs =

                        26
                       111
                     299.3
                    832.82
                   2134.51
                   5487.91
                  20195.91
                  66366.11
                 180579.11
                 440113.11
                1179505.83


LAcs =

                        17
                       201
                    796.02
                   2257.26
                   5524.14
                 14308.006
                 50338.554
                161177.814
                346129.814
               1052134.454
               3434630.854
*/


// some values if order is different.
min_num_nodesRCA = 0; max_num_nodesRCA = 10000000; min_num_nodesLAD = 0; max_num_nodesLAD = 10000000; 
min_num_nodesLCX = 0; max_num_nodesLCX = 10000000; 

if(stopOrder==11){ min_num_nodesRCA = 26; max_num_nodesRCA = 26; min_num_nodesLAD = 17; max_num_nodesLAD = 17; 
min_num_nodesLCX = -1; max_num_nodesLCX = -1; 
}
if(stopOrder==10){ min_num_nodesRCA = 88; max_num_nodesRCA = 133; min_num_nodesLAD = 60; max_num_nodesLAD = 100; 
min_num_nodesLCX = 14; max_num_nodesLCX = 14; 
}
if(stopOrder==9){ min_num_nodesRCA = 240; max_num_nodesRCA = 360; min_num_nodesLAD = 243; max_num_nodesLAD = 365; 
min_num_nodesLCX = 64; max_num_nodesLCX = 100; 
}
if(stopOrder==8){ min_num_nodesRCA = 666; max_num_nodesRCA = 1000; 
// min_num_nodesLAD = 665; max_num_nodesLAD = 1000; 
min_num_nodesLAD = 665; max_num_nodesLAD = 1200; 
// min_num_nodesLCX = 236; max_num_nodesLCX = 354; 
min_num_nodesLCX = 236; max_num_nodesLCX = 500; 
}

if(stopOrder==7){ min_num_nodesRCA = 1700; max_num_nodesRCA = 2560; 
// min_num_nodesLAD = 1750; max_num_nodesLAD = 2625; 
min_num_nodesLAD = 1750; max_num_nodesLAD = 3500; 
// min_num_nodesLCX = 552; max_num_nodesLCX = 829; 
min_num_nodesLCX = 1552; max_num_nodesLCX = 2000; 
}
if(stopOrder==6){ min_num_nodesRCA = 4390; /* max_num_nodesRCA = 6585;  */ max_num_nodesRCA = 15000; 
// min_num_nodesLAD = 4223; max_num_nodesLAD = 6335; 
min_num_nodesLAD = 4223+4000; max_num_nodesLAD = 16000+10000; 
// min_num_nodesLCX = 1629; max_num_nodesLCX = 2444; 
min_num_nodesLCX = 2629; max_num_nodesLCX = 6000; 
}
if(stopOrder==5){ min_num_nodesRCA = 16156; max_num_nodesRCA = 24235; 
// min_num_nodesLAD = 14543; max_num_nodesLAD = 21800; 
min_num_nodesLAD = 14543; max_num_nodesLAD = 35800; 
min_num_nodesLCX = 5413; max_num_nodesLCX = 8120; 
}

// I am doing this only for case 4, all other cases to follow.
// decide on the epi length, transmural length, total number of nodes (=segments+1)
// there is a potential 20% pruning of stuff that cannot be stuffed into the anatomy.

#include "cases1_3.c"

if(whichCase==4){

if(stopOrder>=9){
	// in this case, it means nothing.
	trans_len_min = 0.0;
	trans_len_max = 0.0;	
}
else
if(stopOrder>=5){
	/* These values are given by Kaivomovitz in his appendix. It didnt work for me.
	trans_len_min = 9.0;
	trans_len_max = 14.5;	
	*/
	// pick a narrow range to make the topology a bit more selective.
	trans_len_min = 30;
	trans_len_max = 40;
}
	
	
// RCA root data is -1 and LA root data is -2, in case I need it later.	

	epi_len_min = 120.0+30;
	epi_len_max = 192.0+60;
	
	// RCA
	do{	
	for(temp_y = 0; temp_y < billion; temp_y++){ nodeData[temp_y] = temp_y + 1;   }	

	order = 11	; // this is as written in Kassab, C arrays rearranged.
	startOrder = order;
//	x0 = 61.0;    y0 = 46.0; z0 = MARGIN;
	x0 = 54.75;  y0 = 65.0; z0 = MARGIN;
	
	if(rootRCA!=NULL){ deleteTree(rootRCA); rootRCA = NULL; }
	
	if(debug)
	for(temp_y=0; temp_y < 12; temp_y++) printf("%d %f\n", temp_y, elementDiameters_meanRCA[temp_y]/2.0);
	
	rootRCA = makeKassabTree((int)(segmentElementRatio_meanRCA[order]), nodeData, order, segmentElementRatio_meanRCA, segmentElementRatio_SDRCA, cuprob_connRCA, segmentLengths_meanRCA,  segmentLengths_SDRCA, elementDiameters_meanRCA, -1, 1, stopOrder, startOrder, amx, x0, y0, z0, -1 , -1, 0.0, -1, segLabel);
	// count and reconstruct the tree till you get enough but not too much number of nodes.
	counted = 0; counted = binarytree_count(rootRCA);

	// remove right of the RCA root.
	if(rootRCA->right!=NULL){ deleteTree(rootRCA->right); rootRCA->right = NULL; }
	
	// count number of nodes of order X coming off element of order Y.
	counted_atOrder = 0; counted_atOrder = binarytree_count_atOrder(rootRCA, 11, 10);
	
	if(debug)
	printf("Number of nodes in RCA, case 4, tree: %d\n", counted );
	output = fopen("nodesGenRCA.data","w"); fprintf(output,"%d\n", counted); fclose(output);
	maxLength_epi = maxLength_transmural = 0.0; print_pathsL(rootRCA, startOrder, 9);
	output = fopen("pathLength.data","r"); while(fscanf(output, "%lf", &viz_leng)!=EOF) if(maxLength_epi < viz_leng) maxLength_epi = viz_leng;  fclose(output);                
	strcpy( command, "rm -f pathLength.data" ); system(command);
	print_pathsL(rootRCA, 8, stopOrder); output=fopen("pathLength.data","r"); 
	while(fscanf(output,"%lf",&viz_leng)!=EOF) if(maxLength_transmural<viz_leng) maxLength_transmural=viz_leng;fclose(output);
	strcpy( command, "rm -f pathLength.data" ); system(command);
	if(debug)
		printf("max epi length of RCA tree is: %f; length of transmural %f; number of O(10): %d \n", maxLength_epi, maxLength_transmural, counted_atOrder);
	
	printf("RCA tree: %d %f %f %f; %f %f %f %d %d %d %d\n", stopOrder, epi_len_min, maxLength_epi, epi_len_max, trans_len_min, maxLength_transmural, trans_len_max, min_num_nodesRCA, counted, max_num_nodesRCA, counted_atOrder );
	}while( (maxLength_epi<epi_len_min||maxLength_epi>epi_len_max) || ( (stopOrder < 9) && ( (maxLength_transmural < trans_len_min) || (maxLength_transmural > trans_len_max) ) ) || ( (min_num_nodesRCA > counted) || (max_num_nodesRCA < counted) ) ||  (counted_atOrder < 5 || counted_atOrder > 8) );
// }while(1==0);

/*_________________________________________________________________________________________*/	
	// nodeData has now changed. Copy it to summat.	
	for(temp_y = 0; temp_y < billion; temp_y++){ usednodeData[temp_y] = nodeData[temp_y];   }
/*_________________________________________________________________________________________*/

	epi_len_min = 100.0+35;
	epi_len_max = 160.0+35;

if(stopOrder>=9){
	// in this case, it means nothing.
	trans_len_min = 0.0;
	trans_len_max = 0.0;	
}
else
if(stopOrder>=5){
	/* These values are given by Kaivomovitz in his appendix. It didnt work for me.
	trans_len_min = 9.0;
	trans_len_max = 14.5;	
	*/
	// pick a narrow range to make the topology a bit more selective.
	trans_len_min = 30;
	trans_len_max = 45;
}
	
	
	// LA: LAD.
do{
	
for(temp_y = 0; temp_y < billion; temp_y++){ nodeData[temp_y] = usednodeData[temp_y];   }
 order = 11; startOrder = order;
 	x0 = 61.5;  y0 = 36.0; z0 = 2.0;
if(rootLA!=NULL){ deleteTree(rootLA); rootLA = NULL; }
rootLA = makeKassabTree((int)(segmentElementRatio_meanLAD[order]), nodeData, order, segmentElementRatio_meanLAD, segmentElementRatio_SDLAD, cuprob_connLAD, segmentLengths_meanLAD,  segmentLengths_SDLAD, elementDiameters_meanLAD, -2, 2, stopOrder, startOrder, amx, x0, y0, z0, -1 , -1, 0.0, -1, segLabel);

	// remove right of the LA root.
	if(rootLA->right!=NULL){ deleteTree(rootLA->right); rootLA->right = NULL; }

// now the idea is to attach LCX at the root->right location.
/*
	if(stopOrder < 11){
		if(rootLA->right!=NULL){ deleteTree(rootLA->right);  }
		order = 10; startOrder = order;
	 	x0 = 62.0;  y0 = 31.0; z0 = 3.0;
		rootLA->right = makeKassabTree((int)(segmentElementRatio_meanLCX[order]), nodeData, order, segmentElementRatio_meanLCX, segmentElementRatio_SDLCX, cuprob_connLCX, segmentLengths_meanLCX, segmentLengths_SDLCX, elementDiameters_meanLCX, rootLA->data, 3, stopOrder, startOrder , amx , x0, y0, z0, -1 , -1, 0.0, rootLA->strahler, segLabel);
	}
*/

counted = 0; counted = binarytree_count(rootLA); /* printf("Number of nodes in LAD, case 4, tree: %d\n", counted ); */

	// count number of nodes of order X coming off element of order Y.
	counted_atOrder = 0; counted_atOrder = binarytree_count_atOrder(rootLA, 11, 10);


output = fopen("nodesGenLAD.data","w"); fprintf(output,"%d\n", counted); fclose(output);
maxLength_epi = maxLength_transmural = 0.0;
print_pathsL(rootLA, startOrder, 9);
output = fopen("pathLength.data","r"); while(fscanf(output, "%lf", &viz_leng)!=EOF) if(maxLength_epi < viz_leng) maxLength_epi = viz_leng;  fclose(output);                
strcpy( command, "rm -f pathLength.data" ); system(command);
print_pathsL(rootLA, 8, stopOrder);
output=fopen("pathLength.data","r"); while(fscanf(output,"%lf",&viz_leng)!=EOF) if(maxLength_transmural<viz_leng) maxLength_transmural=viz_leng;fclose(output);
strcpy( command, "rm -f pathLength.data" ); system(command);
if(debug)
printf("max epi length of LAD tree is: %f; length of transmural %f \n", maxLength_epi, maxLength_transmural);

printf("LAD tree: %d %f %f %f; %f %f %f %d %d %d %d\n", stopOrder, epi_len_min, maxLength_epi, epi_len_max, trans_len_min, maxLength_transmural, trans_len_max, min_num_nodesLAD, counted, max_num_nodesLAD, counted_atOrder );

 } while( (maxLength_epi<epi_len_min||maxLength_epi>epi_len_max) || 	      ( (stopOrder < 9) && ( (maxLength_transmural < trans_len_min) || (maxLength_transmural > trans_len_max) ) ) || ( (min_num_nodesLAD > counted) || (max_num_nodesLAD < counted) ) ||  (counted_atOrder < 5 || counted_atOrder > 8)   );
 
// }while(1==0);

// LCX 
	epi_len_min = 100.0;
	epi_len_max = 160.0;
/*_________________________________________________________________________________________*/	
	// nodeData has now changed. Copy it to summat.	
	for(temp_y = 0; temp_y < billion; temp_y++){ usednodeData[temp_y] = nodeData[temp_y];   }
/*_________________________________________________________________________________________*/

if(stopOrder < 11)
do{
	
for(temp_y = 0; temp_y < billion; temp_y++){ nodeData[temp_y] = usednodeData[temp_y];   }

if(rootLA->left->left->left->left->left->right!=NULL) { deleteTree(rootLA->left->left->left->left->left->right); rootLA->left->left->left->left->left->right = NULL;  }
order = 10; startOrder = order;
x0 = 62.0;  y0 = 31.0; z0 = 2.0;

rootLA->left->left->left->left->left->right = makeKassabTree((int)(segmentElementRatio_meanLCX[order]), nodeData, order, segmentElementRatio_meanLCX, segmentElementRatio_SDLCX, cuprob_connLCX, segmentLengths_meanLCX, segmentLengths_SDLCX, elementDiameters_meanLCX, rootLA->left->left->left->left->left->data, 3, stopOrder, startOrder , amx , x0, y0, z0, -1 , -1, 0.0, rootLA->left->left->left->left->left->strahler, segLabel);

counted = 0; counted = binarytree_count(rootLA->left->left->left->left->left->right); /* printf("Number of nodes in LAD, case 4, tree: %d\n", counted ); */
output = fopen("nodesGenLCX.data","w"); fprintf(output,"%d\n", counted); fclose(output);
maxLength_epi = maxLength_transmural = 0.0;
print_pathsL(rootLA->left->left->left->left->left->right, startOrder, 9);
output = fopen("pathLength.data","r"); while(fscanf(output, "%lf", &viz_leng)!=EOF) if(maxLength_epi < viz_leng) maxLength_epi = viz_leng;  fclose(output);                
strcpy( command, "rm -f pathLength.data" ); system(command);
print_pathsL(rootLA->left->left->left->left->left->right, 8, stopOrder);
output=fopen("pathLength.data","r"); while(fscanf(output,"%lf",&viz_leng)!=EOF) if(maxLength_transmural<viz_leng) maxLength_transmural=viz_leng;fclose(output);
strcpy( command, "rm -f pathLength.data" ); system(command);
if(debug)
printf("max epi length of LCX tree is: %f; length of transmural %f \n", maxLength_epi, maxLength_transmural);

	counted_atOrder = 0; counted_atOrder = binarytree_count_atOrder(rootLA->left->left->left->left->left->right, 10, 9);

printf("LCX tree: %d %f %f %f; %f %f %f %d %d %d %d\n", stopOrder, epi_len_min, maxLength_epi, epi_len_max, trans_len_min, maxLength_transmural, trans_len_max, min_num_nodesLCX, counted, max_num_nodesLCX, counted_atOrder );

 } while( (maxLength_epi<epi_len_min||maxLength_epi>epi_len_max) || ( (stopOrder < 9) && ( (maxLength_transmural < trans_len_min) || (maxLength_transmural > trans_len_max) ) ) || ( (min_num_nodesLCX > counted) || (max_num_nodesLCX < counted) )   );
// }while(1==0); 

 // and the onoff-assigned didnt work. So do that in a new recursion.
onOffAssignedRecurse(rootRCA); onOffAssignedRecurse(rootLA); 
// revise the lengths, resistances of the nodes where you already have coordinates.
x0 = -1;  y0 = -1; z0 = -1; reviseLengthResistance(rootRCA, x0, y0, z0);  
x0 = -1;  y0 = -1; z0 = -1; reviseLengthResistance(rootLA, x0, y0, z0); 
	
} // end of whichCase = 4.

printf("constructed trees.\n");

//=====================================================================================================================
// calculate resistance, perfusion, and write a data file.
if(whichCase>0 && whichCase< 4){
	// assign downResistance.
	totalTreeResistance = binaryDownResistance3D(root);
	printf("down resistances: %10.10f %10.10f\n", totalTreeResistance, root->downResistance);
	// pressure and perfusion.
	root->pressure 	= Pin; // g/(m-s)
	Pout = Pin - Q_perf * totalTreeResistance;
	root->perfusion = (Pin - Pout)/(root->downResistance); // In case of individual trees, this works out to Q_perf: see line above.
	root->theta = 0.0;
	perfusionPressure3D(root, Pout);
	writeFile(root);
	printf("%f %f\n", Q_perf, root->perfusion);
}

if(whichCase==4){
// redo lengths here somewhere.	
	doLengths(rootRCA); doLengths(rootLA);
	
	totalTreeResistanceRCA  = binaryDownResistance3D(rootRCA); totalTreeResistanceLA     = binaryDownResistance3D(rootLA);
	rootRCA->pressure 	= Pin; rootLA->pressure = Pin; // Units: kg/(m-s2)

		if(totalTreeResistanceRCA + totalTreeResistanceLA > 0.0)
		totalTreeResistance = totalTreeResistanceRCA * totalTreeResistanceLA / (totalTreeResistanceRCA + totalTreeResistanceLA);
		else
		printf("your LA or RCA or both resistances are nada, check!");

		Pout = Pin - Q_perf * totalTreeResistance; // this Pout is the pressure at the terminals: since they are SN 5 or 6, I dont know this before.

// I potentially havent understood how to assign the pressure and perfusion boundary conditions!!!!!! So Keeping it simple for now.
// rootRCA->perfusion = Q_perf * totalTreeResistanceLA    / totalTreeResistance ;  rootLA->perfusion = Q_perf * totalTreeResistanceRCA / totalTreeResistance ;
rootRCA->perfusion = (Pin - Pout)/ totalTreeResistanceRCA;  rootLA->perfusion = (Pin - Pout) / totalTreeResistanceLA ;
rootRCA->theta 	   = 0.0; 												 rootLA->theta  	 = 0.0;
printf("In main: %f %f %f %f %f %f\n", Q_perf, rootRCA->perfusion, rootLA->perfusion,  rootRCA->perfusion + rootLA->perfusion, Pin, Pout );
printf("In main, resistances: %f %f %f\n", rootRCA->downResistance, rootLA->downResistance , totalTreeResistance );

perfusionPressure3D(rootRCA, Pout);   									perfusionPressure3D(rootLA, Pout);
// writeFile(rootRCA);													writeFile(rootLA);
} // end of whichCase == 4.

 
 // find out which O(11) and O(10) in LCX nodes are NOT defined.
//mainNodesNotDefined(rootLA);
//mainNodesNotDefined(rootRCA);



//=====================================================================================================================

// do this when testing initial branches. Do not do this at run time - the nodes are written twice, cannot read it into my visualisation.
//writeTreeSoFar(rootRCA);
//writeTreeSoFar(rootLA);

//  if(1==0){

// for(thisOrder = 11; thisOrder >= stopOrder; thisOrder--){
	// optimisation and everything is here.
	xu_vectorp.x = xu_vectorp.y = xu_vector.z = -1000.0;
	// fractionforOrder = 0.01; // assuming we always start from orders more than 10 in the root.
	// at least the root or roots must have coordinates and be assignedOrNot = 1 to start with.
	setCoordinatesofTree(rootRCA, rootRCA, rootLA, xu_vectorp, stopOrder, thisOrder);
	xu_vectorp.x = xu_vectorp.y = xu_vector.z = -1000.0;
	// fractionforOrder = 0.01; // assuming we always start from orders more than 10 in the root.
	// at least the root or roots must have coordinates and be assignedOrNot = 1 to start with.
	setCoordinatesofTree(rootLA,    rootRCA, rootLA, xu_vectorp, stopOrder, thisOrder);
// } // end of this order.

// }
 
// =========== end of optimisation. =========================================================================================================
// ==========================================================================================================================================

// RCA part. =================================================================================================================================
int temp_me_data;

// if this pruning is successful, then you have to write the rawSegs again.
 pruneTree2(rootRCA); pruneTree2(rootLA); // this does not work yet.
strcpy( command, "mv rawSegsRCA.data rawSegsRCA.data.old" ); system(command);
strcpy( command, "mv rawSegsLA.data rawSegsLA.data.old" ); system(command);

// you need to recalculate.
if(whichCase==4){
// redo lengths here somewhere.	
	doLengths(rootRCA); doLengths(rootLA);
	
	totalTreeResistanceRCA  = binaryDownResistance3D(rootRCA); totalTreeResistanceLA     = binaryDownResistance3D(rootLA);
	rootRCA->pressure 	= Pin; rootLA->pressure = Pin; // Units: kg/(m-s2)

		if(totalTreeResistanceRCA + totalTreeResistanceLA > 0.0)
		totalTreeResistance = totalTreeResistanceRCA * totalTreeResistanceLA / (totalTreeResistanceRCA + totalTreeResistanceLA);
		else
		printf("your LA or RCA or both resistances are nada, check!");

		Pout = Pin - Q_perf * totalTreeResistance; // this Pout is the pressure at the terminals: since they are SN 5 or 6, I dont know this before.

// I potentially havent understood how to assign the pressure and perfusion boundary conditions!!!!!! So Keeping it simple for now.
// rootRCA->perfusion = Q_perf * totalTreeResistanceLA    / totalTreeResistance ;  rootLA->perfusion = Q_perf * totalTreeResistanceRCA / totalTreeResistance ;
rootRCA->perfusion = (Pin - Pout)/ totalTreeResistanceRCA;  rootLA->perfusion = (Pin - Pout) / totalTreeResistanceLA ;
rootRCA->theta 	   = 0.0; 												 rootLA->theta  	 = 0.0;
printf("After pruning: %f %f %f %f %f %f\n", Q_perf, rootRCA->perfusion, rootLA->perfusion,  rootRCA->perfusion + rootLA->perfusion, Pin, Pout );
printf("After pruning::: %f %f %f\n", rootRCA->downResistance, rootLA->downResistance , totalTreeResistance );

perfusionPressure3D(rootRCA, Pout);   									perfusionPressure3D(rootLA, Pout);
} // end of whichCase == 4.


// make a new copy of raw segs.
 writeTreeSoFar(rootRCA);  writeTreeSoFar(rootLA);

// assign temp_z.
temp_z = 0;
output = fopen("rawSegsRCA.data","r");

if(output==NULL){
printf("no rawsegs RCA file. moving on.\n");
}
else
{
	
	
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF) temp_z++;
fclose(output);

printf("Number of nodes in RCA sub-tree: %d\n", temp_z);

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
temp_z = 0;
output = fopen("rawSegsLA.data","r");
if(output==NULL){
printf("no raw LA file here, moving on.\n");
}
else
{
	
while(fscanf(output, "%lf %lf %lf %d %d %lf %lf %lf %lf %d %d", &tx, &ty, &tz, &intp1, &intp2, &tRad, &viz_resistivity, &viz_pressure, &viz_perfusion, &order, &temp_order)!=EOF) temp_z++;
fclose(output);

printf("Number of nodes in LA sub-tree: %d\n", temp_z);

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

// assuming we got this far, prune it.
//pruneTree(rootRCA);
//pruneTree(rootLA);


// RCA tree
sprintf(str, "inorderRCA.txt");
printInorder(rootRCA, str);
sprintf(str, "preorderRCA.txt");
printPreorder(rootRCA, str);

// LA tree.
sprintf(str, "inorderLA.txt");
printInorder(rootLA, str);
sprintf(str, "preorderLA.txt");
printPreorder(rootLA, str);

printf("did in order and preorder in some way - it has -1's.\n'");
// now read it in.


/*
Experiment codes here.
In all cases, you have to work out the resistance, perfusion, and pressure and write a new distribution file similar to rawSegs
1) reduction of large artery and small artery radius: 5%, 10%, 20%, 25%.
2) Occlusion of RCA, LAD or LCX (one at a time)
3) effect of values of HCT: 0.15, 0.20, 0.25, 0.3, 0.35, 0.4, 0.45
*/
	
	
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ------------------------------------------------------------------------
     Free work space.
   ---------------------------------------------------------------------- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
}


