// these are the 3 cases in topology construction I am not using at production stage. So just put the code into a separate file.

/* top order to stopOrder 4 can be generated easily. order 4 segments are already ~ 0.1 mm, i.e. 5 times the resolution of the CT. */
// the nodeData  array needs to be reset at each new call of makeKassabTree.
if(whichCase==1){
	
do{	
for(temp_y = 0; temp_y < billion; temp_y++){ nodeData[temp_y] = temp_y + 1;   }	

order = 11	; // this is as written in Kassab, C arrays rearranged.
startOrder = order;

if(debug)
printf("%d in main.\n", (int)segmentElementRatio_meanRCA[order]);
// last entry is the parent. Parent of root is -1 i.e. non-existent.
// make sure you make the parents something else when constructing whole tree.
if(root!=NULL) deleteTree(root);
root = makeKassabTree((int)(segmentElementRatio_meanRCA[order]), nodeData, order, segmentElementRatio_meanRCA, segmentElementRatio_SDRCA, cuprob_connRCA, segmentLengths_meanRCA,  segmentLengths_SDRCA, elementDiameters_meanRCA, -1, whichCase, stopOrder, startOrder, amx, -1.0, -1.0, -1.0, -1, -1, 0.0, -1, segLabel);
// count and reconstruct the tree till you get enough but not too much number of nodes.
printf("Number of nodes in this tree: %d\n", binarytree_count(root) );
counted = 0;
counted = binarytree_count(root);

output = fopen("nodesGenRCA.data","a+"); fprintf(output,"%d\n", counted); fclose(output);

maxLength_epi = maxLength_transmural = 0.0;
// epicardial vessels/elements total path length.
print_pathsL(root, startOrder, 9);
// open the data file and work out the maxLength_epi
output = fopen("pathLength.data","r"); while(fscanf(output, "%lf", &viz_leng)!=EOF) if(maxLength_epi < viz_leng) maxLength_epi = viz_leng;  fclose(output);                
// get rid of this file immediately.
strcpy( command, "rm -f pathLength.data" ); system(command);
// transmural vessels/elements total path length.
print_pathsL(root, 8, 5);
// open the data file and work out the maxLength_epi
output=fopen("pathLength.data","r"); while(fscanf(output,"%lf",&viz_leng)!=EOF) if(maxLength_transmural<viz_leng) maxLength_transmural=viz_leng;fclose(output);
// get rid of this file NOW.
strcpy( command, "rm -f pathLength.data" ); system(command);
if(debug)
printf("max epi length of RCA tree is: %f; length of transmural %f \n", maxLength_epi, maxLength_transmural);

if(debug)
printf("finished making RCA tree.\n");
// } while( (counted < (9000-2000) || counted > (9000+2000) ) || (maxLength_epi< 120.0 || maxLength_epi > 192.0 ) || maxLength_transmural < 5 || maxLength_transmural > 40 ); // this is for SN from 11 to 5: 9000 pm 700.



} while( (maxLength_epi < 20.0 || maxLength_epi > 192.0 ) ); 

} // end of whichCase = 1.

// whichCase = 2;
if(whichCase==2){
	
do{
	
for(temp_y = 0; temp_y < billion; temp_y++){ nodeData[temp_y] = temp_y + 1;   }
 order = 11; startOrder = order;
// if(debug)
printf("%d in main.\n", (int)segmentElementRatio_meanLAD[order]);
if(root!=NULL) deleteTree(root);
root = makeKassabTree((int)(segmentElementRatio_meanLAD[order]), nodeData, order, segmentElementRatio_meanLAD, segmentElementRatio_SDLAD, cuprob_connLAD, segmentLengths_meanLAD,  segmentLengths_SDLAD, elementDiameters_meanLAD, -1, whichCase, stopOrder, startOrder, amx, -1.0, -1.0, -1.0, -1, -1, 0.0, -1, segLabel);
printf("Number of nodes in this tree: %d\n", binarytree_count(root) );
counted = 0;
counted = binarytree_count(root);

output = fopen("nodesGenLAD.data","a+"); fprintf(output,"%d\n", counted); fclose(output);

maxLength_epi = maxLength_transmural = 0.0;
// epicardial vessels/elements total path length.
print_pathsL(root, startOrder, 9);
// open the data file and work out the maxLength_epi
output = fopen("pathLength.data","r"); while(fscanf(output, "%lf", &viz_leng)!=EOF) if(maxLength_epi < viz_leng) maxLength_epi = viz_leng;  fclose(output);                
// get rid of this file immediately.
strcpy( command, "rm -f pathLength.data" ); system(command);
// transmural vessels/elements total path length.
print_pathsL(root, 8, 5);
// open the data file and work out the maxLength_epi
output=fopen("pathLength.data","r"); while(fscanf(output,"%lf",&viz_leng)!=EOF) if(maxLength_transmural<viz_leng) maxLength_transmural=viz_leng;fclose(output);
// get rid of this file NOW.
strcpy( command, "rm -f pathLength.data" ); system(command);
if(debug)
printf("max epi length of LAD tree is: %f; length of transmural %f \n", maxLength_epi, maxLength_transmural);

if(debug)
printf("finished making LAD tree.\n");
 } while(counted < (8000-2000) || counted > (8000+2000) || (maxLength_epi< 100.0 || maxLength_epi > 160.0 )  || maxLength_transmural < 5 || maxLength_transmural > 40 );  // this is for SN from 11 to 5: 8000 pm 1100.

} // end of whichCase = 2.

// whichCase = 3;
if(whichCase==3){

do{

for(temp_y = 0; temp_y < billion; temp_y++){ nodeData[temp_y] = temp_y + 1;   }	
order = 10;
startOrder = order;
// if(debug)
printf("%d in main.\n", (int)segmentElementRatio_meanLCX[order]);	
if(root!=NULL) deleteTree(root);
root = makeKassabTree((int)(segmentElementRatio_meanLCX[order]), nodeData, order, segmentElementRatio_meanLCX, segmentElementRatio_SDLCX, cuprob_connLCX, segmentLengths_meanLCX, segmentLengths_SDLCX, elementDiameters_meanLCX, -1, whichCase, stopOrder, startOrder, amx , -1.0, -1.0, -1.0, -1, -1, 0.0, -1, segLabel);
printf("Number of nodes in this tree: %d\n", binarytree_count(root) );
counted = 0;
counted = binarytree_count(root);

output = fopen("nodesGenLCX.data","a+");  fprintf(output,"%d\n", counted); fclose(output);

maxLength_epi = maxLength_transmural = 0.0;
// epicardial vessels/elements total path length.
print_pathsL(root, startOrder, 9);
// open the data file and work out the maxLength_epi
output = fopen("pathLength.data","r"); while(fscanf(output, "%lf", &viz_leng)!=EOF) if(maxLength_epi < viz_leng) maxLength_epi = viz_leng;  fclose(output);                
// get rid of this file immediately.
strcpy( command, "rm -f pathLength.data" ); system(command);
// transmural vessels/elements total path length.
print_pathsL(root, 8, 5);
// open the data file and work out the maxLength_epi
output=fopen("pathLength.data","r"); while(fscanf(output,"%lf",&viz_leng)!=EOF) if(maxLength_transmural<viz_leng) maxLength_transmural=viz_leng;fclose(output);
// get rid of this file NOW.
strcpy( command, "rm -f pathLength.data" ); system(command);
if(debug)
printf("max epi length of LCX tree is: %f; length of transmural %f \n", maxLength_epi, maxLength_transmural);


if(debug)
printf("finished making LCX tree.\n");
} while(counted < (3000-1000) || counted > (3000+1000) || (maxLength_epi< 100.0 || maxLength_epi > 160.0 )  || maxLength_transmural < 5 || maxLength_transmural > 40 ); // this is for SN from 10 to 5.

}


