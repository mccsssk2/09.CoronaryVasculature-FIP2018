/*
EPI ENDO average over all 500 instances. This should be integrated with srkvascExpts.c.
17 Nov 2017.
*/

#include "srkVasculature.h"
// Kirk cannot take anything more than 10^9.

int main(int argc,char **argv){
	
#include "declarations.c"	
int fileNum;

double epiendoRV[2][6], epiendoLV[2][12];

for(o1=0;o1<6;o1++){ epiendoRV[0][o1] = 0; 	epiendoRV[1][o1] = 0; }
for(o1=0;o1<12;o1++){ epiendoLV[0][o1] = 0; 	epiendoLV[1][o1] = 0; }

for(fileNum = 1; fileNum<=500; fileNum++){

		sprintf(str, "epiendoRV0_%d.data", fileNum);
		output = fopen(str,"r");
		if(output==NULL){ printf("could not open RV file, %d. exiting.\n", fileNum); exit(0); }
		while(fscanf(output, "%lf %lf", &viz_leng, &viz_perfusion)!=EOF){
			cn = (int)(viz_leng);
			epiendoRV[0][cn] = viz_leng;
			epiendoRV[1][cn] = epiendoRV[1][cn] + viz_perfusion;
		}


		sprintf(str, "epiendoLV0_%d.data", fileNum);
		output = fopen(str,"r");
		if(output==NULL){ printf("could not open LV file, %d. exiting.\n", fileNum); exit(0); }
		while(fscanf(output, "%lf %lf", &viz_leng, &viz_perfusion)!=EOF){
			cn = (int)(viz_leng);
			epiendoLV[0][cn] = viz_leng;
			epiendoLV[1][cn] = epiendoLV[1][cn] + viz_perfusion;
		}
} // end of fileNum loop.

// average it up and write to a final.

for(o1=0;o1<6;o1++){ epiendoRV[1][o1] = epiendoRV[1][o1]/(double)500.0; }
for(o1=0;o1<12;o1++){ epiendoLV[1][o1] = epiendoLV[1][o1]/(double)500.0; }

// write
output = fopen("rv.dataT","w");
for(o1=0;o1<6;o1++){ fprintf(output, "%d %f %f\n", o1, epiendoRV[0][o1]+0.5, epiendoRV[1][o1]); }
fclose(output);
output = fopen("lv.dataT","w");
for(o1=0;o1<12;o1++){ fprintf(output, "%d %f %f\n", o1, epiendoLV[0][o1]+0.5, epiendoLV[1][o1]); }
fclose(output);

return 0;
}
