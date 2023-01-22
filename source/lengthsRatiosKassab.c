/* !!* *!! MAKE SURE THE WHICHCASE VARIABLE IS ASSIGNED BEFORE YOU RUN THIS SEGMENT OF THE CODE.
*/

// RCA, table 1 and table 5.
	segmentElementRatio_meanRCA[0] =  0; 
	segmentElementRatio_meanRCA[1] = 1.88; 
	segmentElementRatio_meanRCA[2] =      1.88; 
	segmentElementRatio_meanRCA[3] =       2.2; 
	segmentElementRatio_meanRCA[4] =      2.3; 
	segmentElementRatio_meanRCA[5] =       2.0; 
	segmentElementRatio_meanRCA[6] =       2.3; 
	segmentElementRatio_meanRCA[7] =   3.23; 
	segmentElementRatio_meanRCA[8] =   4.68; 
	segmentElementRatio_meanRCA[9] =   5.38; 
	segmentElementRatio_meanRCA[10] =    8.5; 
	segmentElementRatio_meanRCA[11] =      26+7; 

	 segmentElementRatio_SDRCA[0]      = 0; 	
	 segmentElementRatio_SDRCA[1]      = 0.99; 
	 segmentElementRatio_SDRCA[2]      =     1.0; 
	 segmentElementRatio_SDRCA[3]      =         1.2; 
	 segmentElementRatio_SDRCA[4]      =       1.8; 
	 segmentElementRatio_SDRCA[5]      =       0.9; 
	 segmentElementRatio_SDRCA[6]      =       1.3; 
	 segmentElementRatio_SDRCA[7]      =    2.1; 
	 segmentElementRatio_SDRCA[8]      =     2.7; 
	 segmentElementRatio_SDRCA[9]      =      3.6; 
	 segmentElementRatio_SDRCA[10]      =      0.0; 
	 segmentElementRatio_SDRCA[11]      =    0.0;  

 // lengths of segments and standard deviations in mm.
segmentLengths_meanRCA[0] =  0.0  ;   segmentLengths_SDRCA[0] =   0.0;
segmentLengths_meanRCA[1] =  0.069  ;   segmentLengths_SDRCA[1] =   0.046;
segmentLengths_meanRCA[2] =  0.083  ;   segmentLengths_SDRCA[2] =   0.070;
segmentLengths_meanRCA[3] =  0.085  ;   segmentLengths_SDRCA[3] =   0.061;
segmentLengths_meanRCA[4] =  0.118  ;   segmentLengths_SDRCA[4] =   0.113;
segmentLengths_meanRCA[5] =  0.449  ;   segmentLengths_SDRCA[5] =   0.350;
segmentLengths_meanRCA[6] =  0.748  ;   segmentLengths_SDRCA[6] =   0.654;
segmentLengths_meanRCA[7] =  0.986  ;   segmentLengths_SDRCA[7] =   0.810;
segmentLengths_meanRCA[8] =  1.26   ;   segmentLengths_SDRCA[8] =   1.10;
segmentLengths_meanRCA[9] =  1.62   ;   segmentLengths_SDRCA[9] =   1.31;
segmentLengths_meanRCA[10] =  1.89   ;   segmentLengths_SDRCA[10] =   1.38;
segmentLengths_meanRCA[11] =  3.24   ;   segmentLengths_SDRCA[11] =   2.09;

	 elementLengths_meanRCA[0] 	         = 0 ; // mm
	 elementLengths_meanRCA[1] 	         = 0.125 ; // mm
	 elementLengths_meanRCA[2] 	         =  0.141; // mm
	 elementLengths_meanRCA[3] 	         =    0.178; // mm
	 elementLengths_meanRCA[4] 	         =  0.253; // mm
	 elementLengths_meanRCA[5] 	         =  0.545; // mm
	 elementLengths_meanRCA[6] 	         =    1.64 ; // mm
	 elementLengths_meanRCA[7] 	         =  3.13; // mm
	 elementLengths_meanRCA[8] 	         =    5.99; // mm
	 elementLengths_meanRCA[9] 	         =   9.06; // mm
	 elementLengths_meanRCA[10] 	         =    16.1; // mm
	 elementLengths_meanRCA[11] 	         =    84.1; // this is mm

	 elementLengths_SDRCA[0] 		   = 0.0; // mm
	 elementLengths_SDRCA[1] 		   = 0.084; // mm
	 elementLengths_SDRCA[2] 		   =  0.103; // mm
	 elementLengths_SDRCA[3] 		   =  0.105; // mm
	 elementLengths_SDRCA[4] 		   =  0.174; // mm
	 elementLengths_SDRCA[5] 		   =   0.451; // mm
	 elementLengths_SDRCA[6] 		   =     1.13; // mm
	 elementLengths_SDRCA[7] 		   =    2.11; // mm
	 elementLengths_SDRCA[8] 		   =    3.53; // mm
	 elementLengths_SDRCA[9] 		   =    5.56; // mm
	 elementLengths_SDRCA[10] 		   =    13.3; // mm
	 elementLengths_SDRCA[11] 		   =    0.0; // mm

	 elementDiameters_meanRCA[0] 	   = 0; // micrometers
	 elementDiameters_meanRCA[1] 	   = 9.3; // micrometers
	 elementDiameters_meanRCA[2] 	   =       12.8; // micrometers
	 elementDiameters_meanRCA[3] 	   =       17.7; // micrometers
	 elementDiameters_meanRCA[4] 	   =    28.6; // micrometers
	 elementDiameters_meanRCA[5] 	   =     63.1; // micrometers
	 elementDiameters_meanRCA[6] 	   =      132; // micrometers
	 elementDiameters_meanRCA[7] 	   =    256; // micrometers
	 elementDiameters_meanRCA[8] 	   =   428; // micrometers
	 elementDiameters_meanRCA[9] 	   =    706; // micrometers
	 elementDiameters_meanRCA[10] 	   =    1302; // micrometers
	 elementDiameters_meanRCA[11] 	   =  3218; // micrometers

// LAD, table 2 and table 5.
// I am using segmentElementRatio_mean as number of nodes I have to generate.
// for n segments, I need n+1 nodes.
segmentElementRatio_meanLAD[0] = 0; 
segmentElementRatio_meanLAD[1] = 2.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[2] =       1.79; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[3] =       2.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[4] =      2.28; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[5] =       2.02; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[6] =       2.23; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[7] =   3.89; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[8] =   4.69; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[9] =   6.06; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[10] =    9.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
segmentElementRatio_meanLAD[11] =      17+5; // micrometers

 segmentElementRatio_SDLAD[0]      = 0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[1]      = 1.4; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[2]      =  0.95; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[3]      =   1.1; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[4]      =   1.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[5]      =   1.2; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[6]      =   1.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[7]      =   2.1; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[8]      =   3.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[9]      =   4.2; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[10]      =   7.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLAD[11]      =   0.0;  // micrometers

 elementLengths_meanLAD[0] 	          = 0.0 ; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[1] 	          = 0.115 ; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[2] 	          =  0.136; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[3] 	          =    0.149; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[4] 	          =    0.353; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[5] 	          =     0.502; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[6] 	          =    1.31 ; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[7] 	          =    3.54; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[8] 	          =    4.99; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[9] 	          =   9.03; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[10] 	          =    20.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLAD[11] 	          =  47.9; // this is mm

 elementLengths_SDLAD[0] 		   = 0.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[1] 		   = 0.066; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[2] 		   =  0.088; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[3] 		   =  0.104; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[4] 		   =  0.154; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[5] 		   =  0.349; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[6] 		   =  0.914; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[7] 		   =  2.11; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[8] 		   =  3.02; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[9] 		   =  6.13; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[10] 		   =  17.9; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLAD[11] 		   =  0.0; // mm

segmentLengths_meanLAD[0] = 0.0 ;   segmentLengths_SDLAD[0] =  0.0;
segmentLengths_meanLAD[1] = 0.056 ;   segmentLengths_SDLAD[1] =  0.038;
segmentLengths_meanLAD[2] = 0.072 ;   segmentLengths_SDLAD[2] =  0.045;
segmentLengths_meanLAD[3] = 0.072 ;   segmentLengths_SDLAD[3] =  0.049;
segmentLengths_meanLAD[4] = 0.112 ;   segmentLengths_SDLAD[4] =  0.10;
segmentLengths_meanLAD[5] = 0.4542 ;   segmentLengths_SDLAD[5] =  0.33;
segmentLengths_meanLAD[6] = 0.609 ;   segmentLengths_SDLAD[6] =  0.48;
segmentLengths_meanLAD[7] = 0.9202 ;   segmentLengths_SDLAD[7] =  0.79;
segmentLengths_meanLAD[8] = 1.09 ;   segmentLengths_SDLAD[8] =  0.83;
segmentLengths_meanLAD[9] = 1.54 ;   segmentLengths_SDLAD[9] =  1.25;
segmentLengths_meanLAD[10] = 2.26 ;   segmentLengths_SDLAD[10] =  1.56;
segmentLengths_meanLAD[11] = 2.82 ;   segmentLengths_SDLAD[11] =  1.96;

 elementDiameters_meanLAD[0] 	   = 0.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[1] 	   = 9.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[2] 	   = 12.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[3] 	   = 17.7; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[4] 	   = 30.5; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[5] 	   = 66.2; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[6] 	   = 139; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[7] 	   = 308; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[8] 	   = 462; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[9] 	   = 714; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[10] 	   = 1573; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLAD[11] 	   = 3176; // micrometers
	
// LCX, table 3 and table 5.
 segmentElementRatio_meanLCX[0] = 0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[1] = 2.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[2] =       1.79; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[3] =       2.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[4] =      2.06; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[5] =       2.20; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[6] =       2.11; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[7] =   2.75; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[8] =     4.22; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[9] =   6.6; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[10] =    14.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_meanLCX[11] =      -1; // micrometers

 segmentElementRatio_SDLCX[0]      = 0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[1]      = 1.4; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[2]      =        0.95; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[3]      =      1.1; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[4]      =        1.2; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[5]      =        1.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[6]      =         1.1; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[7]      =       1.6; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[8]      =      2.4; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[9]      =     4.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[10]      =     4.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 segmentElementRatio_SDLCX[11]      =      -1;  // micrometers

 elementLengths_meanLCX[0] 	          = 0 ; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[1] 	          = 0.115 ; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[2] 	          =  0.136; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[3] 	          =    0.149; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[4] 	          =    0.405; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[5] 	          =    0.908; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[6] 	          =    1.83 ; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[7] 	          =   4.22; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[8] 	          =    6.98; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[9] 	          =   21; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[10] 	          =    49; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_meanLCX[11] 	          =         -1; // this is mm

segmentLengths_meanLCX[0] =  0.0 ; segmentLengths_SDLCX[0] =    0.0; 
segmentLengths_meanLCX[1] =  0.056 ; segmentLengths_SDLCX[1] =    0.038; 
segmentLengths_meanLCX[2] =  0.072 ; segmentLengths_SDLCX[2] =    0.045; 
segmentLengths_meanLCX[3] =  0.072 ; segmentLengths_SDLCX[3] =    0.049; 
segmentLengths_meanLCX[4] =  0.190 ; segmentLengths_SDLCX[4] =    0.097; 
segmentLengths_meanLCX[5] =  0.615 ; segmentLengths_SDLCX[5] =    0.508; 
segmentLengths_meanLCX[6] =  1.11 ; segmentLengths_SDLCX[6] =    0.983; 
segmentLengths_meanLCX[7] =  1.6 ; segmentLengths_SDLCX[7] =    1.33; 
segmentLengths_meanLCX[8] =  1.78 ; segmentLengths_SDLCX[8] =    1.46; 
segmentLengths_meanLCX[9] =  3.18 ; segmentLengths_SDLCX[9] =    2.41; 
segmentLengths_meanLCX[10] =  3.54 ; segmentLengths_SDLCX[10] =    2.0; 
segmentLengths_meanLCX[11] =  -1 ; segmentLengths_SDLCX[11] =    -1; 

 elementLengths_SDLCX[0] 		   = 0.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[1] 		   = 0.066; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[2] 		   =  0.088; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[3] 		   =  0.104; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[4] 		   =  0.170; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[5] 		   =      0.763; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[6] 		   =    1.34; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[7] 		   =     2.26; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[8] 		   =    3.92; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[9] 		   =   15.6; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[10] 		   =   0.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementLengths_SDLCX[11] 		   =      -1; // mm

 elementDiameters_meanLCX[0] 	   = 0.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[1] 	   = 9.0; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[2] 	   =       12.3; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[3] 	   =       17.7; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[4] 	   =    27.5; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[5] 	   =         73.2; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[6] 	   =       139; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[7] 	   =    279; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[8] 	   =      462; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[9] 	   =   961; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[10] 	   =  2603; // SE ratio and Diameter: micrometers, length mean and SD in mm
 elementDiameters_meanLCX[11] 	   =    -1; // micrometers

// change lengths to micrometers.
// SE ratio and Diameter: micrometers, length mean and SD in mm

double scaling;

// Comparision of data from Dodge for human to Kassab 1993 data seems to indicate that there is a scaling of 1.1 to 1.2 in all lengths.
// scaling = 1.1; // from Dodge et al.

for(temp_x = 0; temp_x < MAXORDERNUM; temp_x++){
	if(temp_x<11) scaling = 1.0; else scaling = 1.4;
	 elementDiameters_meanRCA[temp_x] 	   =    elementDiameters_meanRCA[temp_x]/1000.0;
	 elementDiameters_meanLAD[temp_x] 	   =    elementDiameters_meanLAD[temp_x]/1000.0;
	 elementDiameters_meanLCX[temp_x] 	   	   =    elementDiameters_meanLCX[temp_x]/1000.0;
	 
	segmentLengths_meanRCA[temp_x] = scaling * segmentLengths_meanRCA[temp_x];
	segmentLengths_meanLAD[temp_x] = scaling * segmentLengths_meanLAD[temp_x];
	segmentLengths_meanLCX[temp_x] = scaling * segmentLengths_meanLCX[temp_x];
}

