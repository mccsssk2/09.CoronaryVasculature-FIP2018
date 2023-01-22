for(colum=0; colum<12; colum++)
	for(row=0; row<12; row++){
		connRCA[row][colum] = 0.0;
		connLAD[row][colum] = 0.0;
		connLCX[row][colum] = 0.0;

		prob_connRCA[row][colum] = 0.0;
		prob_connLAD[row][colum] = 0.0;
		prob_connLCX[row][colum] = 0.0;		

		cuprob_connRCA[row][colum] = 0.0;
		cuprob_connLAD[row][colum] = 0.0;
		cuprob_connLCX[row][colum] = 0.0;
	}

// RCA connectivity matrix. Table 6, page H360 Morphometry of coronary arteries.
connRCA[0][1] = 2.75;  // row 0, column 1
connRCA[0][2] = 0.674; // row 0, column 2
connRCA[0][3] = 0.151; // 
connRCA[0][4] = 0.040; // 
connRCA[0][5] = 0; // 
connRCA[0][6] = 0; // 
connRCA[0][7] = 0; // 
connRCA[0][8] = 0; // 
connRCA[0][9] = 0; // 
connRCA[0][10] = 0; // 
connRCA[0][11] = 0; // 

connRCA[1][1] = 0.131; // row 1, column 1
connRCA[1][2] = 2.13; // row 1, column 2
connRCA[1][3] = 0.802; // 
connRCA[1][4] = 0.300; // 
connRCA[1][5] = 0.008; // 
connRCA[1][6] = 0.004; // 
connRCA[1][7] = 0; // 
connRCA[1][8] = 0; // 
connRCA[1][9] = 0; // 
connRCA[1][10] = 0; // 
connRCA[1][11] = 0; // 

connRCA[2][1] = 0.0;  // row 2, column 1
connRCA[2][2] = 0.080; // row 2, column 2
connRCA[2][3] = 2.15; // 
connRCA[2][4] = 0.700; // 
connRCA[2][5] = 0.159; // 
connRCA[2][6] = 0.020; // 
connRCA[2][7] = 0.037; // 
connRCA[2][8] = 0; // 
connRCA[2][9] = 0; // 
connRCA[2][10] = 0; // 
connRCA[2][11] = 0; // 

connRCA[3][1] = 0.0;  // row 3, column 1
connRCA[3][2] = 0.0;  // row 3, column 2
connRCA[3][3] = 0.070; // 
connRCA[3][4] = 2.12; // 
connRCA[3][5] = 0.688; // 
connRCA[3][6] = 0.314; // 
connRCA[3][7] = 0.244; // 
connRCA[3][8] = 0.143; // 
connRCA[3][9] = 0.059; // 
connRCA[3][10] = 0; // 
connRCA[3][11] = 0; // 

connRCA[4][1] = 0.0;  // row 4, column 1
connRCA[4][2] = 0.0;  // row 4, column 2
connRCA[4][3] = 0.0;  // 
connRCA[4][4] = 0.160; // 
connRCA[4][5] = 2.13; // 
connRCA[4][6] = 0.444; // 
connRCA[4][7] = 0.645; // 
connRCA[4][8] = 0.468; // 
connRCA[4][9] = 0.324; // 
connRCA[4][10] = 0; // 
connRCA[4][11] = 0; // 

connRCA[5][1] = 0.0;  // row 5, column 1
connRCA[5][2] = 0.0;  // row 5, column 2
connRCA[5][3] = 0.0;  // 
connRCA[5][4] = 0.0;  // 
connRCA[5][5] = 0.344; // 
connRCA[5][6] = 2.43; // 
connRCA[5][7] = 1.43; // 
connRCA[5][8] = 1.51; // 
connRCA[5][9] = 0.853; // 
connRCA[5][10] = 0.727; // 
connRCA[5][11] = 0; // 

connRCA[6][1] = 0.0;  // row 6, column 1
connRCA[6][2] = 0.0;  // row 6, column 2
connRCA[6][3] = 0.0;  // 
connRCA[6][4] = 0.0;  // 
connRCA[6][5] = 0.0;  // 
connRCA[6][6] = 0.167; // 
connRCA[6][7] = 2.14; // 
connRCA[6][8] = 1.37; // 
connRCA[6][9] = 1.68; // 
connRCA[6][10] = 1.80; // 
connRCA[6][11] = 1; // 

connRCA[7][1] = 0.0;  // row 7, column 1
connRCA[7][2] = 0.0;  // row 7 column 2
connRCA[7][3] = 0.0;  // 
connRCA[7][4] = 0.0;  // 
connRCA[7][5] = 0.0;  // 
connRCA[7][6] = 0.0;  // 
connRCA[7][7] = 0.130; // 
connRCA[7][8] = 2.32; // 
connRCA[7][9] = 1.23; // 
connRCA[7][10] = 2.80; // 
connRCA[7][11] = 11; //
 
connRCA[8][1] = 0.0;  // row 8, column 1
connRCA[8][2] = 0.0;  // row 8 column 2
connRCA[8][3] = 0.0;  // 
connRCA[8][4] = 0.0;  // 
connRCA[8][5] = 0.0;  // 
connRCA[8][6] = 0.0;  // 
connRCA[8][7] = 0.0;  // 
connRCA[8][8] = 0.099; // 
connRCA[8][9] = 2.62; // 
connRCA[8][10] = 2.00; // 
connRCA[8][11] = 8; // 

connRCA[9][1] = 0.0;  // row 9, column 1
connRCA[9][2] = 0.0;  // row 9, column 2
connRCA[9][3] = 0.0;  // 
connRCA[9][4] = 0.0;  // 
connRCA[9][5] = 0.0;  // 
connRCA[9][6] = 0.0;  // 
connRCA[9][7] = 0.0;  // 
connRCA[9][8] = 0.0;  // 
connRCA[9][9] = 0.059; // 
connRCA[9][10] = 2.40; // 
connRCA[9][11] = 5; // 

connRCA[10][1] = 0.0;  // row 10, column 1
connRCA[10][2] = 0.0;  // row 10, column 2
connRCA[10][3] = 0.0;  // 
connRCA[10][4] = 0.0;  // 
connRCA[10][5] = 0.0;  // 
connRCA[10][6] = 0.0;  // 
connRCA[10][7] = 0.0;  // 
connRCA[10][8] = 0.0;  // 
connRCA[10][9] = 0.0;  // 
connRCA[10][10] = 0.400; // 
connRCA[10][11] = 6; // 

// LAD connectivity matrix. Table 7, page H360 Morphometry of coronary arteries.
connLAD[0][1] = 3.18; // Row 0, Column 1 
connLAD[0][2] = 0.675; // Row 0, Column 2
connLAD[0][3] = 0.148; // 
connLAD[0][4] = 0; // 
connLAD[0][5] = 0; // 
connLAD[0][6] = 0; // 
connLAD[0][7] = 0; // 
connLAD[0][8] = 0; // 
connLAD[0][9] = 0; // 
connLAD[0][10] = 0; // 
connLAD[0][11] = 0; // 

connLAD[1][1] =0.144 ; // Row 1, Column 1
connLAD[1][2] = 2.04; //  Row 1, Column 2
connLAD[1][3] = 0.630; // 
connLAD[1][4] = 0.071; // 
connLAD[1][5] = 0; // 
connLAD[1][6] = 0; // 
connLAD[1][7] = 0; // 
connLAD[1][8] = 0; // 
connLAD[1][9] = 0; // 
connLAD[1][10] = 0; // 
connLAD[1][11] = 0; // 

connLAD[2][1] = 0.0;  // Row 2, Column 1
connLAD[2][2] =0.094 ; // Row 2, Column 2
connLAD[2][3] = 2.24; // 
connLAD[2][4] = 1.50; // 
connLAD[2][5] = 0.063; // 
connLAD[2][6] = 0.094; // 
connLAD[2][7] = 0.023; // 
connLAD[2][8] = 0; // 
connLAD[2][9] = 0; // 
connLAD[2][10] = 0; // 
connLAD[2][11] = 0; // 

connLAD[3][1] = 0.0;  // Row 3, Column 1
connLAD[3][2] = 0.0;  // Row 3, Column 2
connLAD[3][3] = 0.074; // 
connLAD[3][4] = 2.14; // 
connLAD[3][5] = 0.381; // 
connLAD[3][6] = 0.098; // 
connLAD[3][7] = 0.120; // 
connLAD[3][8] =0.092 ; // 
connLAD[3][9] = 0.030; // 
connLAD[3][10] = 0; // 
connLAD[3][11] = 0; // 

connLAD[4][1] = 0.0;  // Row 4, Column 1
connLAD[4][2] = 0.0;  // Row 4, Column 2
connLAD[4][3] = 0.0;  // 
connLAD[4][4] = 0.143 ; // 
connLAD[4][5] = 2.25 ; // 
connLAD[4][6] = 0.425; // 
connLAD[4][7] = 0.380; // 
connLAD[4][8] = 0.428; // 
connLAD[4][9] = 0.303; // 
connLAD[4][10] = 0.167; // 
connLAD[4][11] = 0; // 

connLAD[5][1] = 0.0;  // Row 5, Column 1
connLAD[5][2] = 0.0;  // Row 5, Column 2
connLAD[5][3] = 0.0;  // 
connLAD[5][4] = 0.0;  // 
connLAD[5][5] =0.238 ; // 
connLAD[5][6] = 0.250; // 
connLAD[5][7] = 1.91; // 
connLAD[5][8] = 1.56; // 
connLAD[5][9] = 1.36; // 
connLAD[5][10] = 0.667 ; // 
connLAD[5][11] = 0; // 

connLAD[6][1] = 0.0;  // Row 6, Column 1
connLAD[6][2] = 0.0;  // Row 6, Column 2
connLAD[6][3] = 0.0;  // 
connLAD[6][4] = 0.0;  // 
connLAD[6][5] = 0.0;  // 
connLAD[6][6] = 0.155; // 
connLAD[6][7] = 2.50; // 
connLAD[6][8] = 1.58; // 
connLAD[6][9] = 1.48; // 
connLAD[6][10] = 1.17; // 
connLAD[6][11] = 0; // 

connLAD[7][1] = 0.0;  // Row 7, Column 1
connLAD[7][2] = 0.0;  // Row 7, Column 2
connLAD[7][3] = 0.0;  // 
connLAD[7][4] = 0.0;  // 
connLAD[7][5] = 0.0;  // 
connLAD[7][6] = 0.0;  // 
connLAD[7][7] = 0.116; // 
connLAD[7][8] = 2.09; // 
connLAD[7][9] = 1.30; // 
connLAD[7][10] = 2.00; // 
connLAD[7][11] = 2; // 

connLAD[8][1] = 0.0;  // Row 8, Column 1
connLAD[8][2] = 0.0;  // Row 8, Column 2
connLAD[8][3] = 0.0;  // 
connLAD[8][4] = 0.0;  // 
connLAD[8][5] = 0.0;  // 
connLAD[8][6] = 0.0;  // 
connLAD[8][7] = 0.0;  // 
connLAD[8][8] = 0.061; // 
connLAD[8][9] = 2.50; // 
connLAD[8][10] = 2.50; // 
connLAD[8][11] = 3; // 

connLAD[9][1] = 0.0;  // Row 9, Column 1
connLAD[9][2] = 0.0;  // Row 9, Column 2
connLAD[9][3] = 0.0;  // 
connLAD[9][4] = 0.0;  // 
connLAD[9][5] = 0.0;  // 
connLAD[9][6] = 0.0;  // 
connLAD[9][7] = 0.0;  // 
connLAD[9][8] = 0.0;  // 
connLAD[9][9] = 0.121; // 
connLAD[9][10] = 3.33; // 
connLAD[9][11] = 8; // 

connLAD[10][1] = 0.0;  // Row 10, Column 1
connLAD[10][2] = 0.0;  // Row 10, Column 2
connLAD[10][3] = 0.0;  // 
connLAD[10][4] = 0.0;  // 
connLAD[10][5] = 0.0;  // 
connLAD[10][6] = 0.0;  // 
connLAD[10][7] = 0.0;  // 
connLAD[10][8] = 0.0;  // 
connLAD[10][9] = 0.0;  // 
connLAD[10][10] = 0.100; // 
connLAD[10][11] = 5; // 

connLAD[11][1] = 0.0;  // Row 11, Column 1
connLAD[11][2] = 0.0;  // Row 11, Column 2
connLAD[11][3] = 0.0;  // 
connLAD[11][4] = 0.0;  // 
connLAD[11][5] = 0.0;  // 
connLAD[11][6] = 0.0;  // 
connLAD[11][7] = 0.0;  // 
connLAD[11][8] = 0.0;  // 
connLAD[11][9] = 0.0;  // 
connLAD[11][10] = 0.0;  // 
connLAD[11][11] = 0; // 


// LCX connectivity matrix. Table 8, page H361 Morphometry of coronary arteries.
connLCX[0][1] = 3.18; // row 0, column 1
connLCX[0][2] = 0.675; // row 0, column 2
connLCX[0][3] = 0.148;  //
connLCX[0][4] = 0; //
connLCX[0][5] = 0; //
connLCX[0][6] = 0; //
connLCX[0][7] = 0; //
connLCX[0][8] = 0; //
connLCX[0][9] = 0; //
connLCX[0][10] = 0; //

connLCX[1][1] = 0.144; // row 1, column 1 
connLCX[1][2] = 2.04; // row 1, column 2
connLCX[1][3] = 0.630; //
connLCX[1][4] = 0.071; //
connLCX[1][5] = 0; //
connLCX[1][6] = 0; //
connLCX[1][7] = 0; //
connLCX[1][8] = 0; //
connLCX[1][9] = 0; //
connLCX[1][10] = 0; //

connLCX[2][1] = 0.0;  // row 2, column 1 
connLCX[2][2] = 0.094; // row 2, column 2 
connLCX[2][3] = 2.24; //
connLCX[2][4] = 1.50; //
connLCX[2][5] = 0.150; //
connLCX[2][6] = 0; //
connLCX[2][7] = 0; //
connLCX[2][8] = 0; //
connLCX[2][9] = 0; //
connLCX[2][10] = 0; //

connLCX[3][1] = 0.0;  // row 3, column 1 
connLCX[3][2] = 0.0;  // row 3, column 2 
connLCX[3][3] = 0.074; //
connLCX[3][4] = 2.14; //
connLCX[3][5] = 0.150; //
connLCX[3][6] = 0.025; //
connLCX[3][7] = 0.011; //
connLCX[3][8] = 0; //
connLCX[3][9] = 0; //
connLCX[3][10] = 0; //

connLCX[4][1] = 0.0;  // row 4, column 1 
connLCX[4][2] = 0.0;  // row 4, column 2 
connLCX[4][3] = 0.0;  //
connLCX[4][4] = 0.143; //
connLCX[4][5] = 2.85; //
connLCX[4][6] = 0.385; //
connLCX[4][7] = 0.179; //
connLCX[4][8] = 0.109; //
connLCX[4][9] = 0; //
connLCX[4][10] = 0; //

connLCX[5][1] = 0.0;  // row 5, column 1 
connLCX[5][2] = 0.0;  // row 5, column 2 
connLCX[5][3] = 0.0;  //
connLCX[5][4] = 0.0;  //
connLCX[5][5] = 0.100; //
connLCX[5][6] = 2.51; //
connLCX[5][7] = 1.13; //
connLCX[5][8] = 0.956; //
connLCX[5][9] = 0.444; //
connLCX[5][10] = 0; //

connLCX[6][1] = 0.0;  // row 6, column 1 
connLCX[6][2] = 0.0;  // row 6, column 2 
connLCX[6][3] = 0.0;  //
connLCX[6][4] = 0.0;  //
connLCX[6][5] = 0.0;  //
connLCX[6][6] = 0.213; //
connLCX[6][7] = 2.33; //
connLCX[6][8] = 2.02; //
connLCX[6][9] = 1.78; //
connLCX[6][10] = 2; //

connLCX[7][1] = 0.0;  // row 7, column 1 
connLCX[7][2] = 0.0;  // row 7, column 2 
connLCX[7][3] = 0.0;  //
connLCX[7][4] = 0.0;  //
connLCX[7][5] = 0.0;  //
connLCX[7][6] = 0.0;  //
connLCX[7][7] = 0.168; //
connLCX[7][8] = 2.04; //
connLCX[7][9] = 2.00; //
connLCX[7][10] = 4; //

connLCX[8][1] = 0.0;  // row 8, column 1 
connLCX[8][2] = 0.0;  // row 8, column 2 
connLCX[8][3] = 0.0;  //
connLCX[8][4] = 0.0;  //
connLCX[8][5] = 0.0;  //
connLCX[8][6] = 0.0;  //
connLCX[8][7] = 0.0;  //
connLCX[8][8] = 0.196; //
connLCX[8][9] = 4.00; //
connLCX[8][10] = 4; //

connLCX[9][1] = 0.0;  // row 9, column 1 
connLCX[9][2] = 0.0;  // row 9, column 2 
connLCX[9][3] = 0.0;  //
connLCX[9][4] = 0.0;  //
connLCX[9][5] = 0.0;  //
connLCX[9][6] = 0.0;  //
connLCX[9][7] = 0.0;  //
connLCX[9][8] = 0.0; //
connLCX[9][9] = 0.111; //
connLCX[9][10] = 8; //

connLCX[10][1] = 0.0;  // row 10, column 1 
connLCX[10][2] = 0.0;  // row 10, column 2 
connLCX[10][3] = 0.0;  //
connLCX[10][4] = 0.0;  //
connLCX[10][5] = 0.0;  //
connLCX[10][6] = 0.0;  //
connLCX[10][7] = 0.0;  // 
connLCX[10][8] = 0.0;  //
connLCX[10][9] = 0.0;  //
connLCX[10][10] = 0; //

// now you want probabilities. Take sum of each column, divide entries of that column by that sum. Go to next column and do the same.
for(colum=0; colum<12; colum++){
// totals	
	total_col1 = total_col2 = total_col3 = 0.0;
	for(row=0; row<12; row++){
		total_col1 =  total_col1 	+ connRCA[row][colum];
		total_col2 = total_col2 	+ connLAD[row][colum];
		total_col3 = total_col3 	+ connLCX[row][colum]; 
	}
// normalised
	for(row=0; row<12; row++){
		if(total_col1>0.0){
			prob_connRCA[row][colum] = connRCA[row][colum]/total_col1;
		}

		if(total_col2>0.0){
			prob_connLAD[row][colum] = connLAD[row][colum]/total_col2;
		}
		
		if(total_col3>0.0){
			prob_connLCX[row][colum] = connLCX[row][colum]/total_col3;
		}
	}	
}

// now cumulate the probability.
for(colum=0; colum<12; colum++)	
	for(row2=0; row2<12; row2++){
	total_col1 = total_col2 = total_col3 = 0.0;
		for(row=0; row<=row2; row++){
			total_col1 = total_col1   + prob_connRCA[row][colum];
			total_col2 = total_col2 	+ prob_connLAD[row][colum];
			total_col3 = total_col3 	+ prob_connLCX[row][colum]; 
		}
		cuprob_connRCA[row2][colum] = total_col1;
		cuprob_connLAD[row2][colum] = total_col2;
		cuprob_connLCX[row2][colum] = total_col3;		
	}


	// assymmetry matrix.
// Koimovitz and Kassab paper, Figure 1A.
amx[0] = amx[1] = amx[2] = amx[3] = amx[4] = 1.0;
amx[5] = 0.81; amx[6] = 0.875; amx[7] = 0.925;
amx[8] = amx[9] = amx[10] = amx[11] = 0.96;
 	
	
// print it out.
 if(debug){
printf("RCA probability matrix:\n");
for(row=0; row<12; row++){
for(colum=0; colum<12; colum++)		
		printf("%2.3f\t", cuprob_connRCA[row][colum] * 100);
		printf("\n");
}

printf("LAD probability matrix:\n");
for(row=0; row<12; row++){
for(colum=0; colum<12; colum++)		
		printf("%2.3f\t", cuprob_connLAD[row][colum] * 100) ;
		printf("\n");
}

printf("LCX probability matrix:\n");
for(row=0; row<12; row++){
for(colum=0; colum<12; colum++)		
		printf("%2.3f\t", prob_connLCX[row][colum] * 100);
		printf(" Row: %d ", row);
		printf("\n");
}


 } // end of if debug.
