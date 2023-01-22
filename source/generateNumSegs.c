// generate the right branch number of segments. this is called when we are not at the end of an element yet.
	/* right node is where you attach a suitable new element. */
				do{
				U1 = genrand_real1(); order2 = -1;
					// w.r.t. notation in connectivityKassab.c: colum is order. row is order2 what has to be chosen.
					// pick the order of branch you are attaching here: 
					if(U1<cuprob_conn[0][order] ){order2 = 0; break; }  // for the first row.
				// lower bound of the interval is what you need.
				// for smaller trees: revise the maxOrderNum to maxorder+1.
				if(whichCase==3) maxOrderNum = MAXORDERNUM-1; else   maxOrderNum = MAXORDERNUM; 
				for(row=1; row< startOrder; row++) if(cuprob_conn[row-1][order]< U1 && U1 <= cuprob_conn[row][order] ) order2 = row;
				} while(order2 < 0 || order2 > order || order2 >= maxOrderNum ); // order2 > order means order2 can be equal to order (many elements). order2 >= order for testing.


			mean  							= segmentElementRatio_mean[order2];
			sigma 							= segmentElementRatio_SD[order2];						
			numberOfSegmentsInElement2 	= numSegsF(mean, sigma);

