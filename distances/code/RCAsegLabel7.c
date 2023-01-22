if(whichCase==1 && order == 10 && segLabel == 7){
	if(numberOfSegmentsInElement==(9)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 4.4; tNode->y_cood = 27; tNode->z_cood = 24; tNode->onOff = 1; 
		segLabel=7; 
		if(debug) 
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-1)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 4.8; tNode->y_cood = 24; tNode->z_cood = 23; tNode->onOff = 1; 
		segLabel=7;
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-2)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 5.4; tNode->y_cood = 22; tNode->z_cood = 22; tNode->onOff = 1; 
		segLabel=7;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-3)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 6.5; tNode->y_cood = 19; tNode->z_cood = 22; tNode->onOff = 1; 
		segLabel=7;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-4)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 7.5; tNode->y_cood = 17; tNode->z_cood = 20; tNode->onOff = 1; 
		segLabel=7;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-5)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 9.3; tNode->y_cood = 14; tNode->z_cood = 19; tNode->onOff = 1; 
		segLabel=7;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-6)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 10.7; tNode->y_cood = 12; tNode->z_cood = 18; tNode->onOff = 1;
		segLabel=7;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-7)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 13; tNode->y_cood = 10.2; tNode->z_cood = 16; tNode->onOff = 1;
		segLabel=7;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-8)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 15; tNode->y_cood = 8.5; tNode->z_cood = 16; tNode->onOff = 1;
		segLabel=0;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}

}
