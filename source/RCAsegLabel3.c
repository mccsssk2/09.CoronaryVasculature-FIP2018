if(whichCase==1 && order==10 && segLabel==3){
	
	if(numberOfSegmentsInElement==(9)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 49; tNode->y_cood = 73; tNode->z_cood = 11; tNode->onOff = 1; 
		segLabel=3;
		if(debug)
		printf("%d \n", numberOfSegmentsInElement);
	}
	else if(numberOfSegmentsInElement==(9-1)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 48.25; tNode->y_cood = 73; tNode->z_cood = 16; tNode->onOff = 1; 
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-2)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 47.5; tNode->y_cood = 73; tNode->z_cood = 19; tNode->onOff = 1; 
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-3)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 46.5; tNode->y_cood = 73; tNode->z_cood = 22; tNode->onOff = 1; 
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-4)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 44.5; tNode->y_cood = 74; tNode->z_cood = 25; tNode->onOff = 1; 
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-5)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 42.5; tNode->y_cood = 73.5; tNode->z_cood = 30; tNode->onOff = 1; 
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-6)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 42; tNode->y_cood = 73; tNode->z_cood = 32; tNode->onOff = 1; 
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-7)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 40.2; tNode->y_cood = 72.2; tNode->z_cood = 35; tNode->onOff = 1;
		segLabel=3;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-8)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 38; tNode->y_cood = 70; tNode->z_cood = 40; tNode->onOff = 1;
		segLabel=3;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}

} // end of RCA segLabel 3.
