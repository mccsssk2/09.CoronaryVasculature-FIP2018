if(whichCase==1 && order == 10 && segLabel == 5){

	if(numberOfSegmentsInElement==(9)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 12; tNode->y_cood = 63; tNode->z_cood = 26; tNode->onOff = 1; 
		segLabel=5;
		if(debug)
		printf("%d \n", numberOfSegmentsInElement);
	}
	else if(numberOfSegmentsInElement==(9-1)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 11; tNode->y_cood = 62; tNode->z_cood = 24; tNode->onOff = 1; 
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-2)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 9.5; tNode->y_cood = 60; tNode->z_cood = 22; tNode->onOff = 1; 
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-3)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 7.75; tNode->y_cood = 57; tNode->z_cood = 20; tNode->onOff = 1; 
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-4)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 6; tNode->y_cood = 53; tNode->z_cood = 18; tNode->onOff = 1; 
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-5)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 5; tNode->y_cood = 50; tNode->z_cood = 16; tNode->onOff = 1; 
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-6)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 4.5; tNode->y_cood = 48; tNode->z_cood = 14; tNode->onOff = 1; 
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-7)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 3.75; tNode->y_cood = 45; tNode->z_cood = 12; tNode->onOff = 1;
		segLabel=5;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-8)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 3.2; tNode->y_cood = 40; tNode->z_cood = 11; tNode->onOff = 1;
		segLabel=5;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	
}
