if(whichCase==2 && order == 10 && segLabel == 6){

	if(numberOfSegmentsInElement==(9)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 50.5; tNode->y_cood = 40; tNode->z_cood = 54; tNode->onOff = 1; 
		segLabel=6;
		if(debug)
		printf("%d \n", numberOfSegmentsInElement);
	}
	else if(numberOfSegmentsInElement==(9-1)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 48.5; tNode->y_cood = 42; tNode->z_cood = 55; tNode->onOff = 1; 
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-2)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 48.2; tNode->y_cood = 44.2; tNode->z_cood = 53; tNode->onOff = 1; 
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-3)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 47; tNode->y_cood = 47; tNode->z_cood = 51; tNode->onOff = 1; 
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-4)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 46; tNode->y_cood = 50; tNode->z_cood = 50.5; tNode->onOff = 1; 
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-5)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 44.5; tNode->y_cood = 53; tNode->z_cood = 50.5; tNode->onOff = 1; 
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-6)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 45; tNode->y_cood = 55; tNode->z_cood = 49; tNode->onOff = 1; 
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-7)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 45.25; tNode->y_cood = 57; tNode->z_cood = 47; tNode->onOff = 1;
		segLabel=6;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-8)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 45.25; tNode->y_cood = 59; tNode->z_cood = 45; tNode->onOff = 1;
		segLabel=6;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}

}
