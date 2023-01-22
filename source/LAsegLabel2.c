if(whichCase==2 && order==10 && segLabel==2){

	if(numberOfSegmentsInElement==(9) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 36; tNode->y_cood = 49; tNode->z_cood = 5.0; tNode->onOff = 1; 
		segLabel=2;
		if(debug)
		printf("%d \n", numberOfSegmentsInElement);
	}
	else if(numberOfSegmentsInElement==(9-1)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 36; tNode->y_cood = 49; tNode->z_cood = 8.0; tNode->onOff = 1; 
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-2)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 34; tNode->y_cood = 49.5; tNode->z_cood = 11.0; tNode->onOff = 1; 
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-3)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 33; tNode->y_cood = 49.5; tNode->z_cood = 13.0; tNode->onOff = 1; 
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-4)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 30; tNode->y_cood = 49.5; tNode->z_cood = 15; tNode->onOff = 1; 
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-5)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 28; tNode->y_cood = 49; tNode->z_cood = 18; tNode->onOff = 1; 
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-6)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 27; tNode->y_cood = 48.75; tNode->z_cood = 20.0; tNode->onOff = 1; 
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-7)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 25; tNode->y_cood = 47.5; tNode->z_cood = 22.0; tNode->onOff = 1;
		segLabel=2;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-8)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 23; tNode->y_cood = 45.7; tNode->z_cood = 24; tNode->onOff = 1;
		segLabel=2;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
}
