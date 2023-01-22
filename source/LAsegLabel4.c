if(whichCase==2 && order == 10 && segLabel == 4){

	if(numberOfSegmentsInElement==(9)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 60.75; tNode->y_cood = 38; tNode->z_cood = 9; tNode->onOff = 1; 
		segLabel=4;
		if(debug)
		printf("%d \n", numberOfSegmentsInElement);
	}
	else if(numberOfSegmentsInElement==(9-1)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 60.75; tNode->y_cood = 42; tNode->z_cood = 11; tNode->onOff = 1; 
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-2)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 60.2; tNode->y_cood = 44; tNode->z_cood = 14; tNode->onOff = 1; 
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-3)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 59.75; tNode->y_cood = 46; tNode->z_cood = 17; tNode->onOff = 1; 
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-4)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 59; tNode->y_cood = 46; tNode->z_cood = 21; tNode->onOff = 1; 
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-5)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 58.5; tNode->y_cood = 47; tNode->z_cood = 23; tNode->onOff = 1; 
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-6)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 57.6; tNode->y_cood = 48; tNode->z_cood = 26; tNode->onOff = 1; 
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-7)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 56.5; tNode->y_cood = 50; tNode->z_cood = 28; tNode->onOff = 1;
		segLabel=4;		
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
	}
	else if(numberOfSegmentsInElement==(9-8)    ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 55.5; tNode->y_cood = 52; tNode->z_cood = 31; tNode->onOff = 1;
		segLabel=4;	// dont call this after this node.
		if(debug)	
		printf("%d \n", numberOfSegmentsInElement);	
}

}
