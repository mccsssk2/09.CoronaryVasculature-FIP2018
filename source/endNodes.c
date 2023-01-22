			if(whichCase==1 && order == 10 && segLabel==1 ){
					tNode->assignedOrNot = 1;  tNode->x_cood = 44; tNode->y_cood = 57.75; tNode->z_cood = 22; tNode->onOff = 1;
					segLabel=10000;	// dont call this after this node.
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			}
			
			if(whichCase==1 && order==10 && segLabel==3){
					tNode->assignedOrNot = 1;  tNode->x_cood = 35.5; tNode->y_cood = 68.5; tNode->z_cood = 43; tNode->onOff = 1;
					segLabel=100000;
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			}

			if(whichCase==1 && order == 10 && segLabel == 5){
					tNode->assignedOrNot = 1;  tNode->x_cood = 3.2; tNode->y_cood = 38; tNode->z_cood = 11; tNode->onOff = 1;
					segLabel=10000000;
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			}

			if(whichCase==1 && order == 10 && segLabel == 7){
					tNode->assignedOrNot = 1;  tNode->x_cood = 17; tNode->y_cood = 6.75; tNode->z_cood = 15; tNode->onOff = 1;
					segLabel=100000;
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			} 
			
			if(whichCase==2 && order==10 && segLabel==2){
					tNode->assignedOrNot = 1;  tNode->x_cood = 22.5; tNode->y_cood = 45.5; tNode->z_cood = 26; tNode->onOff = 1;
					segLabel=0;	// dont call this after this node.
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			} 
			
			if(whichCase==2 && order == 10 && segLabel == 4){
					tNode->assignedOrNot = 1;  tNode->x_cood = 55; tNode->y_cood = 51.5; tNode->z_cood = 33; tNode->onOff = 1;
					segLabel=0;	// dont call this after this node.
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			}  
			
			if(whichCase==2 && order == 10 && segLabel == 6){
					tNode->assignedOrNot = 1;  tNode->x_cood = 44; tNode->y_cood = 61; tNode->z_cood = 45; tNode->onOff = 1;
					segLabel=0;	// dont call this after this node.
					if(debug)	
					printf("%d \n", numberOfSegmentsInElement);	
			}			

