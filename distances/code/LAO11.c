if(whichCase==2 && order == 11 && segLabel<1){ // ellipse.

	if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order])) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 32; tNode->y_cood = 50; tNode->z_cood = 2; tNode->onOff = 1;	
	}
else	if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order])-1) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 37; tNode->y_cood = 49; tNode->z_cood = 2; tNode->onOff = 1;	
	}
else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order])-2) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 42; tNode->y_cood = 47; tNode->z_cood = 2; tNode->onOff = 1;	
	}
else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order])-3) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 46; tNode->y_cood = 43.5; tNode->z_cood = 2; tNode->onOff = 1;	
	}
else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order])-4) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 50; tNode->y_cood = 40; tNode->z_cood = 2; tNode->onOff = 1;	
	}
else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order])-5) ){
		tNode->assignedOrNot = 1;  tNode->x_cood = 61.5; tNode->y_cood = 36; tNode->z_cood = 2; tNode->onOff = 1;
		t0 = atan2( (tNode->z_cood - LV_z_c) * LV_X_RADIUS, (tNode->x_cood - LV_x_c) * LV_Z_RADIUS );		
	}
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 6)){
  tNode->assignedOrNot = 1; tNode->x_cood = 61.995789; tNode->y_cood = 36.000000; tNode->z_cood = 3.172806;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 7)){
  tNode->assignedOrNot = 1; tNode->x_cood = 61.911180; tNode->y_cood = 36.000000; tNode->z_cood = 7.382513;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 8)){
  tNode->assignedOrNot = 1; tNode->x_cood = 61.798163; tNode->y_cood = 36.000000; tNode->z_cood = 10.106279;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 9)){
  tNode->assignedOrNot = 1; tNode->x_cood = 61.372430; tNode->y_cood = 36.000000; tNode->z_cood = 16.242961;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 10)){
  tNode->assignedOrNot = 1; tNode->x_cood = 60.913107; tNode->y_cood = 36.000000; tNode->z_cood = 20.671365;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 11)){
  tNode->assignedOrNot = 1; tNode->x_cood = 60.222874; tNode->y_cood = 36.000000; tNode->z_cood = 25.734643;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 12)){
  tNode->assignedOrNot = 1; tNode->x_cood = 59.754304; tNode->y_cood = 36.000000; tNode->z_cood = 28.573234;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 13)){
  tNode->assignedOrNot = 1; tNode->x_cood = 58.286678; tNode->y_cood = 36.000000; tNode->z_cood = 35.733463;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 14)){
  tNode->assignedOrNot = 1; tNode->x_cood = 56.815360; tNode->y_cood = 36.000000; tNode->z_cood = 41.335731;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 15)){
  tNode->assignedOrNot = 1; tNode->x_cood = 53.109283; tNode->y_cood = 36.000000; tNode->z_cood = 51.738773;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 16)){
  tNode->assignedOrNot = 1; tNode->x_cood = 52.937464; tNode->y_cood = 36.000000; tNode->z_cood = 52.132610;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 17)){
  tNode->assignedOrNot = 1; tNode->x_cood = 51.479460; tNode->y_cood = 36.000000; tNode->z_cood = 55.236351;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 18)){
  tNode->assignedOrNot = 1; tNode->x_cood = 50.112940; tNode->y_cood = 36.000000; tNode->z_cood = 57.801383;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 19)){
  tNode->assignedOrNot = 1; tNode->x_cood = 48.780128; tNode->y_cood = 36.000000; tNode->z_cood = 60.025797;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 20)){
  tNode->assignedOrNot = 1; tNode->x_cood = 47.343292; tNode->y_cood = 36.000000; tNode->z_cood = 62.152201;  tNode->onOff = 1; 
} 
 else if(numberOfSegmentsInElement==((int)(segmentElementRatio_mean[order]) - 21)){
  tNode->assignedOrNot = 1; tNode->x_cood = 45; tNode->y_cood = 36.000000; tNode->z_cood = 64;  tNode->onOff = 1; 
} 
else
{
tNode->assignedOrNot = 0;  tNode->x_cood = -1; tNode->y_cood = -1; tNode->z_cood = -1; tNode->onOff = 0;
}

}
