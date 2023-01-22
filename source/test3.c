// test the orientation code.
#include<stdio.h>
#include<math.h>

int main(){

double p0_x, p0_y, p0_z, p1_x, p1_y; // this is just mpx, mpy, tx, ty
// segment 2, all other segments except for segmentj1 and segmentj2.
double p2_x, p2_y, p3_x, p3_y, p1_z, p2_z, p3_z;
double q1_x, q1_y, q1_z, q2_x, q2_y, q2_z, q3_x, q3_y, q3_z;

int intersected;
int o1, o2, o3, o4;
	
	// the two segments are (p1,q1) and (p2,q2)
 p1_x	= 0.0;
 p1_y	= 50.0;
 p1_z	= 0.0;
 q1_x	= 100.0;
 q1_y	= 50.0;
 q1_z	= 0.0;
 p2_x	= 0.0;
 p2_y	= 49.99999999999999;
 p2_z	= 0.0;
 q2_x	= 0.0;
 q2_y	= 100.0;
 q2_z	= 0.0;

 // the four triplets are: (p1, q1, p2), (p1, q1, q2), (p2, q2, p1), (p2, q2, q1)
// below is the Z directional single component of the vectors (q1 - p1) x (q1 - p2) etc.
	if( (q1_x - p1_x)*(q1_y - p2_y) - (q1_x - p2_x)*(q1_y - p1_y) > 0 ) o1 = 1; else o1 = -1;
	if( (q1_x - p1_x)*(q1_y - q2_y) - (q1_x - q2_x)*(q1_y - p1_y) > 0 ) o2 = 1; else o2 = -1;
	if( (q2_x- p2_x)*(q2_y - p1_y) - (q2_x - p1_x)*(q2_y - p2_y) > 0 ) o3 = 1; else o3 = -1;
	if( (q2_x- p2_x)*(q2_y - q1_y) - (q2_x - q1_x)*(q2_y - p2_y) > 0 ) o4 = 1; else o4 = -1;
	printf("%d %d %d %d\n", o1, o2, o3, o4);
	if( (o1/o2<0)&&(o3/o4<0) ) intersected++;

printf("%d\n", intersected);

	return 0;
}
