% equations for ellipsoid.
%
clear all;
clear all;
close all;
close all;

 LV_X_RADIUS= 30;
 LV_Y_RADIUS= 30;
 LV_Z_RADIUS= 70;
 LV_THICKNESS= 12;

 RV_X_RADIUS= 30;
 RV_Y_RADIUS= 51;
 RV_Z_RADIUS= 60;
 RV_THICKNESS= 6;

 MARGIN=2.0;

	LV_x_c = MARGIN + LV_X_RADIUS;
	LV_y_c = MARGIN + LV_Y_RADIUS;
	LV_z_c = MARGIN;


fractionforOrder=1.0;

a = LV_X_RADIUS-fractionforOrder*LV_THICKNESS;  b = LV_Y_RADIUS-fractionforOrder*LV_THICKNESS; c = LV_Z_RADIUS-fractionforOrder*LV_THICKNESS;

t = linspace(0,pi,22);
p = linspace(0,pi,22);
[T,P] = meshgrid(t,p);

% x = a*cos(T).*cos(P);
% y = b*cos(T).*sin(P);
% z = c*sin(T);

x = a * cos(P).*cos(T) + LV_x_c; y = b * cos(P).*sin(T) + LV_y_c; z  = c*sin(P) + LV_z_c;	
figure,surf(x,y,z);
shading interp;
zlabel("height");
