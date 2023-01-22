% number of nodes at each order based on Kassab Tables 5 and 9.

clear all;
clear all;
close all;
close all;

se = load('SERatios.EXPT');
se = flipud(se);
numElements = load('numVesselElements.EXPT');

RCAnumnodes = se(:,2).*numElements(:,2);


lcalcxse1 = se(:,4);
lcalcxne1 = numElements(:,4);

LADnumnodes = lcalcxse1 .* lcalcxne1;

format longG;
RCAcs = cumsum(RCAnumnodes);
LADcs = cumsum(LADnumnodes);


RCAcs - 0.2*RCAcs;
RCAcs + 0.2*RCAcs;

LADcs - 0.2*LADcs
LADcs + 0.2*LADcs

rca1 = RCAnumnodes(4:end);
la1 = LADnumnodes(4:end);

rca1 = rca1/rca1(1);
la1 = la1/la1(1);



