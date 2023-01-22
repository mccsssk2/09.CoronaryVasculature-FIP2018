
vasc: srkVasculature.o
	gcc -o vasc srkVasculature.o -g3 -lm
#	./vasc 4

viz:	vizVasc.o
	gcc -o viz vizVasc.c -lm -g3 -lm

dis: srkvascDistance.o
	gcc -o mes srkvascDistance.o -g -lm

expt: srkvascExpts.o
	gcc -o expt srkvascExpts.o -g -lm

veryclean:
	rm -rf *.o *.vts *.dat geom srk3d *.vtk *.info *.bin fibs *.xyz *.txt *.data NodeData.data.*

