main = interpolate_xy.f90
OBJmain = interpolatexy.o
FFLAGS = -O3 -ffree-line-length-0

 
FC = gfortran
LIB=-llapack -lblas

interpolate : ${OBJmain} 
	${FC} ${FFLAGS} -o interpolate_xy interpolate_xy.o ${LIB}

${OBJmain} : ${@:.o=.f90} ${main}
	${FC} ${FFLAGS} -c ${@:.o=.f90}