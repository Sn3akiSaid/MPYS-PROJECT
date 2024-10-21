main = interpolate.f90
OBJmain = interpolate.o
FFLAGS = -O3 -ffree-line-length-0

 
FC = gfortran
LIB=-llapack -lblas

interpolate : ${OBJmain} 
	${FC} ${FFLAGS} -o interpolate interpolate.o ${LIB}

${OBJmain} : ${@:.o=.f90} ${main}
	${FC} ${FFLAGS} -c ${@:.o=.f90}

