main = wanr2k.f90
OBJmain = wanr2k.o
FFLAGS = -O3 -ffree-line-length-0

 
FC = gfortran
LIB=-llapack -lblas

wanr2k : ${OBJmain} 
	${FC} ${FFLAGS} -o wanr2k wanr2k.o ${LIB}

${OBJmain} : ${@:.o=.f90} ${main}
	${FC} ${FFLAGS} -c ${@:.o=.f90}

