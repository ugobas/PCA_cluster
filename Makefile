PRG= PCA_cluster
CC = gcc

CFLAGS=-g
#CFLAGS=-O2
SRC= PCA_cluster.c expmax.c expmax_aux.c expmax_order.c cluster_score.c \
     	Initialize_expmax.c random3.c	\
     	choldc.c average_linkage.c allocate.c \
	diagonalize.c  jacobi.c tqli.c \
	tred2.c eigsrt.c nrutil.c pythag.c
OBJ= PCA_cluster.o expmax.o expmax_order.o expmax_aux.o cluster_score.o \
        Initialize_expmax.o random3.o	\
     	choldc.o average_linkage.o allocate.o \
	diagonalize.o  jacobi.o tqli.o\
	tred2.o eigsrt.o nrutil.o pythag.o


$(PRG): $(OBJ)
	$(CC) $(CFLAGS) -o $(PRG) $(OBJ) -lm -static
	rm -f *.o
	
	echo "Executable file generated:" $(PRG)
	

$(OBJ) : $(HEAD) $(SRC)
	$(CC) $(CFLAGS) -c  $(SRC)

