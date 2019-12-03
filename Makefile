FC := gfortran #This is the variable for the fortran compiler to be used
FFLAGS := -o #These are the flags that we use while compilation
RM := rm -f
OBJ := sod_shock.o
PROG := sod_shock
OUTDIR := fortran_op

.PHONY: all clean run


all: $(PROG)

clean: 
	$(RM) $(PROG) *.o $(OUTDIR)/*

run: $(PROG)
	./$(PROG) || { echo "$(PROG) failed"; exit 1; };

%.o: %.f90
	$(FC) -c -o $@ $^

$(PROG): $(OBJ) $(OUTDIR)
	$(FC) $(FFLAGS) $@ $<


$(OUTDIR):
	mkdir -p $@
