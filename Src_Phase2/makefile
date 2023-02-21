FC = ifort
LDFLAGS = -O -fPIC -vec-report0
LBLAS = -L/mnt/n7fs/ens/tp_guivarch/opt/OpenBLAS/ -lopenblas
OBJS  = main.o subspace_iter.o
MODS  = m_subspace_iter.o

main:  $(MODS) $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $^ $(LBLAS)

%.o: %.f90
	$(FC) $(LDFLAGS) -c $< $(LBLAS)

clean:
	(rm -f *.o *.mod main)
