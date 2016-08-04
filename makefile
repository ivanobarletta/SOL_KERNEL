MAIN  = SOL_KERNEL4_ORCA2
FC    = mpiifort

#HDF5 = /users/home/opt-intel_2015.3.187/hdf5/hdf5-1.8.15-patch1
INC_PATH1 = -I/users/home/sco116/petsc/petsc-3.7.1/linux-gnu-intel/include/
INC_PATH2 = -I/users/home/sco116/petsc/petsc-3.7.1/include/
LIB_PATH1 = -L/users/home/sco116/petsc/petsc-3.7.1/linux-gnu-intel/lib
LIB_PATH2 = -L/users/home/sco116/petsc/petsc-3.7.1/

NCDF_INC = -I${NETCDF}/include
NCDF_LIB = -L${NETCDF}/lib -lnetcdff -lnetcdf
HDF5_INC = -I${HDF5}/include
HDF5_LIB = -L${HDF5}/lib -lhdf5_hl -lhdf5
XIOS_INC = -I${XIOS}inc
XIOS_LIB = -L${XIOS}lib -lxios

INC_PATH  = $(INC_PATH1) $(INC_PATH2) $(NCDF_INC) $(HDF5_INC) $(XIOS_INC)
LIB_PATH  = $(LIB_PATH1) $(NCDF_LIB) $(HDF5_LIB) $(XIOS_LIB)

#FCFLAGS = -r8 -O3 -xHost -fp-model source -traceback -g -warn -check $(NCDF_INC) $(HDF5_INC) $(XIOS_INC)
#LDFLAGS = $(FCFLAGS) $(NCDF_LIB) $(HDF5_LIB) $(XIOS_LIB) -lstdc++ -lz -lgpfs -lcurl  

FCFLAGS = -r8 -O3 -xHost -fp-model source -traceback -g -warn -check $(INC_PATH)
LDFLAGS = $(FCFLAGS) $(LIB_PATH) -lstdc++ -lz -lgpfs -lcurl -lpetsc 

OBJECTS = nc4interface.o wrk_nemo.o utils.o sol_kernel4.o 

SOURCES = $(OBJECTS:.o,.F90)

$(MAIN): $(OBJECTS)
	$(FC) *.o -o $(MAIN) $(LDFLAGS)

%.o: %.F90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm *.o *.mod

 
