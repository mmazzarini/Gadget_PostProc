#HDF5_DISABLE_VERSION_CHECK=1

$(info $)
$(info $)
$(info $)
$(info $)
$(info $   +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + $)
$(info $)
$(info $)
$(info $    Gadget_PostProc GNU Makefile by M. Mazzarini (ARI, Heidelberg, 2017)     $)
$(info $)
$(info $)
$(info $          Thanks to Manuel Arca Sedda for helping write the Makefile	     $) 					     	
$(info $)
$(info $)
$(info $   +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + $)
$(info $)
$(info $)
$(info $)
$(info $)

CC = g++ -std=c++11
FFLAGS += -O3 -W -Wall 

EXEC = GadgetPostProc.x
 
PATH_FILE = $(EXEC)

SOURCE += functions_utilities_Matteo.cpp Hdf5_particles_calculate_density.cpp functions_on_mass_distribution.cpp functions_make_arrays.cpp functions_Hdf5_writeHDF5.cpp functions_readHDF5.cpp readDensities.cpp Gadget_postProc.cpp

OBJECTS = $(SOURCE:.cpp=.o)

HDF5LIBS = -L$(HDF5_LIB_DIR)
HDF5_INCL = -I$(HDF5_INC_DIR)
INCL += $(HDF5_INCL) 
LIBS += $(HDF5LIBS)

LIBS +=  -lhdf5 -lhdf5_cpp -lz

build: $(EXEC)

$(EXEC): $(OBJECTS)
	$(CC) $(FFLAGS) $(OBJECTS) -o $(EXEC) $(LIBS) $(INCL) 
	$(info $)
	$(info $)
	$(info $)
	$(info $)
	$(info $)

clean:
	$(info $ Cleaning directory...$)	
	$(info $)
	$(info $)
	-rm GadgetPostProc.x *.o
	$(info $)	
	$(info $)
	$(info $)
	$(info $)	
	$(info $)
	$(info $)	
	$(info $)

