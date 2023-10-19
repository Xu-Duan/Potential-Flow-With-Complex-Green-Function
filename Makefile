# By DX from SJTU
TARGET = PAP
COMPILER = gfortran
F90FILES = Constants.f90 Mesh.f90 GreenFunction.f90 BodyCondition.f90 Output.f90 Main.f90
OFILES = $(patsubst %.f90, %.o, $(F90FILES))
Flags = -llapack -lblas

# Main make target is the executable
all: $(TARGET)

# To build the executable we require all object files
$(TARGET): $(OFILES)
# Linking the object files to create the target executable
	$(COMPILER) -o $(TARGET) $(OFILES) $(Flags)

# Recipe how to build object files, which require the *.f90 files
%.o: %.f90
# Compiling one f90 file ($<) to create the corresponding .o file ($@)
	$(COMPILER) -c $< -o $@

clean: 
	rm *.o *.mod

cleanall: 
	rm *.o *.mod PAP