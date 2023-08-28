# Typically this makefile is only called from make in the parent directory,
# in which the variables FC and EXTRA_COMPILE_FLAGS are set.

$(shell mkdir -p Release)

LM_f_SRC = $(shell find Sources/LIBSTELL_minimal -name '*.f' | sed "s|&\./||")
LM_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(LM_f_SRC))))
LM_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(LM_f_SRC))))
LM_f90_SRC = $(shell find Sources/LIBSTELL_minimal -name '*.f90' | sed "s|&\./||")
LM_f90_OBJ = $(addprefix Release/,$(subst .f90,.o,$(notdir $(LM_f90_SRC))))
LM_f90_MOD = $(addprefix Release/,$(subst .f90,.mod,$(notdir $(LM_f90_SRC))))
IC_f_SRC = $(shell find Sources/Initialization_Cleanup -name '*.f' | sed "s|&\./||")
IC_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(IC_f_SRC))))
IC_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(IC_f_SRC))))
IC_f90_SRC = $(shell find Sources/Initialization_Cleanup -name '*.f90' | sed "s|&\./||")
IC_f90_OBJ = $(addprefix Release/,$(subst .f90,.o,$(notdir $(IC_f90_SRC))))
IC_f90_MOD = $(addprefix Release/,$(subst .f90,.mod,$(notdir $(IC_f90_SRC))))
H_f_SRC = $(shell find Sources/Hessian -name '*.f' | sed "s|&\./||")
H_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(H_f_SRC))))
H_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(H_f_SRC))))
G_f_SRC = $(shell find Sources/General -name '*.f' | sed "s|&\./||")
G_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(G_f_SRC))))
G_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(G_f_SRC))))
G_f90_SRC = $(shell find Sources/General -name '*.f90' | sed "s|&\./||")
G_f90_OBJ = $(addprefix Release/,$(subst .f90,.o,$(notdir $(G_f90_SRC))))
G_f90_MOD = $(addprefix Release/,$(subst .f90,.mod,$(notdir $(G_f90_SRC))))
IO_f_SRC = $(shell find Sources/Input_Output -name '*.f' | sed "s|&\./||")
IO_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(IO_f_SRC))))
IO_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(IO_f_SRC))))
NV_f_SRC = $(shell find Sources/NESTOR_vacuum -name '*.f' | sed "s|&\./||")
NV_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(NV_f_SRC))))
NV_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(NV_f_SRC))))
S_f_SRC = $(shell find Sources/Splines -name '*.f' | sed "s|&\./||")
S_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(S_f_SRC))))
S_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(S_f_SRC))))
TS_f_SRC = $(shell find Sources/TimeStep -name '*.f' | sed "s|&\./||")
TS_f_OBJ = $(addprefix Release/,$(subst .f,.o,$(notdir $(TS_f_SRC))))
TS_f_MOD = $(addprefix Release/,$(subst .f,.mod,$(notdir $(TS_f_SRC))))
TS_f90_SRC = $(shell find Sources/TimeStep -name '*.f90' | sed "s|&\./||")
TS_f90_OBJ = $(addprefix Release/,$(subst .f90,.o,$(notdir $(TS_f90_SRC))))
TS_f90_MOD = $(addprefix Release/,$(subst .f90,.mod,$(notdir $(TS_f90_SRC))))

VMEC_TARGET = Release/libvmec.a
VMEC_EXEC = xvmec2000

# The following lines are used in OSX
SCALAPACK_HOME ?= /usr/local/opt/scalapack
FFTW_HOME ?= /usr/local/opt/fftw
HDF5_HOME ?= /usr/local/opt/hdf5
NETCDF_HOME ?= /usr/local/opt/netcdf
LIBS    = -L${SCALAPACK_HOME}/lib -lscalapack
# GCC compiler options
PRECOMP = -cpp -DLINUX -DMPI_OPT -DNETCDF
FC_COMPILE_FLAGS = -I${NETCDF_HOME}/include -I${FFTW_HOME}/include -I Release 
FLINKER = mpif90 -O3 -march=native -ffree-line-length-none


.PHONY: all vmec_clean

all: $(VMEC_TARGET)
exec: $(VMEC_EXEC)

include VMEC2000.dep

%.o: %.f90
Release/%.o: Sources/LIBSTELL_minimal/%.f90
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f90
Release/%.o: Sources/General/%.f90
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f90
Release/%.o: Sources/Initialization_Cleanup/%.f90
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f90
Release/%.o: Sources/TimeStep/%.f90
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/LIBSTELL_minimal/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/General/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/Initialization_Cleanup/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/Hessian/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/Input_Output/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/NESTOR_vacuum/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/Splines/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

%.o: %.f
Release/%.o: Sources/TimeStep/%.f
	$(FLINKER) $(FC_COMPILE_FLAGS) $(PRECOMP) -c -o $@ $<

$(VMEC_TARGET): $(LM_f_OBJ) $(LM_f90_OBJ) $(G_f_OBJ) $(G_f90_OBJ) $(IC_f_OBJ) $(IC_f90_OBJ) $(H_f_OBJ) $(IO_f_OBJ) $(NV_f_OBJ) $(S_f_OBJ) $(TS_f_OBJ) $(TS_f90_OBJ)
	ar rcs $@ $^ 

$(VMEC_EXEC): $(LM_f_OBJ) $(LM_f90_OBJ) $(G_f_OBJ) $(G_f90_OBJ) $(IC_f_OBJ) $(IC_f90_OBJ) $(H_f_OBJ) $(IO_f_OBJ) $(NV_f_OBJ) $(S_f_OBJ) $(TS_f_OBJ) $(TS_f90_OBJ)
	$(FLINKER) $(FC_COMPILE_FLAGS) -o $@ $^ -L${NETCDF_HOME}/lib -lnetcdf -lnetcdff \
	-L${FFTW_HOME}/lib -lfftw3\
	-L${SCALAPACK_HOME}/lib -lscalapack

vmec_clean:
	rm -f Release/*.o Release/*.mod $(VMEC_TARGET)
