include Makefile.inc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

SRCCONTROL =	\
	Control/dyncall.f90      \
	Control/forc_atm.f90     \
	Control/init_arrays.f90  \
	Control/init_lqbc.f90    \
	Control/init_pars.f90    \
	Control/output.f90   

SRCSERVICE =	    \
	Service/input_output_data.f90     \
	Service/oc_algorithms.f90         \
	Service/atm2oc_interpol.f90       \
	Service/basin_parameters.f90      \
	Service/bc_time_interpol.f90      \
	Service/grid_construction.f90     \
	Service/read_write_parameters.f90 \
	Service/rw_ctl_file.f90           \
	Service/time_tools.f90           

SRCFUNCTION =	\
	Function/depth.f90   \
	Function/parallel_ssh_solve.f90	\
	Function/seaice.f90  \
	Function/fluxes.f90  \
	Function/mixing.f90  \
	Function/tem_sal.f90 \
	Function/velocity.f90   	 	

SRCMODULES = 	\
	Modules/mod_constants.f90            \
	Modules/mod_math.f90                 \
	Modules/mod_mpi_parallel_tools.f90   \
	Modules/mod_atm2oc_interpol.f90      \
	Modules/mod_atm_pars.f90             \
	Modules/mod_atm_variables.f90        \
	Modules/mod_basin_grid.f90           \
	Modules/mod_ocean_bc.f90             \
	Modules/mod_ocean_variables.f90      \
	Modules/mod_time_integration.f90 

####################################################

all: inmsom inmsom_clean

#inmsom:
#order is important
#	$(FC) -c $(FCFLAGS) $(SRCMODULES) $(SRCSERVICE) $(SRCFUNCTION) $(SRCCONTROL) inmsom_head.f90      \
#	$(SRCLIB)
#	$(FC) -o inmsom.exe  *.o *.a  $(FCFLAGS)	
#clean:
#	$(RM) *.o *.mod 

inmsom:
	-${FLINKER} -o inmsom.exe $(SRCMODULES) $(SRCSERVICE) $(SRCFUNCTION) $(SRCCONTROL) ${PETSC_KSP_LIB} inmsom_head.f90 

inmsom_clean:
	$(RM) *.o *.mod

run_seq:
	-@${MPIEXEC} -n 1 ./inmsom.exe ${RUN_PETSC_OPTS} ${PETSC_KSP_OPTS} ${PETSC_PC_OPTS}
	
run:
	-@${MPIEXEC} -n ${RUN_PROCS} ./inmsom.exe ${RUN_PETSC_OPTS} ${PETSC_KSP_OPTS} ${PETSC_PC_OPTS}
