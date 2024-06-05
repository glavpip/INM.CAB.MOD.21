MAKE = make

# For gfortran:
#FFLAGS = -cpp -w -ffree-line-length-0 -O3 -Wtabs -fopenmp ${PETSC_FC_INCLUDES} -I./Inc
#CFLAGS = ${PETSC_CC_INCLUDES}

# For ifort:
FFLAGS = -cpp -w -ffree-line-length-0 -O3 -Wtabs -fopenmp -mcmodel medium ${PETSC_FC_INCLUDES} -I./Inc
CFLAGS = ${PETSC_CC_INCLUDES}

RM = rm -rf

# PETSc common options
# example: -draw_pause -1
RUN_PETSC_OPTS = -draw_pause -1

# PETSc solver configuration
# example: -ksp_view -ksp_plot_eigenvalues_explicitly
PETSC_KSP_OPTS = -ksp_type gmres -ksp_atol 1.0e-16 -ksp_rtol 1.0e-08

# PETSc precond configuration
# example:
PETSC_PC_OPTS = -pc_type asm -pc_asm_overlap 1 -sub_pc_type ilu -sub_pc_factor_levels 8

# Number of parallel procs for "make run"
RUN_PROCS = 4
