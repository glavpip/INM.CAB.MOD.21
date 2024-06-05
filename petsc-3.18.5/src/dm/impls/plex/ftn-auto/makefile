
#requiresdefine   'PETSC_HAVE_FORTRAN'
ALL: lib
CPPFLAGS = ${NETCFD_INCLUDE} ${EXODUSII_INCLUDE}
CFLAGS   =
FFLAGS   =
SOURCEC  = plexdistributef.c plexmetricf.c plexreorderf.c plexcgnsf.c plexfemf.c plexfvmf.c plexindicesf.c plexnaturalf.c plexorientf.c plexpartitionf.c plexpointf.c plexprojectf.c plexrefinef.c plexsubmeshf.c plextreef.c plexexodusiif.c plexegadsf.c plexgeometryf.c plexf.c plexcheckinterfacef.c plexcreatef.c plexgmshf.c plexinterpolatef.c
SOURCEF  =
SOURCEH  = 
DIRS     =
LIBBASE  = libpetscdm
LOCDIR   = src/dm/impls/plex/ftn-auto/
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test
