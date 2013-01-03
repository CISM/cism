module mpi_mod
! This module wraps the external mpi module

#ifndef NO_MPIMOD
  use mpi
#endif

  implicit none

#ifdef NO_MPIMOD
#include <mpif.h>
#endif

  public

end module mpi_mod
