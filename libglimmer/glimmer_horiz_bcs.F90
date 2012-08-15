!Stores routines to support the maintainance of the desired horizontal boundary conditions on the global domain.

module glimmer_horiz_bcs
  implicit none !If a variable isn't defined, throw a compiler error
  private       !Must declare all public routines in the header. You're welcome.

  public :: horiz_bcs_unstag_scalar


  interface horiz_bcs_unstag_scalar
     module procedure horiz_bcs_unstag_scalar_real8_3d
  end interface


contains


  subroutine horiz_bcs_unstag_scalar_real8_3d( a )
    use parallel, only: this_rank, north, south, east, west, local_ewn, local_nsn
    implicit none
    real(8),dimension( : , : , : ) :: a
    if (this_rank < west ) then
      write(*,*) this_rank, 'west'
    endif
  end subroutine horiz_bcs_unstag_scalar_real8_3d


end module glimmer_horiz_bcs

