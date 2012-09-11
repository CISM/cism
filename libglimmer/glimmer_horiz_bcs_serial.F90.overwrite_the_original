!This is a skeleton dummy module for the serial case. It allows serial compilation.

module glimmer_horiz_bcs
  implicit none !If a variable isn't defined, throw a compiler error
  private       !Must declare all public routines in the header.

  !Enforce boundary conditions for a variety of variable types
  public :: horiz_bcs_unstag_scalar   !Unstaggered scalar variables
  public :: horiz_bcs_stag_vector_ew  !Staggered ew-direction vector
  public :: horiz_bcs_stag_vector_ns  !Staggered ns-direction vector
  public :: horiz_bcs_stag_scalar   !Unstaggered scalar variables

  interface horiz_bcs_unstag_scalar
    module procedure horiz_bcs_unstag_scalar_integer_2d
    module procedure horiz_bcs_unstag_scalar_real8_2d
    module procedure horiz_bcs_unstag_scalar_real8_3d
  end interface

  interface horiz_bcs_stag_vector_ew
    module procedure horiz_bcs_stag_vector_ew_real8_3d
  end interface

  interface horiz_bcs_stag_vector_ns
    module procedure horiz_bcs_stag_vector_ns_real8_3d
  end interface

  interface horiz_bcs_stag_scalar
    module procedure horiz_bcs_stag_scalar_integer_2d
    module procedure horiz_bcs_stag_scalar_real8_2d
  end interface


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_2d( a )
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_unstag_scalar_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_integer_2d( a )
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_unstag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_unstag_scalar_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ew_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_stag_vector_ew_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ns_real8_3d( a )
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_stag_vector_ns_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_integer_2d( a )
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_stag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_scalar_real8_2d( a )
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    write(*,*) 'Using a dummy file for serial building. Horizontal BCs are not being enforced!'
    write(*,'(A,I5)') 'Printed from line: ', __LINE__
  end subroutine horiz_bcs_stag_scalar_real8_2d

end module glimmer_horiz_bcs

