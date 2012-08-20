!Stores routines to support the maintainance of the desired horizontal boundary conditions on the global domain.

!Currently, these routines only support closed BC's with a free-slip transverse velocity. Once these are tested,
!more will be added with namelist variables.

!Each routine has been tested with a field to make sure that at least, they aren't indexing out of bounds. Testing
!for correctness is to come soon.

module glimmer_horiz_bcs
  implicit none !If a variable isn't defined, throw a compiler error
  private       !Must declare all public routines in the header.

  !Enforce boundary conditions for a variety of variable types
  public :: horiz_bcs_unstag_scalar   !Unstaggered scalar variables
  public :: horiz_bcs_stag_vector_ew  !Staggered ew-direction vector
  public :: horiz_bcs_stag_vector_ns  !Staggered ns-direction vector

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


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_2d( a )
    use parallel, only: this_rank, north, south, east, west, own_ewn, own_nsn, lhalo, uhalo
    implicit none
    real(8),dimension(:,:), intent(inout) :: a
    integer :: i

    !If this process is on a domain boundary, loop over the halo (ghost cells) and populate them

    if ( this_rank > north ) then   !I am on the north boundary
      do i = 1 , uhalo
        a(:,lhalo+own_nsn+i) = a(:,lhalo+own_nsn+1-i)
      enddo
    endif
    if ( this_rank > east ) then    !I am on the east boundary
      do i = 1 , uhalo
        a(lhalo+own_ewn+i,:) = a(lhalo+own_ewn+1-i,:)
      enddo
    endif
    if ( this_rank < south ) then   !I am on the south boundary
      do i = 1 , lhalo
        a(:,lhalo+1-i) = a(:,lhalo+i)
      enddo
    endif
    if ( this_rank < west ) then    !I am on the west boundary
      do i = 1 , lhalo
        a(lhalo+1-i,:) = a(lhalo+i,:)
      enddo
    endif
  end subroutine horiz_bcs_unstag_scalar_real8_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_integer_2d( a )
    use parallel, only: this_rank, north, south, east, west, own_ewn, own_nsn, lhalo, uhalo
    implicit none
    integer,dimension(:,:), intent(inout) :: a
    integer :: i
    if ( this_rank > north ) then
      do i = 1 , uhalo
        a(:,lhalo+own_nsn+i) = a(:,lhalo+own_nsn+1-i)
      enddo
    endif
    if ( this_rank > east ) then
      do i = 1 , uhalo
        a(lhalo+own_ewn+i,:) = a(lhalo+own_ewn+1-i,:)
      enddo
    endif
    if ( this_rank < south ) then
      do i = 1 , lhalo
        a(:,lhalo+1-i) = a(:,lhalo+i)
      enddo
    endif
    if ( this_rank < west ) then
      do i = 1 , lhalo
        a(lhalo+1-i,:) = a(lhalo+i,:)
      enddo
    endif
  end subroutine horiz_bcs_unstag_scalar_integer_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_unstag_scalar_real8_3d( a )
    use parallel, only: this_rank, north, south, east, west, own_ewn, own_nsn, lhalo, uhalo
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i
    if ( this_rank > north ) then
      do i = 1 , uhalo
        a(:,:,lhalo+own_nsn+i) = a(:,:,lhalo+own_nsn+1-i)
      enddo
    endif
    if ( this_rank > east ) then
      do i = 1 , uhalo
        a(:,lhalo+own_ewn+i,:) = a(:,lhalo+own_ewn+1-i,:)
      enddo
    endif
    if ( this_rank < south ) then
      do i = 1 , lhalo
        a(:,:,lhalo+1-i) = a(:,:,lhalo+i)
      enddo
    endif
    if ( this_rank < west ) then
      do i = 1 , lhalo
        a(:,lhalo+1-i,:) = a(:,lhalo+i,:)
      enddo
    endif
  end subroutine horiz_bcs_unstag_scalar_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ew_real8_3d( a )
    use parallel, only: this_rank, north, south, east, west, own_ewn, own_nsn, &
                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( this_rank > north ) then
      do i = 1 , nhalo
        a(:,:,shalo+ns_npts+i) = a(:,:,shalo+ns_npts-i)
      enddo
    endif
    if ( this_rank > east ) then
      do i = 1 , ehalo
        a(:,whalo+ew_npts+i,:) = -a(:,whalo+ew_npts-i,:)
      enddo
    endif
    if ( this_rank < south ) then
      do i = 1 , shalo
        a(:,:,shalo+1-i) = a(:,:,shalo+1+i)
      enddo
    endif
    if ( this_rank < west ) then
      do i = 1 , whalo
        a(:,whalo+1-i,:) = -a(:,whalo+1+i,:)
      enddo
    endif
    !Normal velocities are zero at boundaries
    a(:,whalo+1      ,:) = 0.D0 
    a(:,whalo+ew_npts,:) = 0.D0
    !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
    !southern boundary, however, and that must be interpolated.
    !For this, some assumptions must be made: shalo >= 1 & nsn >= 2. Using three data points allows the
    !inclusion of curvature in this interpolation using shalo, shalo+2, and shalo+3
    a(:,:,shalo+1) = a(:,:,shalo) / 3.D0 + a(:,:,shalo+2) - a(:,:,shalo+3) / 3.D0

  end subroutine horiz_bcs_stag_vector_ew_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine horiz_bcs_stag_vector_ns_real8_3d( a )
    use parallel, only: this_rank, north, south, east, west, own_ewn, own_nsn, &
                        staggered_nhalo, staggered_ehalo, staggered_shalo, staggered_whalo
    implicit none
    real(8),dimension(:,:,:), intent(inout) :: a
    integer :: i
    integer :: nhalo, ehalo, shalo, whalo !Sizes of halos
    integer :: ew_npts, ns_npts           !Number of staggered points in ew and ns directions
    nhalo = staggered_nhalo
    ehalo = staggered_ehalo
    shalo = staggered_shalo-1 !Technically, domains on the south or west boundaries do not "own" the south or 
    whalo = staggered_whalo-1 !west edges. So we must specify these values. So, these halos are reduced by 1.
    ew_npts = own_ewn+1 !number of physical domain points this process is responsible for in ew direction.
    ns_npts = own_nsn+1 !number of physical domain points this process is responsible for in ns direction.

    !Physical domain is mirrored at all boundaries and negated at normal boundaries
    if ( this_rank > north ) then
      do i = 1 , nhalo
        a(:,:,shalo+ns_npts+i) = -a(:,:,shalo+ns_npts-i)
      enddo
    endif
    if ( this_rank > east ) then
      do i = 1 , ehalo
        a(:,whalo+ew_npts+i,:) = a(:,whalo+ew_npts-i,:)
      enddo
    endif
    if ( this_rank < south ) then
      do i = 1 , shalo
        a(:,:,shalo+1-i) = -a(:,:,shalo+1+i)
      enddo
    endif
    if ( this_rank < west ) then
      do i = 1 , whalo
        a(:,whalo+1-i,:) = a(:,whalo+1+i,:)
      enddo
    endif
    !Normal velocities are zero at boundaries
    a(:,:,shalo+1      ) = 0.D0 
    a(:,:,shalo+ns_npts) = 0.D0
    !For slip transverse BC's, the northern boundary is left alone. The velocity solve does not treat the
    !western boundary, however, and that must be interpolated.
    !For this, some assumptions must be made: whalo >= 1 & ewn >= 2. Using three data points allows the
    !inclusion of curvature in this interpolation using whalo, whalo+2, and whalo+3
    a(:,whalo+1,:) = a(:,whalo,:) / 3.D0 + a(:,whalo+2,:) - a(:,whalo+3,:) / 3.D0

  end subroutine horiz_bcs_stag_vector_ns_real8_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module glimmer_horiz_bcs

