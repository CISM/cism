!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glimmer_restart_gcm

!BOP
! !MODULE: glimmer_restart_gcm

! !DESCRIPTION:
!  Contains routines for specialized glimmer restarts called by gcm's
!
! !REVISION HISTORY:
!
! !USES:

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glimmer_read_restart_gcm

!----------------------------------------------------------------------
!
!   module variables
!
!----------------------------------------------------------------------

!EOP
!BOC
!EOC
!***********************************************************************
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glimmer_read_restart_gcm
! !INTERFACE:

   subroutine glimmer_read_restart_gcm(model, restart_filename)

    use glide_types
    implicit none
    type(glide_global_type), intent(inout) :: model
    character(*),            intent(in   ) :: restart_filename

    ! local variables
    type(glimmer_nc_input),  pointer :: ic => null()

!-----------------------------------------------------------------------

    ! create the input unit
    allocate(ic)
    ic%get_time_slice = 1
    ic%nc%filename    = trim(restart_filename)
    ic%nc%vars        = ' hot '
    ic%nc%hotstart    = .true.
    ic%nc%vars_copy   = ic%nc%vars

    ! add the input unit to the model
    ! note that the model will do the actual reading of data
    model%funits%in_first => ic

  end subroutine glimmer_read_restart_gcm

end module glimmer_restart_gcm
