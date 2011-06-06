module glam_Basal_Proc	

!!! NOTE: This is a module under development that is not supported in the current release of CISM !!!

use glide_types
use glimmer_paramets, only : dp,sp,vel0,tau0,thk0,tim0
use glimmer_physcon,  only : grav, rhow, rhos, scyr
use glimmer_log,      only : write_log

  implicit none;save

!!Variables
  real (kind = dp), dimension(:,:,:), allocatable:: dy  !u,etill,
  real (kind = dp), dimension(:,:), allocatable:: minTauf_init,Hwater_init
  logical, dimension(:,:),   allocatable::tillmask
  integer, parameter :: unin = 90
  
contains
  
  subroutine Basal_Proc_init(ewn,nsn,basalproc,ntem)
    implicit none
    
    !Arguments
    integer,intent(in) :: ewn,nsn
	real (kind = sp),intent (in) :: ntem
    type(glide_basalproc),intent(inout) :: basalproc
    
    !Variables
    real (kind = dp), dimension (ewn-1,nsn-1,basalproc%tnodes) :: por,Ztot,N
    real (kind = dp), dimension (ewn-1,nsn-1) :: stagHwater
    integer :: x,y,i
    character(len=512) :: message

    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
    stop

  end subroutine Basal_Proc_init
  
  
  subroutine Basal_Proc_driver (ewn,      nsn,      upn,  				&
								dt,ubas,vbas,  what, bmlt, basalproc)
								
	use glide_grids, only: stagvarb								
  implicit none
    
    !Arguments
    integer, intent (in) ::ewn, nsn, upn, what
    real (kind = sp), intent(in) :: dt
    real (kind = dp), dimension(:,:), intent (in) :: ubas,vbas
    real (kind = dp), dimension(:,:), intent (in) :: bmlt
    type(glide_basalproc),intent(inout) :: basalproc
    
    !Variables
    real (kind = dp), dimension (ewn-1,nsn-1) :: Ub,stagHwater,stagbmlt
    real (kind = dp) :: f1
    integer :: i
    
    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
    stop

  end subroutine Basal_Proc_driver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  subroutine Basal_Proc_final(basalproc)
	type(glide_basalproc),intent(inout) :: basalproc
  ! Deallocate till variables

    deallocate(dy)
  !  deallocate(basalproc%etill)
  !  deallocate(basalproc%u)
  deallocate (basalproc%till_dz)
    
  end subroutine Basal_Proc_final
 
end module glam_Basal_Proc

  
