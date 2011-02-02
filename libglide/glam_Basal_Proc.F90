module glam_Basal_Proc	
	
	
use glide_types
use glimmer_paramets, only : dp,sp,vel0,tau0_glam,thk0,tim0
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

!allocate basal processes variables
	allocate (dy(ewn-1,nsn-1,basalproc%tnodes-1));		dy=0.6d0
	allocate (minTauf_init(ewn-1,nsn-1));				minTauf_init=5000
	allocate (Hwater_init(ewn-1,nsn-1)); 				Hwater_init=3
	allocate (tillmask(ewn-1,nsn-1)); 					tillmask=.true.
	allocate (basalproc%till_dz(basalproc%tnodes)); basalproc%till_dz=1.0d0



    if (basalproc%till_hot.eq.1) then
    
       !From restart file, the following variables are known: basalproc%u, basalproc%minTauf and basalproc%etill
       basalproc%minTauf=basalproc%minTauf*tau0_glam
       por=basalproc%etill/(1+basalproc%etill)
       dy=basalproc%Zs*(1+basalproc%etill(:,:,1:basalproc%tnodes-1))/(basalproc%tnodes-1)
       stagHwater=0.0
       do i=2,basalproc%tnodes
          stagHwater(:,:)=stagHwater+dy(:,:,i-1)*(por(:,:,i-1)+por(:,:,i))/2
       enddo

       
       
    else if (basalproc%till_hot.eq.0) then

       basalproc%minTauf=basalproc%minTauf*tau0_glam
       N(:,:,1)=basalproc%minTauf/basalproc%fric
       do i=2,basalproc%tnodes
          N(:,:,i)=N(:,:,1)
       end do
       basalproc%etill=basalproc%etillo-basalproc%Comp*log10(N/basalproc%No)
       por=basalproc%etill/(1+basalproc%etill)
       dy=basalproc%Zs*(1+basalproc%etill(:,:,1:basalproc%tnodes-1))/(basalproc%tnodes-1) 
       Ztot(:,:,1)=0.0
       stagHwater=0.0
       do i=2,basalproc%tnodes
          Ztot(:,:,i)=Ztot(:,:,i-1)+dy(:,:,i-1)
          stagHwater(:,:)=stagHwater+dy(:,:,i-1)*(por(:,:,i-1)+por(:,:,i))/2
       enddo
       basalproc%u=(rhos-rhow)*(1-por)*grav*Ztot-N

      
    end if
    
    	!Define tillmask based on values of Tauf in the initial file
		!Wherever Tauf has been set to 10Pa, tillmask=.false.,
		!e.g., we will not calculate till properties
		where (basalproc%minTauf.eq.10.0) tillmask=.false.
		minTauf_init=basalproc%minTauf
		Hwater_init=stagHwater
    	!print*,'mean Tauf init=',sum(minTauf_init)/((ewn-1)*(nsn-1))
    	!Calculate Hwater on normal grid - using zero gradient as BC
		call stag2norm(ewn,nsn,stagHwater,basalproc%Hwater)	    
        
        write(message,*) 'Till layer has been initialized'
       	call write_log(message) 
       	print*,'minTauf=',sum(basalproc%minTauf)/((ewn-1)*(nsn-1))
!       	print*,'dy=',sum(dy)/((ewn-1)*(nsn-1)*(basalproc%tnodes-1))
!       	print*,'u=',sum(basalproc%u)/((ewn-1)*(nsn-1)*basalproc%tnodes)
!       	print*,'etill=',sum(basalproc%etill)/((ewn-1)*(nsn-1)*basalproc%tnodes)
!       	print*,'Hwater=',sum(stagHwater)/((ewn-1)*(nsn-1))
       	
       	
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

    
    !Calculate basal melt rate on staggered grid
  	call stagvarb(bmlt, stagbmlt,ewn,nsn)
	!Re-Scale bmlt (scale2d_f1=scyr*thk0/tim0):
	stagbmlt=stagbmlt*(scyr*thk0/tim0)
    
    
!    !Calculate the magnitude of basal velocity, in m/yr
	Ub=scyr*vel0* (sqrt(ubas(:,:)**2+vbas(:,:)**2))

!	print*,'begin basal_proc_driver'
!	print*,'mean Tauf=',sum(basalproc%minTauf)/((ewn-1)*(nsn-1))
!	print*,'mean bmlt=',sum(stagbmlt)/((ewn-1)*(nsn-1))
!	print*,'mean Ub=',sum(Ub)/((ewn-1)*(nsn-1))
	
    select case(what)

    case(1)
    call Till_FullRes  (ewn,nsn,dt,stagbmlt,Ub,basalproc%tnodes,			&
    					basalproc%Kh,basalproc%Cv, basalproc%etillo,        &
    					basalproc%fric, basalproc%Comp, basalproc%No,  		&
    					basalproc%Zs, stagHwater,basalproc%minTauf,			&
    					basalproc%u,basalproc%etill)

    case(2)
    call Till_FastCalc (dt,stagbmlt,basalproc%aconst,basalproc%bconst,      &
    					basalproc%Zs,basalproc%minTauf,stagHwater,basalproc%etill) 
    
    end select

	!Calculate Hwater on normal grid - using zero gradient as BC
	call stag2norm(ewn,nsn,stagHwater,basalproc%Hwater)
	print*,'ENDING WITH mean Tauf=',sum(basalproc%minTauf)/((ewn-1)*(nsn-1))
	! Scale minTauf
	basalproc%minTauf=basalproc%minTauf/tau0_glam


 
    
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
  
  
  
  subroutine stag2norm(ewn,nsn,stagvar,normvar)
  use glimmer_paramets, only : dp	
  implicit none
  
  integer, intent(in) :: ewn,nsn
  real (kind = dp), intent(in), dimension(:,:) :: stagvar
  real (kind = dp), intent(out), dimension (:,:) :: normvar

    normvar(1:ewn,1:nsn)=0.0
    normvar(2:ewn-1,2:nsn-1)=(stagvar(1:ewn-2,2:nsn-1)+stagvar(2:ewn-1,2:nsn-1)+ & 
    						stagvar(1:ewn-2,1:nsn-2)+stagvar(2:ewn-1,1:nsn-2)) / 4.0d0

    !Apply zero-gradient to field on normal grid - using swapbnmelt
    call swapbnmelt(0,size(normvar(1,:)),normvar(1,:),normvar(2,:),normvar(ewn,:),normvar(ewn-1,:))
    call swapbnmelt(0,size(normvar(:,1)),normvar(:,1),normvar(:,2),normvar(:,nsn),normvar(:,nsn-1))

  end subroutine stag2norm	

  subroutine swapbnmelt(bc,lgth,a,b,c,d)
	use glimmer_paramets, only : dp        
    implicit none
    
    integer, intent(in) :: bc,lgth
    real (kind = dp), intent(in), dimension(lgth) :: b, d
    real (kind = dp), intent(out), dimension (lgth):: a, c
    
    if (bc==0) then
       a = b
       c = d
    end if
    return
    
  end subroutine swapbnmelt


subroutine Till_FullRes(ewn,nsn,dt,bdot,Ub,						&
						tnodes,Kh,Cv,etillo,fric,Comp,No,Zs,		&
						Hwater,minTauf,u,etill)
 !Largely follows Bougamont et al. 2003 (JGR), with explicit till layer properties calculations 
 !Includes vertical mixing as in Christoffersen et al. 2003
	
	use glimmer_paramets, only : dp,sp
	use glimmer_physcon,  only : grav, rhow, rhos, scyr

  implicit none
  
  integer, intent(in) :: ewn,nsn,tnodes
  real (kind=sp), intent(in):: dt
  real (kind=dp),intent(in) ::Kh,Cv,etillo,fric,Comp,No,Zs
  real (kind = dp), dimension(:,:),intent(in) :: bdot ! basal melt rate m.yr
  real (kind = dp), dimension(:,:),intent(in) :: Ub   ! velocity m/yr
  real (kind = dp), dimension(:,:,:),intent(inout) :: u  !Excess pore pressure (Pa)
  real (kind = dp), dimension(:,:,:),intent(inout) :: etill !Till void ratio (ND)
  real (kind = dp), dimension(:,:),intent(out) :: minTauf   ! 
  real (kind = dp), dimension(:,:),intent(out) :: Hwater   ! 

  !Local variables
  real, parameter :: f=1e-3
  real (kind = dp), dimension(ewn-1,nsn-1,tnodes) :: uold,Tauf,deltaU,por
  real (kind = dp), dimension(ewn-1,nsn-1,tnodes) :: vw   !vertical water flow, m/s
  real (kind = dp), dimension(ewn-1,nsn-1) :: du  
  integer :: i

  !Boundary conditions at the ice/till interface
  

  uold=u 
  du=dy(:,:,1)*(bdot/scyr)*rhow*grav/Kh  ! Removed *dt so that du has indeed unit Pa
  u(:,:,1)=uold(:,:,1) + (Cv*scyr*dt/dy(:,:,1)**2)*(uold(:,:,2)-uold(:,:,1)+du)
  
  i=2
  do while (i.lt.tnodes)
     u(:,:,i)=uold(:,:,i) + (Cv*scyr*dt/(dy(:,:,i)**2))*(uold(:,:,i+1)-2*uold(:,:,i)+uold(:,:,i-1)) 
     i=i+1
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set up two possible lower boundary conditions:

  ! 1. ZERO-FLUX assuming an imaginary nodes at (tnodes+1) 
  u(:,:,tnodes)=uold(:,:,tnodes)+(Cv*scyr*dt/(dy(:,:,tnodes-1)**2))*(uold(:,:,tnodes-1)-uold(:,:,tnodes))

  ! 2. CONSTANT FLUX assuming an imaginary nodes at (tnodes+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!Calculate vertical velocity from darcy equation:
  vw(:,:,1)=bdot(:,:)/scyr       !Upper BC, vertical vel=bdot, in m/sec
  vw(:,:,2:tnodes)=Kh/(rhow*grav)*(1/dy(:,:,1:tnodes-1))*(u(:,:,1:tnodes-1)-u(:,:,2:tnodes))
  
  
!From there, new void ratio distribution, flux in - flux out

	etill(:,:,1:tnodes-1)=etill(:,:,1:tnodes-1)+((vw(:,:,1:tnodes-1)-vw(:,:,2:tnodes))*scyr)*dt/(Zs/(tnodes-1))
	etill(:,:,tnodes)=etill(:,:,tnodes-1)
  where (etill.lt.0.15)
     etill=0.15
  end where
  
  Tauf=No*fric*10**(-(etill-etillo)/Comp)
  minTauf=minval(Tauf(:,:,1:tnodes-1),3)  
  dy(:,:,1:tnodes-1)=(Zs/(tnodes-1))*(1+(etill(:,:,1:tnodes-1)+etill(:,:,2:tnodes))/2)
  por=etill/(1+etill)

  Hwater=0.
  do i=2,tnodes
     Hwater=Hwater+dy(:,:,i-1)*((por(:,:,i-1)+por(:,:,i))/2)
  enddo

	!Reset minTauf values where tillmask=.false.
	where (.not.tillmask)
		minTauf=minTauf_init
		Hwater=Hwater_init
	end where

end subroutine Till_FullRes

subroutine Till_FastCalc(dt,bdot,aconst,bconst,Zs,minTauf,Hwater,etill) 
 !Largely follows Bougamont et al. 2003 (JGR)
 !>>>>>>>Till layer is represented with only one node!!  

  use glimmer_paramets, only : dp,sp
  
  implicit none
  
  !Arguments
  real (kind=sp), intent(in):: dt
  real (kind=dp), intent(in):: aconst,bconst,Zs
  real (kind = dp), dimension(:,:),intent(in) :: bdot ! basal melt rate m.yr
  real (kind = dp), dimension(:,:,:),intent(inout) :: etill  !Till void ratio
  real (kind = dp), dimension(:,:),intent(out) :: minTauf ! basal melt rate m.yr
  real (kind = dp), dimension(:,:),intent(out) :: Hwater ! basal melt rate m.yr
 
  etill(:,:,1)=etill(:,:,1)+dt*bdot/Zs
  where (etill.lt.0.15)
     etill=0.15
  end where

  minTauf=aconst*exp(-bconst*etill(:,:,1))
  Hwater=Zs*(etill(:,:,1)/(1+etill(:,:,1)))

	!Reset minTauf values where tillmask=.false.
	where (.not.tillmask)
		minTauf=minTauf_init
		Hwater=Hwater_init
	end where


end subroutine Till_FastCalc

 
end module glam_Basal_Proc

  
