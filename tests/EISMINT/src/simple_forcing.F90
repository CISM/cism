! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_forcing.f90 - part of the GLIMMER ice model       + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module simple_forcing
  !*FD read configuration and generate simple massbalance and 
  !*FD temperature fields

  use glimmer_global, only : sp

  type simple_climate
     ! holds parameters for the simple climate

     integer :: eismint_type = 0
     !*FD select EISMINT experiment
     !*FD \begin{description}
     !*FD \item[{\bf 1}] EISMINT-1 fixed margin
     !*FD \item[{\bf 2}] EISMINT-1 moving margin
     !*FD \end{description}
     real(kind=sp), dimension(2) :: airt = (/ -3.150, 1.e-2 /)  
     !*FD air temperature parameterisation K, K km$^{-3}$
     real(kind=sp), dimension(3) :: nmsb = (/ 0.5, 1.05e-5, 450.0e3 /)
     !*FD mass balance parameterisation m yr$^{-1}$, yr$^{-1}$, m
     real(kind=sp) :: period = 0.
     !*FD EISMINT time-dep climate forcing period, switched off when set to 0
     real(kind=sp) :: mb_amplitude = 0.2
     !*FD EISMINT amplitude of mass balance time-dep climate forcing
  end type simple_climate

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_SIMPLE_FORCING
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_SIMPLE_FORCING
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_SIMPLE_FORCING
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_SIMPLE_FORCING
!MH!#endif

  subroutine simple_initialise(climate,config)
    !*FD initialise simple climate model
    use paramets, only: thk0, acc0, scyr
    use glimmer_config
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

  
    call simple_readconfig(climate,config)
    call simple_printconfig(climate)

    ! scale parameters
    select case(climate%eismint_type)
    case(1)
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
    case(2)
       climate%airt(2) = climate%airt(2) * thk0
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
       climate%nmsb(2) = climate%nmsb(2) / (acc0 * scyr)
    case(3)
       climate%nmsb(1) = climate%nmsb(1) / (acc0 * scyr)
       climate%nmsb(2) = climate%nmsb(2) / (acc0 * scyr)
    end select
       
  end subroutine simple_initialise

  subroutine simple_readconfig(climate, config)
    !*FD read configuration
    use glimmer_log
    use glimmer_config
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    real(kind=sp), dimension(:), pointer :: dummy

    call GetSection(config,section,'EISMINT-1 fixed margin')
    if (associated(section)) then
       climate%eismint_type = 1
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       climate%airt = (/-34.15, 8.e-8/)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,1)
       climate%nmsb = (/0.3, 0.0, 0.0/)
       if (associated(dummy)) then
          climate%nmsb(1) = dummy(1)
       end if
       call GetValue(section,'period',climate%period)
       call GetValue(section,'mb_amplitude',climate%mb_amplitude)
       return       
    end if
    call GetSection(config,section,'EISMINT-1 moving margin')
    if (associated(section)) then
       climate%eismint_type = 2
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       call GetValue(section,'period',climate%period)
       climate%mb_amplitude = 100000.
       call GetValue(section,'mb_amplitude',climate%mb_amplitude)
       return
    end if
    call GetSection(config,section,'EISMINT-2')
    if (associated(section)) then
       climate%eismint_type = 3
       dummy=>NULL()
       call GetValue(section,'temperature',dummy,2)
       if (associated(dummy)) then
          climate%airt = dummy
          deallocate(dummy)
          dummy=>NULL()
       else
          climate%airt = (/-35., 1.67e-5/)
       end if
       call GetValue(section,'massbalance',dummy,3)
       if (associated(dummy)) then
          climate%nmsb = dummy
          deallocate(dummy)
          dummy=>NULL()
       end if
       return
    end if    

    call write_log('No EISMINT forcing selected',GM_FATAL)
  end subroutine simple_readconfig

  subroutine simple_printconfig(climate)
    !*FD print simple climate configuration
    use glimmer_log
    implicit none
    type(simple_climate) :: climate   !*FD structure holding climate info
    character(len=100) :: message

    call write_log_div
    select case(climate%eismint_type)
    case(1)
       call write_log('EISMINT-1 fixed margin configuration')
       call write_log('------------------------------------')
       write(message,*) 'temperature  : ',climate%airt(1)
       call write_log(message)
       write(message,*) '               ',climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',climate%nmsb(1)
       call write_log(message)
       write(message,*) 'period       : ',climate%period
       call write_log(message)
       if (climate%period .gt.0) then
          write(message,*) 'mb amplitude : ',climate%mb_amplitude
          call write_log(message)
       end if
    case(2)
       call write_log('EISMINT-1 moving margin configuration')
       call write_log('-------------------------------------')
       write(message,*) 'temperature  : ',climate%airt(1)
       call write_log(message)
       write(message,*) '               ',climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',climate%nmsb(1)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(2)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(3)
       call write_log(message)
       write(message,*) 'period       : ',climate%period
       call write_log(message)
       if (climate%period .gt.0) then
          write(message,*) 'mb amplitude : ',climate%mb_amplitude
          call write_log(message)
       end if
    case(3)
       call write_log('EISMINT-2')
       call write_log('---------')
       write(message,*) 'temperature  : ',climate%airt(1)
       call write_log(message)
       write(message,*) '               ',climate%airt(2)
       call write_log(message)
       write(message,*) 'massbalance  : ',climate%nmsb(1)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(2)
       call write_log(message)
       write(message,*) '               ',climate%nmsb(3)
       call write_log(message)
    end select
    call write_log('')
  end subroutine simple_printconfig

  subroutine simple_massbalance(climate,model,time)
    !*FD calculate simple mass balance
    use glimmer_global, only:rk
    use glide_types
    use paramets, only : len0, acc0, scyr
    use physcon, only : pi
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance
    real(kind=rk), intent(in) :: time                !*FD current time

    ! local variables
    integer  :: ns,ew
    real :: dist, ewct, nsct, grid, rel
    real :: periodic_bc = 1.

    ewct = real(model%general%ewn+1) / 2.0
    nsct = real(model%general%nsn+1) / 2.0
    grid = model%numerics%dew * len0

    periodic_bc = real(1-model%options%periodic_ew)

    select case(climate%eismint_type)
    case(1)
       ! EISMINT-1 fixed margin
       model%climate%acab(:,:) = climate%nmsb(1)
       if (climate%period.ne.0) then
          model%climate%acab(:,:) = model%climate%acab(:,:) + climate%mb_amplitude * sin(2.*pi*time/climate%period)/ (acc0 * scyr)
       end if
    case(2)
       ! EISMINT-1 moving margin       
       if (climate%period.ne.0) then
          rel = climate%nmsb(3) + climate%mb_amplitude*sin(2.*pi*time/climate%period)
       else
          rel = climate%nmsb(3)
       end if

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic_bc*(real(ew) - ewct)**2 + (real(ns) - nsct)**2)
             model%climate%acab(ew,ns) = min(climate%nmsb(1), climate%nmsb(2) * (rel - dist))
          end do
       end do
    case(3)
       ! EISMINT-2
       rel = climate%nmsb(3)

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic_bc*(real(ew) - ewct)**2 + (real(ns) - nsct)**2)
             model%climate%acab(ew,ns) = min(climate%nmsb(1), climate%nmsb(2) * (rel - dist))
          end do
       end do       
    end select
  end subroutine simple_massbalance

  subroutine simple_surftemp(climate,model,time)
    !*FD calculate simple air surface temperature
    use glide_types
    use glimmer_global, only:rk
    use paramets, only : len0
    use physcon, only : pi
    implicit none
    type(simple_climate) :: climate         !*FD structure holding climate info
    type(glide_global_type) :: model        !*FD model instance
    real(kind=rk), intent(in) :: time                !*FD current time

    ! local variables
    integer  :: ns,ew
    real :: dist, ewct, nsct, grid
    real :: periodic_bc = 1.

    ewct = real(model%general%ewn+1) / 2.0
    nsct = real(model%general%nsn+1) / 2.0
    grid = model%numerics%dew * len0

    periodic_bc = real(1-model%options%periodic_ew)

    select case(climate%eismint_type)
    case(1)
       ! EISMINT-1 fixed margin
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * max(periodic_bc*abs(real(ew) - ewct),abs(real(ns) - nsct))*1e-3
             model%climate%artm(ew,ns) = climate%airt(1) + climate%airt(2) * dist*dist*dist
          end do
       end do
       if (climate%period.ne.0) then
          model%climate%artm(:,:) = model%climate%artm(:,:) + 10.*sin(2.*pi*time/climate%period)
       end if
    case(2)
       ! EISMINT-1 moving margin
       model%climate%artm(:,:) = climate%airt(1) - model%geometry%thck(:,:) * climate%airt(2)
       if (climate%period.ne.0) then
          model%climate%artm(:,:) = model%climate%artm(:,:) + 10.*sin(2.*pi*time/climate%period)
       end if
    case(3)
       ! EISMINT-2
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             dist = grid * sqrt(periodic_bc*(real(ew) - ewct)**2 + (real(ns) - nsct)**2)
             model%climate%artm(ew,ns) = climate%airt(1)+climate%airt(2) * dist
          end do
       end do       
    end select
  end subroutine simple_surftemp
end module simple_forcing
