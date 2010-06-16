! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_types.f90 - part of the Glimmer-CISM ice model     + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_types

  !*FD Holds type definitions for the derived types used by each 
  !*FD instance of the ice model. Originally, each of these types
  !*FD was a module containing variables, which were used as containers
  !*FD for global variables. However, the need to allow for multiple
  !*FD ice model instances meant that the nested derived types were instituted
  !*FD instead. However, there is probably one too many levels in this scheme. 
  !*FD It would be better if the different types here were contained in the 
  !*FD higher-level instance type (\texttt{glimmer\_params}), rather than 
  !*FD the intermediate model type (\texttt{glimmer\_global\_type}). 
  !*FD 
  !*FD Note that this \emph{is} now where the defaults are defined for these
  !*FD variables.
 
  use glimmer_sparse
  use glimmer_global
  use glimmer_ncdf
  use isostasy_types
  use profile
  use glimmer_coordinates
  use glimmer_map_types, pi_dummy=>pi

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_general

    !*FD Holds fundamental parameters of the ice model geometry.

    integer :: ewn = 0  !*FD The number of grid-points in the E-W direction.
    integer :: nsn = 0  !*FD The number of grid-points in the N-S direction.
    integer :: upn = 1  !*FD The number of vertical levels in the model.

    type(coordsystem_type) :: ice_grid  !*FD coordinate system of the ice grid
    type(coordsystem_type) :: velo_grid !*FD coordinate system of the velocity grid

    real(sp), dimension(:),pointer :: x0 => null() !original x0 grid 
    real(sp), dimension(:),pointer :: y0 => null() !original y0 grid
    real(sp), dimension(:),pointer :: x1 => null() !original x1 grid
    real(sp), dimension(:),pointer :: y1 => null() !original y1 grid
  end type glide_general

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_options

    !*FD Holds user options controlling the methods used in the ice-model
    !*FD integration.

    integer :: whichtemp = 1

    !*FD Method of ice temperature calculation:
    !*FD \begin{description} 
    !*FD \item[0] Set column to surface air temperature
    !*FD \item[1] Do full temperature solution (also find vertical velocity
    !*FD and apparent vertical velocity)
    !*FD \end{description}

    integer :: whichflwa = 0

    !*FD Method for calculating flow factor $A$:
    !*FD \begin{description} 
    !*FD \item[0] \emph{Patterson and Budd} relationship 
    !*FD \item[1] \emph{Patterson and Budd} relationship, 
    !*FD with temperature set to $-10^{\circ}\mathrm{C}$ 
    !*FD \item[2] Set equal to $1\times 10^{-16}\,\mathrm{yr}^{-1}
    !*FD \,\mathrm{Pa}^{-n}$
    !*FD \end{description}

    integer :: whichbwat = 2

    !*FD Basal water depth: 
    !*FD \begin{description} 
    !*FD \item[0] Calculated from local basal water balance 
    !*FD \item[1] as {\bf 0}, including constant horizontal flow 
    !*FD \item[2] Set to zero everywhere 
    !*FD \end{description}

    integer :: whichmarn = 1

    !*FD Ice thickness: 
    !*FD \begin{description} 
    !*FD \item[0] No action 
    !*FD \item[1] Set thickness to zero if floating 
    !*FD \item[2] Set thickness to zero if relaxed bedrock is more 
    !*FD than certain water depth  
    !*FD \item[3] Lose fraction of ice when edge cell
    !*FD \end{description}

    integer :: whichbtrc = 0

    !*FD Basal slip coefficient: 
    !*FD \begin{description}
    !*FD \item[0] Set equal to zero everywhere
    !*FD \item[1] Set (non--zero) constant
    !*FD \item[2] Set to (non--zero) constant where where temperature is at pressure melting point of ice, otherwise to zero
    !*FD \item[3] \texttt{tanh} function of basal water depth 
    !*FD \end{description}

    integer :: whichevol = 0

    !*FD Thickness evolution method:
    !*FD \begin{description}
    !*FD \item[0] Pseudo-diffusion approach 
    !*FD \item[2] Diffusion approach (also calculates velocities) 
    !*FD \end{description}

    integer :: whichwvel = 0

    !*FD Vertical velocities: 
    !*FD \begin{description}
    !*FD \item[0] Usual vertical integration 
    !*FD \item[1] Vertical integration constrained so that 
    !*FD upper kinematic B.C. obeyed 
    !*FD \end{description}

    integer :: whichrelaxed = 0
    !*FD relaxed topography:
    !*FD \begin{description}
    !*FD \item[0] get relaxed topo from separate variable
    !*FD \item[1] first time slice of input topo is relaxed
    !*FD \item[2] first time slice of input topo is in isostatic equilibrium
    !*FD \end{description}

    integer :: hotstart = 0
    !*FD hotstart the model
    !*FD \begin{description}
    !*FD \item[0] normal start-up
    !*FD \item[1] hotstart model from previous run
    !*FD \end{description}

    integer :: periodic_ew = 0
    !*FD \begin{description}
    !*FD \item[0] no periodic EW boundary conditions
    !*FD \item[1] periodic EW boundary conditions
    !*FD \end{description}

    integer :: gthf = 0
    !*FD \begin{description}
    !*FD \item[0] no geothermal heat flux calculations
    !*FD \item[1] calculate gthf using 3d diffusion
    !*FD \end{description}

    integer :: which_sigma = 0
    !*FD \begin{description}
    !*FD \item[0] calculate sigma coordinates
    !*FD \item[1] sigma coordinates are given in external file
    !*FD \item[2] sigma coordinates are given in configuration file
    !*FD \end{description}

    integer :: basal_mbal = 0
    !*FD \begin{description}
    !*FD \item[0] Basal melt rate not included in continuity equation
    !*FD \item[1] Basal melt rate included in continuity equation
    !*FD \end{description}

  end type glide_options

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geometry

    !*FD Holds fields and other information relating to the
    !*FD geometry of the ice sheet and bedrock.

    real, dimension(:,:), pointer :: temporary0 => null()
    !*FD temporary array used for masking velocity grid
    real, dimension(:,:), pointer :: temporary1 => null()
    !*FD temporary array used for masking temperature grid

    real(dp),dimension(:,:),pointer :: thck => null()
    !*FD The thickness of the ice, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: usrf => null()
    !*FD The elevation of the upper ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: lsrf => null() 
    !*FD The elevation of the lower ice surface, divided by \texttt{thk0}.

    real(dp),dimension(:,:),pointer :: topg => null() 
    !*FD The elevation of the topography, divided by \texttt{thk0}.

    real(dp),dimension(:,:,:),pointer :: age => null()
    !*FD The age of a given ice layer, divided by \texttt{tim0}.

    integer, dimension(:,:),pointer :: mask => null()
    !*FD Set to zero for all points where $\mathtt{thck}=0$, otherwise non-zero.
    !*FD the non-zero points are numbered in sequence from the bottom left to the 
    !*FD top right, going along the rows.

    integer, dimension(:,:),pointer :: thkmask => null()
    !*FD see glide_mask.f90 for possible values

    integer :: totpts = 0
    !*FD The total number of points with non-zero thickness

    integer, dimension(4) :: dom   = 0      !*FD I have no idea what this is for.
    logical               :: empty = .true. !*FD I have no idea what this is for.

    real(dp) :: ivol, iarea !*FD ice volume and ice area

  end type glide_geometry

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_geomderv

    !*FD Holds the horizontal and temporal derivatives of the thickness and
    !*FD upper surface elevation, as well as the thickness on the staggered grid.

    real(dp),dimension(:,:),pointer :: dthckdew => null() !*FD E-W derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdew => null() !*FD E-W derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdns => null() !*FD N-S derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdns => null() !*FD N-S derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: dthckdtm => null() !*FD Temporal derivative of thickness.
    real(dp),dimension(:,:),pointer :: dusrfdtm => null() !*FD Temporal derivative of upper surface elevation.
    real(dp),dimension(:,:),pointer :: stagthck => null() !*FD Thickness averaged onto the staggered grid.

  end type glide_geomderv

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velocity

    !*FD Holds the velocity fields in 2D and 3D. At least some of these fields
    !*FD are stored on the displaced grid.

    real(dp),dimension(:,:,:),pointer :: uvel  => null() !*FD 3D $x$-velocity.
    real(dp),dimension(:,:,:),pointer :: vvel  => null() !*FD 3D $y$-velocity.
    real(dp),dimension(:,:,:),pointer :: wvel  => null() !*FD 3D $z$-velocity.
    real(dp),dimension(:,:,:),pointer :: wgrd  => null() !*FD 3D grid vertical velocity.
    real(dp),dimension(:,:)  ,pointer :: uflx  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: vflx  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: diffu => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: total_diffu => null() !*FD total diffusivity
    real(dp),dimension(:,:)  ,pointer :: ubas  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: ubas_tavg  => null()
    real(dp),dimension(:,:)  ,pointer :: vbas  => null() !*FD 
    real(dp),dimension(:,:)  ,pointer :: vbas_tavg  => null() 
    real(dp),dimension(:,:)  ,pointer :: bed_softness => null() !*FD bed softness parameter
    real(dp),dimension(:,:)  ,pointer :: btrc  => null() !*FD  basal traction
    real(dp),dimension(:,:)  ,pointer :: tau_x => null() !*FD basal shear stress, x-dir
    real(dp),dimension(:,:)  ,pointer :: tau_y => null() !*FD basal shear stress, y-dir
  end type glide_velocity

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_climate
     !*FD Holds fields used to drive the model
     real(sp),dimension(:,:),pointer :: acab     => null() !*FD Annual mass balance.
     real(sp),dimension(:,:),pointer :: acab_tavg     => null() !*FD Annual mass balance (time average).
     real(sp),dimension(:,:),pointer :: artm     => null() !*FD Annual mean air temperature
     real(sp),dimension(:,:),pointer :: lati     => null() !*FD Latitudes of model grid points
     real(sp),dimension(:,:),pointer :: loni     => null() !*FD Longitudes of model grid points
     real(sp),dimension(:,:),pointer :: calving  => null() !*FD Calving flux (scaled as mass balance, thickness, etc)
     real(sp) :: eus = 0.                                  !*FD eustatic sea level
  end type glide_climate

  type glide_temper

    !*FD Holds fields relating to temperature.

    real(dp),dimension(:,:,:),pointer :: temp => null() !*FD Three-dimensional temperature field.
    real(dp),dimension(:,:),  pointer :: bheatflx => null() !*FD basal heat flux
    real(dp),dimension(:,:,:),pointer :: flwa => null() !*FD Glenn's $A$.
    real(dp),dimension(:,:),  pointer :: bwat => null() !*FD Basal water depth
    real(dp),dimension(:,:),  pointer :: stagbwat => null() !*FD Basal water depth in velo grid
    real(dp),dimension(:,:),  pointer :: stagbtemp => null() !*FD Basal temperature on velo grid
    real(dp),dimension(:,:),  pointer :: bmlt => null() !*FD Basal melt-rate
    real(dp),dimension(:,:),  pointer :: bmlt_tavg => null() !*FD Basal melt-rate
    real(dp),dimension(:,:),  pointer :: bpmp => null() !*FD Basal pressure melting point
    real(dp),dimension(:,:),  pointer :: stagbpmp => null() !*FD Basal pressure melting point on velo grid
    
    integer  :: niter   = 0      !*FD
    real(sp) :: perturb = 0.0    !*FD
    real(sp) :: grid    = 0.0    !*FD
    integer  :: tpt     = 0      !*FD Pointer to time series data
    logical  :: first1  = .true. !*FD
    logical  :: newtemps = .false. !*FD new temperatures
  end type glide_temper

  type glide_lithot_type
     !*FD holds variables for temperature calculations in the lithosphere

     real(dp),dimension(:,:,:),pointer :: temp => null()    !*FD Three-dimensional temperature field.
     logical, dimension(:,:), pointer :: mask => null()     !*FD whether the point has been ice covered at some time

     integer :: num_dim = 1                                 !*FD either 1 or 3 for 1D/3D calculations

     ! The sparse matrix and linearised arrays
     type(sparse_matrix_type) :: fd_coeff, fd_coeff_slap
     integer :: all_bar_top
     real(dp), dimension(:), pointer :: rhs
     real(dp), dimension(:), pointer :: answer
     real(dp), dimension(:), pointer :: supd,diag,subd

     ! work arrays for solver
     real(dp), dimension(:), pointer :: rwork
     integer, dimension(:), pointer :: iwork
     integer mxnelt

     real(dp), dimension(:), pointer :: deltaz => null()    !*FD array holding grid spacing in z
     real(dp), dimension(:,:), pointer :: zfactors => null()!*FD array holding factors for finite differences of vertical diffu
     real(dp) :: xfactor,yfactor !*FD factors for finite differences of horizontal diffu


     real :: surft = 2.         !*FD surface temperature, used for calculating initial temperature distribution
     real :: mart  = 2.         !*FD sea floor temperature 
     integer :: nlayer = 20     !*FD number of layers in lithosphere
     real :: rock_base = -5000. !*FD depth below sea-level at which geothermal heat gradient is applied
     
     integer :: numt = 0        !*FD number time steps for spinning up GTHF calculations

     real(dp) :: rho_r = 3300.0d0 !*FD The density of lithosphere (kg m$^{-3}$)
     real(dp) :: shc_r = 1000.0d0 !*FD specific heat capcity of lithosphere (J kg$^{-1}$ K$^{-1}$)
     real(dp) :: con_r = 3.3d0    !*FD thermal conductivity of lithosphere (W m$^{-1}$ K$^{-1}$)

     real(dp) :: diffu = 0. !*FD diffusion coefficient

  end type glide_lithot_type

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_funits
    character(fname_length) :: sigfile=''                      !*FD sigma coordinates file
    character(fname_length) :: ncfile=''                       !*FD configuration file for netCDF I/O
    type(glimmer_nc_output),pointer :: out_first=>NULL()       !*FD first element of linked list defining netCDF outputs
    type(glimmer_nc_input), pointer :: in_first=>NULL()        !*FD first element of linked list defining netCDF inputs
  end type glide_funits

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_numerics

    !*FD Parameters relating to the model numerics.
    real(sp) :: tstart = 0.0      !*FD starting time
    real(sp) :: tend   = 20000.0  !*FD end time
    real(sp) :: time   =    0.0   !*FD main loop counter in years
    real(sp) :: tinc   =   20.0   !*FD time step of main loop in years 
    real(sp) :: ntem   =    1.0   !*FD temperature time step (multiplier of main time step)
    real(sp) :: nvel   =    1.0   !*FD velocity time step (multiplier of main time step)
    real(dp) :: alpha  =    0.5d0 !*FD richard suggests 1.5 - was a parameter in original
    real(dp) :: alphas =    0.5d0 !*FD was a parameter in the original
    real(dp) :: thklim =  100.0   
    real(dp) :: mlimit = -200.0d0
    real(dp) :: calving_fraction = 0.8d0
    real(dp) :: dew    =   20.0d3
    real(dp) :: dns    =   20.0d3
    real(dp) :: dt     =    0.0
    real(dp) :: dttem  =    0.0
    real(sp) :: nshlf  =    0.0

    integer  :: timecounter = 0   !*FD count time steps
    
    ! Vertical coordinate ---------------------------------------------------
                                                               
    real(dp),dimension(:),pointer :: sigma => null() !*FD Sigma values for 
                                                     !*FD vertical spacing of 
                                                     !*FD model levels
    integer :: profile_period = 100            !*FD profile frequency
    integer :: ndiag = 1000                    !*FD diagnostic frequency

  end type glide_numerics


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_velowk
    real(dp),dimension(:),  pointer :: depth    => null()
    real(dp),dimension(:),  pointer :: dupsw    => null()
    real(dp),dimension(:),  pointer :: depthw   => null()
    real(dp),dimension(:),  pointer :: suvel    => null()
    real(dp),dimension(:),  pointer :: svvel    => null()
    real(dp),dimension(:,:),pointer :: fslip    => null()
    real(dp),dimension(:,:),pointer :: dintflwa => null()
    real(dp),dimension(:),  pointer :: dups     => null()
    real(dp),dimension(4) :: fact
    real(dp),dimension(4) :: c    = 0.0
    real(dp) :: watwd  = 3.0d0
    real(dp) :: watct  = 10.0d0
    real(dp) :: trc0   = 0.0
    real(dp) :: trcmin = 0.0d0
    real(dp) :: marine = 1.0d0
    real(dp) :: trcmax = 10.0d0
    real(dp) :: btrac_const = 0.0d0
    real(dp) :: btrac_slope = 0.0d0
    real(dp) :: btrac_max = 0.d0
  end type glide_velowk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_pcgdwk
    real(dp),dimension(:),pointer :: pcgval  => null()
    real(dp),dimension(:),pointer :: rhsd    => null()
    real(dp),dimension(:),pointer :: answ    => null()
    integer, dimension(:),pointer :: pcgcol  => null()
    integer, dimension(:),pointer :: pcgrow  => null()
    integer, dimension(2)         :: pcgsize = 0
    real(dp),dimension(4)         :: fc      = 0.0
    real(dp),dimension(6)         :: fc2     = 0.0
    integer :: ct     = 0
  end type glide_pcgdwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_thckwk
     real(dp),dimension(:,:),  pointer :: oldthck   => null()
     real(dp),dimension(:,:),  pointer :: oldthck2  => null()
     real(dp),dimension(:,:),pointer :: float => null()
     real(dp),dimension(:,:,:),pointer :: olds      => null()
     integer  :: nwhich  = 2
     real(sp) :: oldtime = 0.0
     
     real(dp), dimension(:), pointer :: alpha => null()
     real(dp), dimension(:), pointer :: beta  => null()
     real(dp), dimension(:), pointer :: gamma => null()
     real(dp), dimension(:), pointer :: delta => null()

  end type glide_thckwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_tempwk
    real(dp),dimension(:,:,:),pointer :: inittemp => null()
    real(dp),dimension(:,:,:),pointer :: dissip   => null()
    real(dp),dimension(:,:,:),pointer :: initadvt => null()
    real(dp),dimension(:),    pointer :: dupa     => null()
    real(dp),dimension(:),    pointer :: dupb     => null()
    real(dp),dimension(:),    pointer :: dupc     => null()
    real(dp),dimension(:),    pointer :: c1       => null()
    real(dp),dimension(:,:),  pointer :: dups     => null()
    real(dp),dimension(:,:),  pointer :: wphi     => null()
    real(dp),dimension(:,:),  pointer :: bwatu    => null()
    real(dp),dimension(:,:),  pointer :: bwatv    => null()
    real(dp),dimension(:,:),  pointer :: fluxew   => null()
    real(dp),dimension(:,:),  pointer :: fluxns   => null()
    real(dp),dimension(:,:),  pointer :: bint     => null()
    real(dp),dimension(:,:),  pointer :: smth     => null()
    real(dp),dimension(:,:,:),pointer :: hadv_u   => null()
    real(dp),dimension(:,:,:),pointer :: hadv_v   => null()
    real(dp),dimension(4)             :: cons     = 0.0
    real(dp),dimension(4)             :: f        = 0.0
    real(dp),dimension(8)             :: c        = 0.0
    real(dp),dimension(2)             :: slide_f
    real(dp) :: noflow      = -1
    real(dp),dimension(2) :: advconst = 0.0
    real(dp) :: zbed        = 0.0
    real(dp) :: dupn        = 0.0
    real(dp) :: wmax        = 0.0
    real(dp) :: dt_wat      = 0.0
    real(dp) :: watvel      = 0.0
    integer  :: nwat        = 0
  end type glide_tempwk

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_paramets
    real(dp),dimension(5) :: bpar = (/ 2.0d0, 10.0d0, 10.0d0, 0.0d0, 1.0d0 /)
    real(dp) :: btrac_const = 0.d0 ! m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_slope = 0.0d0 ! Pa^{-1} (gets scaled during init)
    real(dp) :: btrac_max = 0.d0  !  m yr^{-1} Pa^{-1} (gets scaled during init)
    real(dp) :: geot   = -5.0d-2  ! W m^{-2}
    real(dp) :: fiddle = 3.0d0    ! -
    real(dp) :: hydtim = 1000.0d0 ! yr^{-1} converted to s^{-1} and scaled, 
                                  ! 0 if no drainage = 0.0d0 * tim0 / scyr
    real(dp) :: bwat_smooth = 0.01d0 ! basal water field smoothing strength
  end type glide_paramets

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type glide_prof_type
     integer :: geomderv
     integer :: hvelos
     integer :: ice_mask1
     integer :: temperature
     integer :: ice_evo
     integer :: ice_mask2
     integer :: isos_water
     integer :: isos
  end type glide_prof_type
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  type glide_phaml
    real(dp),dimension(:,:),pointer :: uphaml => null()
    real(dp),dimension(:,:),pointer :: init_phaml => null()
    real(dp),dimension(:,:),pointer :: rs_phaml => null()
    !maybe put the x/y vectors here too just for simplicity
  end type glide_phaml
  
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  type glide_global_type
    type(glide_general)  :: general
    type(glide_options)  :: options
    type(glide_geometry) :: geometry
    type(glide_geomderv) :: geomderv
    type(glide_velocity) :: velocity
    type(glide_climate)  :: climate
    type(glide_temper)   :: temper
    type(glide_lithot_type) :: lithot
    type(glide_funits)   :: funits
    type(glide_numerics) :: numerics
    type(glide_velowk)   :: velowk
    type(glide_pcgdwk)   :: pcgdwk
    type(glide_thckwk)   :: thckwk
    type(glide_tempwk)   :: tempwk
    type(glide_paramets) :: paramets
    type(glimmap_proj)   :: projection
    type(profile_type)   :: profile
    type(glide_prof_type) :: glide_prof
    type(isos_type)      :: isos
    type(glide_phaml)    :: phaml
  end type glide_global_type

contains

  subroutine glide_allocarr(model)
    
    !*FD Allocates the model arrays, and initialises some of them to zero.
    !*FD These are the arrays allocated, and their dimensions:
    !*FD
    !*FD In \texttt{model\%temper}:
    !*FD \begin{itemize}
    !*FD \item \texttt{temp(upn,0:ewn+1,0:nsn+1))}
    !*FD \item \texttt{bheatflx(ewn,nsn))}
    !*FD \item \texttt{flwa(upn,ewn,nsn))}
    !*FD \item \texttt{bwat(ewn,nsn))}
    !*FD \item \texttt{bmlt(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%velocity}:
    !*FD \begin{itemize}
    !*FD \item \texttt{uvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{vvel(upn,ewn-1,nsn-1))}
    !*FD \item \texttt{wvel(upn,ewn,nsn))}
    !*FD \item \texttt{wgrd(upn,ewn,nsn))}
    !*FD \item \texttt{uflx(ewn-1,nsn-1))}
    !*FD \item \texttt{vflx(ewn-1,nsn-1))}
    !*FD \item \texttt{diffu(ewn,nsn))}
    !*FD \item \texttt{btrc(ewn,nsn))}
    !*FD \item \texttt{ubas(ewn,nsn))}
    !*FD \item \texttt{vbas(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%climate}:
    !*FD \begin{itemize}
    !*FD \item \texttt{acab(ewn,nsn))}
    !*FD \item \texttt{artm(ewn,nsn))}
    !*FD \item \texttt{lati(ewn,nsn))}
    !*FD \item \texttt{loni(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%geomderv}:
    !*FD \begin{itemize}
    !*FD \item \texttt{dthckdew(ewn,nsn))}
    !*FD \item \texttt{dusrfdew(ewn,nsn))}
    !*FD \item \texttt{dthckdns(ewn,nsn))}
    !*FD \item \texttt{dusrfdns(ewn,nsn))}
    !*FD \item \texttt{dthckdtm(ewn,nsn))}
    !*FD \item \texttt{dusrfdtm(ewn,nsn))}
    !*FD \item \texttt{stagthck(ewn-1,nsn-1))}
    !*FD \end{itemize}
  
    !*FD In \texttt{model\%geometry}:
    !*FD \begin{itemize}
    !*FD \item \texttt{thck(ewn,nsn))}
    !*FD \item \texttt{usrf(ewn,nsn))}
    !*FD \item \texttt{lsrf(ewn,nsn))}
    !*FD \item \texttt{topg(ewn,nsn))}
    !*FD \item \texttt{mask(ewn,nsn))}
    !*FD \item \texttt{age(ewn,nsn))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%thckwk}:
    !*FD \begin{itemize}
    !*FD \item \texttt{olds(ewn,nsn,thckwk\%nwhich))}
    !*FD \end{itemize}

    !*FD In \texttt{model\%numerics}:
    !*FD \begin{itemize}
    !*FD \item \texttt{sigma(upn))}
    !*FD \end{itemize}

    use glimmer_log

    implicit none

    type(glide_global_type),intent(inout) :: model

    integer :: ewn,nsn,upn

    ! for simplicity, copy these values...

    ewn=model%general%ewn
    nsn=model%general%nsn
    upn=model%general%upn

    ! Allocate appropriately

    allocate(model%general%x0(ewn-1))!; model%general%x0 = 0.0
    allocate(model%general%y0(nsn-1))!; model%general%y0 = 0.0
    allocate(model%general%x1(ewn))!; model%general%x1 = 0.0
    allocate(model%general%y1(nsn))!; model%general%y1 = 0.0
    allocate(model%temper%temp(upn,0:ewn+1,0:nsn+1)); model%temper%temp = 0.0
    call coordsystem_allocate(model%general%ice_grid, upn, model%temper%flwa)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bheatflx)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bwat)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbwat)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbtemp)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bmlt)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bpmp)
    call coordsystem_allocate(model%general%ice_grid, model%temper%bmlt_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%temper%stagbpmp)

    allocate(model%lithot%temp(1:ewn,1:nsn,model%lithot%nlayer)); model%lithot%temp = 0.0
    call coordsystem_allocate(model%general%ice_grid, model%lithot%mask)

    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%uvel)
    call coordsystem_allocate(model%general%velo_grid, upn, model%velocity%vvel)
    call coordsystem_allocate(model%general%ice_grid, upn, model%velocity%wvel)
    call coordsystem_allocate(model%general%ice_grid, upn, model%velocity%wgrd)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%uflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vflx)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%diffu)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%total_diffu)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%bed_softness)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%btrc)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%ubas_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%vbas_tavg)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%tau_x)
    call coordsystem_allocate(model%general%velo_grid, model%velocity%tau_y)

    call coordsystem_allocate(model%general%ice_grid, model%climate%acab)
    call coordsystem_allocate(model%general%ice_grid, model%climate%acab_tavg)
    call coordsystem_allocate(model%general%ice_grid, model%climate%artm)
    call coordsystem_allocate(model%general%ice_grid, model%climate%lati)
    call coordsystem_allocate(model%general%ice_grid, model%climate%loni)
    call coordsystem_allocate(model%general%ice_grid, model%climate%calving)

    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdew)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dthckdns)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%dusrfdns)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dthckdtm)
    call coordsystem_allocate(model%general%ice_grid, model%geomderv%dusrfdtm)
    call coordsystem_allocate(model%general%velo_grid, model%geomderv%stagthck)
  
    call coordsystem_allocate(model%general%velo_grid, model%geometry%temporary0)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%temporary1)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thck)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%usrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%lsrf)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%topg)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%mask)
    call coordsystem_allocate(model%general%ice_grid, model%geometry%thkmask)
    call coordsystem_allocate(model%general%ice_grid, upn, model%geometry%age)

    allocate(model%thckwk%olds(ewn,nsn,model%thckwk%nwhich))
    model%thckwk%olds = 0.0d0
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck)
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%oldthck2)
    call coordsystem_allocate(model%general%ice_grid, model%thckwk%float)

    ! If we already have sigma, don't reallocate
    if (associated(model%numerics%sigma)) then
       if (size(model%numerics%sigma)/=upn) then
          call write_log('Wrong number of sigma levels given',GM_FATAL)
       end if
    else
       allocate(model%numerics%sigma(upn))
    endif

    ! allocate memory for sparse matrix
    allocate (model%pcgdwk%pcgrow(ewn*nsn*5))
    allocate (model%pcgdwk%pcgcol(ewn*nsn*5+2))
    allocate (model%pcgdwk%pcgval(ewn*nsn*5))
    allocate (model%pcgdwk%rhsd(ewn*nsn))
    allocate (model%pcgdwk%answ(ewn*nsn))

    ! allocate isostasy grids
    call isos_allocate(model%isos,ewn,nsn)
    
    !allocate phaml variables
    call coordsystem_allocate(model%general%ice_grid, model%phaml%init_phaml)
    call coordsystem_allocate(model%general%ice_grid, model%phaml%rs_phaml)
    call coordsystem_allocate(model%general%ice_grid, model%phaml%uphaml)
  end subroutine glide_allocarr

  subroutine glide_deallocarr(model)
    !*FD deallocate model arrays
    implicit none
    type(glide_global_type),intent(inout) :: model

    deallocate(model%general%x0) 
    deallocate(model%general%y0) 
    deallocate(model%general%x1) 
    deallocate(model%general%y1) 

    deallocate(model%temper%temp)
    deallocate(model%temper%flwa)
    deallocate(model%temper%bheatflx)
    deallocate(model%temper%bwat)
    deallocate(model%temper%stagbwat)
    deallocate(model%temper%stagbtemp)
    deallocate(model%temper%bmlt)
    deallocate(model%temper%bmlt_tavg)
    deallocate(model%temper%bpmp)
    deallocate(model%temper%stagbpmp)

    deallocate(model%lithot%temp)
    deallocate(model%lithot%mask)

    deallocate(model%velocity%uvel)
    deallocate(model%velocity%vvel)
    deallocate(model%velocity%wvel)
    deallocate(model%velocity%wgrd)
    deallocate(model%velocity%uflx)
    deallocate(model%velocity%vflx)
    deallocate(model%velocity%diffu)
    deallocate(model%velocity%total_diffu)
    deallocate(model%velocity%bed_softness)
    deallocate(model%velocity%btrc)
    deallocate(model%velocity%ubas)
    deallocate(model%velocity%ubas_tavg)
    deallocate(model%velocity%vbas)
    deallocate(model%velocity%vbas_tavg)
    deallocate(model%velocity%tau_x)
    deallocate(model%velocity%tau_y)

    deallocate(model%climate%acab)
    deallocate(model%climate%acab_tavg)
    deallocate(model%climate%artm)
    deallocate(model%climate%lati)
    deallocate(model%climate%loni)

    deallocate(model%geomderv%dthckdew)
    deallocate(model%geomderv%dusrfdew)
    deallocate(model%geomderv%dthckdns)
    deallocate(model%geomderv%dusrfdns)
    deallocate(model%geomderv%dthckdtm)
    deallocate(model%geomderv%dusrfdtm)
    deallocate(model%geomderv%stagthck)
  
    deallocate(model%geometry%temporary0)
    deallocate(model%geometry%temporary1)
    deallocate(model%geometry%thck)
    deallocate(model%geometry%usrf)
    deallocate(model%geometry%lsrf)
    deallocate(model%geometry%topg)
    deallocate(model%geometry%mask)
    deallocate(model%geometry%thkmask)
    deallocate(model%geometry%age)
    deallocate(model%thckwk%olds)
    deallocate(model%thckwk%oldthck)
    deallocate(model%thckwk%oldthck2)
    deallocate(model%thckwk%float)
    deallocate(model%numerics%sigma)
    
    deallocate(model%pcgdwk%pcgrow,model%pcgdwk%pcgcol,model%pcgdwk%pcgval,model%pcgdwk%rhsd,model%pcgdwk%answ)

    ! allocate isostasy grids
    call isos_deallocate(model%isos)

    !deallocate phaml variables
    deallocate(model%phaml%init_phaml)
    deallocate(model%phaml%rs_phaml)    
    deallocate(model%phaml%uphaml)
    
  end subroutine glide_deallocarr

  ! some accessor functions
  function get_dew(model)
    !*FD return scaled x node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dew
    type(glide_global_type) :: model

    get_dew = model%numerics%dew * len0
  end function get_dew

  function get_dns(model)
    !*FD return scaled y node spacing
    use glimmer_paramets, only : len0
    implicit none
    real(dp) :: get_dns
    type(glide_global_type) :: model

    get_dns = model%numerics%dns * len0
  end function get_dns

  function get_tstart(model)
    !*FD return start time
    implicit none
    real(sp) :: get_tstart
    type(glide_global_type) :: model
    
    get_tstart = model%numerics%tstart
  end function get_tstart

  function get_tend(model)
    !*FD return end time
    implicit none
    real(sp) :: get_tend
    type(glide_global_type) :: model
    
    get_tend = model%numerics%tend
  end function get_tend

  function get_tinc(model)
    !*FD return time increment
    implicit none
    real(sp) :: get_tinc
    type(glide_global_type) :: model
    
    get_tinc = model%numerics%tinc
  end function get_tinc

  function get_ewn(model)
    !*FD get number of nodes in x dir
    implicit none
    integer get_ewn
    type(glide_global_type) :: model

    get_ewn = model%general%ewn
  end function get_ewn

  function get_nsn(model)
    !*FD get number of nodes in y dir
    implicit none
    integer get_nsn
    type(glide_global_type) :: model

    get_nsn = model%general%nsn
  end function get_nsn
  
  subroutine set_time(model,time)
    !*FD Set the model time counter --- useful for
    !*FD fractional year output
    implicit none
    type(glide_global_type) :: model
    real :: time

    model%numerics%time=time
  end subroutine set_time

end module glide_types

