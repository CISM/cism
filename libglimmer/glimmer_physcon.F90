
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_physcon.f90 - part of the GLIMMER ice model      + 
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

module physcon

  !*FD Contains physical constants required by the ice model.

  use glimmer_global, only : dp,sp

  implicit none
  
  save

  real(dp),parameter :: scyr = 31556926.0        !*FD Number of seconds in a year (s). Note
                                                 !*FD that this is for a 365.242 day year, and might
                                                 !*FD need changing.

  real(dp),parameter :: pi = 3.1415926535897d0   !*FD Value of $\pi$.

  real(dp),parameter :: rhoi = 910.0             !*FD The density of ice (kg m$^{-3}$)
  real(dp),parameter :: rhom = 3300.0d0          !*FD The density of magma(?) (kg m$^{-3}$)

  real(dp),parameter :: rhoo = 1028.0d0          !*FD The density of the ocean (kg m$^{-3}$)
  real(dp),parameter :: rhow = 1000.0d0          !*FD The density of fresh water (kg m$^{-3}$)
  real(dp),parameter :: f = - rhoo / rhoi

  real(dp),parameter :: grav = 9.81              !*FD The acceleration due to gravity (m s$^{-2}$)

  integer, parameter :: gn = 3                   !*FD The power dependency of Glenn's flow law.

  real(dp),parameter :: arrmlh = 1.733d3         !*FD Constant of proportionality in Arrhenius relation
                                                 !*FD in \texttt{patebudd}, for $T^{*}\geq263$K.
                                                 !*FD (Pa$^{-3}$ s$^{-1}$) 
  real(dp),parameter :: arrmll = 3.613d-13       !*FD Constant of proportionality in Arrhenius relation
                                                 !*FD in \texttt{patebudd}, for $T^{*}<263$K.
                                                 !*FD (Pa$^{-3}$ s$^{-1}$) 
  real(dp),parameter :: gascon = 8.314d0         !*FD The gas ideal constant $R$ (J mol$^{-1}$ K$^{-1}$)
  real(dp),parameter :: actenh = 139.0d3         !*FD Activation energy in Glenn's flow law for $T^{*}\geq263$K. (J mol$^{-1}$)
  real(dp),parameter :: actenl = 60.0d3          !*FD Activation energy in Glenn's flow law for $T^{*}<263$K. (J mol$^{-1}$)

  real(dp),parameter :: shci = 2009.0d0          !*FD Specific heat capacity of ice (J kg$^{-1}$ K$^{-1}$)
  real(dp),parameter :: lhci = 335.0d3           !*FD Latent heat of melting of ice (J kg$^{-1}$) 
  real(dp),parameter :: coni = 2.1d0             !*FD Thermal conductivity of ice (W m$^{-1}$ K$^{-1}$)

  real(dp),parameter :: pmlt = 9.7456d-8         !*FD Factor for dependence of melting point on pressure (K Pa$^{-1}$)
  real(dp),parameter :: trpt = 273.15d0          !*FD Triple point of water (K)

end module physcon
