
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

!> Contains physical constants required by the ice model.
module glimmer_physcon


  use glimmer_global, only : dp,sp

  implicit none
  
  save

  real(dp),parameter :: scyr = 31556926.0        !< Number of seconds in a year (s). Note that this is for a 365.242 day year, and might need changing.

  real(dp),parameter :: pi = 3.1415926535897d0   !< Value of \f$\pi\f$.

  real(dp),parameter :: rhoi = 910.0             !< The density of ice (kg m<SUP>-3</SUP>)
  real(dp),parameter :: rhom = 3300.0d0          !< The density of magma(?) (kg m<SUP>-3</SUP>)

  real(dp),parameter :: rhoo = 1028.0d0          !< The density of the ocean (kg m<SUP>-3</SUP>)
  real(dp),parameter :: rhow = 1000.0d0          !< The density of fresh water (kg m<SUP>-3</SUP>)
  real(dp),parameter :: f = - rhoo / rhoi

  real(dp),parameter :: grav = 9.81              !< The acceleration due to gravity (m s<SUP>-2</SUP>)

  integer, parameter :: gn = 3                   !< The power dependency of Glenn's flow law.

  real(dp),parameter :: arrmlh = 1.733d3         !< Constant of proportionality in Arrhenius relation
                                                 !< in \texttt{patebudd}, for \f$T^{*}\geq263\f$K.
                                                 !< (Pa<SUP>-3</SUP> s<SUP>-1</SUP>) 
  real(dp),parameter :: arrmll = 3.613d-13       !< Constant of proportionality in Arrhenius relation
                                                 !< in \texttt{patebudd}, for \f$T^{*}<263\f$K.
                                                 !< (Pa<SUP>-3</SUP> s<SUP>-1</SUP>) 
  real(dp),parameter :: gascon = 8.314d0         !< The gas ideal constant \f$R\f$ (J mol<SUP>-1</SUP> K<SUP>-1</SUP>)
  real(dp),parameter :: actenh = 139.0d3         !< Activation energy in Glenn's flow law for \f$T^{*}\geq263\f$K. (J mol<SUP>-1</SUP>)
  real(dp),parameter :: actenl = 60.0d3          !< Activation energy in Glenn's flow law for \f$T^{*}<263\f$K. (J mol<SUP>-1</SUP>)

  real(dp),parameter :: shci = 2009.0d0          !< Specific heat capacity of ice (J kg<SUP>-1</SUP> K<SUP>-1</SUP>)
  real(dp),parameter :: lhci = 335.0d3           !< Latent heat of melting of ice (J kg<SUP>-1</SUP>) 
  real(dp),parameter :: coni = 2.1d0             !< Thermal conductivity of ice (W m<SUP>-1</SUP> K<SUP>-1</SUP>)

  real(dp),parameter :: pmlt = 9.7456d-8         !< Factor for dependence of melting point on pressure (K Pa<SUP>-1</SUP>)
  real(dp),parameter :: trpt = 273.15d0          !< Triple point of water (K)

end module glimmer_physcon
