module nan_mod

! Set parameter for the floating point flag "nan" not-a-number
!
! Based on the similar module in CESM's CLM & CAM
!
  use glimmer_global, only : dp

  implicit none
  save

#ifdef __PGI
! quiet nan for portland group compilers
  real(dp), parameter :: NaN = O'0777700000000000000000'
#else
! signaling nan otherwise
  real(dp), parameter :: NaN = O'0777610000000000000000'
#endif

end module nan_mod
