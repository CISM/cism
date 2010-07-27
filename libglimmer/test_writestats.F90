!> test glimmer_writestats module
!!
!! \author Magnus Hagdorn
!! \date April 2009

program test_writestats
  use glimmer_writestats
  use glimmer_global, only : dp
  implicit none

  call glimmer_write_stats("results","model.conf",1000.2_dp)
end program test_writestats
