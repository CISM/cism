module parallel
  use netcdf

  implicit none

  integer,parameter :: lhalo = 0
  integer,parameter :: uhalo = 0
  logical,parameter :: main_task = .true.
  integer,parameter :: this_rank = 0
  integer,parameter :: tasks = 1

  ! distributed grid
  integer,save :: global_ewn,global_nsn,local_ewn,local_nsn,own_ewn,own_nsn
  integer,save :: ewlb,ewub,nslb,nsub

  ! global IDs
  integer,parameter :: ProcsEW = 1

  ! JEFF Declarations for undistributed variables on main_task.
  ! Later move to separate module?  These are only temporary until code is completely distributed.
  real(8),dimension(:,:,:),allocatable :: gathered_efvs  ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:,:),allocatable :: gathered_efvs2  ! Variable for testing that scatter/gather are inverses
  real(8),dimension(:,:,:),allocatable :: gathered_uvel  ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:,:),allocatable :: gathered_vvel  ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:),allocatable :: gathered_uflx    ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:),allocatable :: gathered_vflx    ! Output var from glam_velo_fordsiapstr(), used often
  real(8),dimension(:,:,:),allocatable :: gathered_velnorm  ! Variable calculated in run_ho_diagnostic(), is this used?
  real(8),dimension(:,:),allocatable :: gathered_thck    ! Used in horizontal_remap_in()
  real(8),dimension(:,:),allocatable :: gathered_stagthck ! Used in horizontal_remap_in()
  real(4),dimension(:,:),allocatable :: gathered_acab    ! Used in horizontal_remap_in()
  real(8),dimension(:,:,:),allocatable :: gathered_temp  ! Used in horizontal_remap_in()
  real(8),dimension(:,:),allocatable :: gathered_dusrfdew  ! Used in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_dusrfdns  ! Used in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_dthckdew  ! Used in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_dthckdns  ! Used in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauxx   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauyy   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauxy   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauscalar   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauxz   ! Calculated in glide_stress()
  real(8),dimension(:,:,:),allocatable :: gathered_tauyz   ! Calculated in glide_stress()
  real(8),dimension(:,:),allocatable :: gathered_topg  ! Bedrock topology, Used in glide_set_mask()
  integer,dimension(:,:),allocatable :: gathered_thkmask  ! Calculated in glide_set_mask()
  real(8),dimension(:,:),allocatable :: gathered_marine_bc_normal  ! Calculated in glide_marine_margin_normal()
  real(8),dimension(:,:,:),allocatable :: gathered_surfvel   ! Used in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_gline_flux   ! Calculated in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_ubas   ! Used in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_vbas   ! Used in calc_gline_flux()
  real(8),dimension(:,:),allocatable :: gathered_relx   ! Used in glide_marinlim()
  real(8),dimension(:,:,:),allocatable :: gathered_flwa   ! Used in glide_marinlim()
  real(4),dimension(:,:),allocatable :: gathered_calving   ! Used in glide_marinlim()
  real(4),dimension(:,:),allocatable :: gathered_backstress   ! Used in glide_marinlim()
  real(8),dimension(:,:),allocatable :: gathered_usrf   ! Used in glide_marinlim()
  logical,dimension(:,:),allocatable :: gathered_backstressmap ! Used in glide_marinlim()
  real(8),dimension(:,:),allocatable :: gathered_tau_x   ! Calculated in calc_basal_shear()
  real(8),dimension(:,:),allocatable :: gathered_tau_y   ! Calculated in calc_basal_shear()
  real(8),dimension(:,:),allocatable :: gathered_lsrf   ! Used in glide_marinlim()

  integer,parameter :: staggered_lhalo = lhalo
  integer,parameter :: staggered_uhalo = 0

  interface broadcast
     module procedure broadcast_character
     module procedure broadcast_integer
     module procedure broadcast_integer_1d
     module procedure broadcast_logical
     module procedure broadcast_real8_1d
  end interface

  interface distributed_gather_var
     module procedure distributed_gather_var_integer_2d
     module procedure distributed_gather_var_logical_2d
     module procedure distributed_gather_var_real4_2d
     module procedure distributed_gather_var_real4_3d
     module procedure distributed_gather_var_real8_2d
     module procedure distributed_gather_var_real8_3d
  end interface

  interface distributed_get_var
     module procedure distributed_get_var_integer_2d
     module procedure distributed_get_var_real4_1d
     module procedure distributed_get_var_real4_2d
     module procedure distributed_get_var_real8_2d
     module procedure distributed_get_var_real8_3d
  end interface

  interface distributed_print
     ! Gathers a distributed variable and writes to file
     module procedure distributed_print_integer_2d
     module procedure distributed_print_real8_2d
     module procedure distributed_print_real8_3d
  end interface

  interface distributed_put_var
     module procedure distributed_put_var_integer_2d
     module procedure distributed_put_var_real4_1d
     module procedure distributed_put_var_real4_2d
     module procedure distributed_put_var_real8_2d
     module procedure distributed_put_var_real8_3d
     module procedure parallel_put_var_real8
  end interface

  interface parallel_def_var
     module procedure parallel_def_var_dimids
     module procedure parallel_def_var_nodimids
  end interface

  interface parallel_get_att
     module procedure parallel_get_att_character
     module procedure parallel_get_att_real4
     module procedure parallel_get_att_real4_1d
     module procedure parallel_get_att_real8
     module procedure parallel_get_att_real8_1d
  end interface

  interface distributed_scatter_var
     module procedure distributed_scatter_var_integer_2d
     module procedure distributed_scatter_var_logical_2d
     module procedure distributed_scatter_var_real4_2d
     module procedure distributed_scatter_var_real4_3d
     module procedure distributed_scatter_var_real8_2d
     module procedure distributed_scatter_var_real8_3d
  end interface

  interface parallel_get_var
     module procedure parallel_get_var_integer_1d
     module procedure parallel_get_var_real4_1d
  end interface

  interface parallel_halo
     module procedure parallel_halo_integer_2d
     module procedure parallel_halo_logical_2d
     module procedure parallel_halo_real4_2d
     module procedure parallel_halo_real8_2d
     module procedure parallel_halo_real8_3d
  end interface

  interface parallel_halo_verify
     module procedure parallel_halo_verify_integer_2d
     module procedure parallel_halo_verify_real8_2d
     module procedure parallel_halo_verify_real8_3d
  end interface

  interface parallel_print
     module procedure parallel_print_integer_2d
     module procedure parallel_print_real8_2d
     module procedure parallel_print_real8_3d
  end interface

  interface parallel_put_att
     module procedure parallel_put_att_character
     module procedure parallel_put_att_real4
     module procedure parallel_put_att_real4_1d
     module procedure parallel_put_att_real8
     module procedure parallel_put_att_real8_1d
  end interface

  interface parallel_put_var
     module procedure parallel_put_var_real4
     module procedure parallel_put_var_real8_1d
  end interface

contains

  subroutine broadcast_character(c)
    implicit none
    character(len=*) :: c
  end subroutine

  subroutine broadcast_integer(i)
    implicit none
    integer :: i
  end subroutine

  subroutine broadcast_integer_1d(a)
    implicit none
    integer,dimension(:) :: a
  end subroutine

  subroutine broadcast_logical(l)
    implicit none
    logical :: l
  end subroutine

  subroutine broadcast_real8_1d(a)
    implicit none
    real(8),dimension(:) :: a
  end subroutine

  function distributed_get_var_integer_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_get_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    ! begin
    if (main_task) distributed_get_var_integer_2d = nf90_get_var(ncid,varid,values(:,:),start)
  end function

  function distributed_get_var_real4_1d(ncid,varid,values,start)
    implicit none
    integer :: distributed_get_var_real4_1d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:) :: values

    ! begin
    if (main_task) distributed_get_var_real4_1d = nf90_get_var(ncid,varid,values(:),start)
  end function

  function distributed_get_var_real4_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_get_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:,:) :: values

    ! begin
    if (main_task) distributed_get_var_real4_2d = nf90_get_var(ncid,varid,values(:,:),start)
  end function

  function distributed_get_var_real8_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_get_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:) :: values

    ! begin
    if (main_task) distributed_get_var_real8_2d = nf90_get_var(ncid,varid,values(:,:),start)
  end function

  function distributed_get_var_real8_3d(ncid,varid,values,start)
    implicit none
    integer :: distributed_get_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:,:) :: values

    ! begin
    if (main_task) distributed_get_var_real8_3d = nf90_get_var(ncid,varid,values(:,:,:),start)
  end function

  subroutine distributed_grid(ewn,nsn)
    implicit none
    integer :: ewn,nsn
    ! begin
    global_ewn = ewn
    global_nsn = nsn

!JEFFCHECK Should these be performed even if tasks==1?    if (tasks >= 2) then
!Feb11 Bug, in simple_glide.F90, ewn and nsn set to local_ewn, _nsn which are zero.
		ewlb = 1
		ewub = global_ewn
		local_ewn = ewub-ewlb+1
		own_ewn = local_ewn-lhalo-uhalo
		ewn = local_ewn

		nslb = 1
		nsub = global_nsn
		local_nsn = nsub-nslb+1
		own_nsn = local_nsn-lhalo-uhalo
		nsn = local_nsn
!    endif
  end subroutine

  function distributed_execution()
     ! Returns if running distributed or not.
     logical distributed_execution

     distributed_execution = .false.
  end function

  subroutine distributed_gather_var_integer_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    implicit none
    integer,dimension(:,:),intent(in) :: values
    integer,dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

    allocate(global_values(size(values,1), size(values,2)))

    global_values(:,:) = values(:,:)
  end subroutine

  subroutine distributed_gather_var_logical_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    implicit none
    logical,dimension(:,:),intent(in) :: values
    logical,dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

    allocate(global_values(size(values,1), size(values,2)))

    global_values(:,:) = values(:,:)
  end subroutine

  subroutine distributed_gather_var_real4_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    implicit none
    real(4),dimension(:,:),intent(in) :: values
    real(4),dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

    allocate(global_values(size(values,1), size(values,2)))

    global_values(:,:) = values(:,:)
  end subroutine

  subroutine distributed_gather_var_real4_3d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    implicit none
    real(4),dimension(:,:,:),intent(in) :: values
    real(4),dimension(:,:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

    allocate(global_values(size(values,1), size(values,2), size(values,3)))

    global_values(:,:,:) = values(:,:,:)
  end subroutine

  subroutine distributed_gather_var_real8_2d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    implicit none
    real(8),dimension(:,:),intent(in) :: values
    real(8),dimension(:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

    allocate(global_values(size(values,1), size(values,2)))

    global_values(:,:) = values(:,:)
  end subroutine

  subroutine distributed_gather_var_real8_3d(values, global_values)
    ! JEFF Gather a distributed variable back to main_task node
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task will store the variable.
    ! If global_values is allocated, then it will be deallocated and reallocated.  It will be unused on other nodes.
    implicit none
    real(8),dimension(:,:,:),intent(in) :: values
    real(8),dimension(:,:,:),allocatable,intent(inout) :: global_values

    if (allocated(global_values)) then
       deallocate(global_values)
    endif

    allocate(global_values(size(values,1), size(values,2), size(values,3)))

    global_values(:,:,:) = values(:,:,:)
  end subroutine

  function distributed_owner(ew,ewn,ns,nsn)
    implicit none
    logical :: distributed_owner
    integer :: ew,ewn,ns,nsn
    ! begin
    distributed_owner = .true.
  end function

  subroutine distributed_print_integer_2d(name,values)
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k

    write(ts,'(i3.3)') tasks
    open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
    if (size(values,1)<local_ewn) then
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    else
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    end if
    close(u)
  end subroutine

  subroutine distributed_print_real8_2d(name,values)
    implicit none
    character(*) :: name
    real(8),dimension(:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k

    write(ts,'(i3.3)') tasks
    open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
    if (size(values,1)<local_ewn) then
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    else
       do j = lbound(values,2),ubound(values,2)
          do i = lbound(values,1),ubound(values,1)
             write(u,*) j,i,values(i,j)
          end do
          write(u,'()')
       end do
    end if
    close(u)
  end subroutine

  subroutine distributed_print_real8_3d(name,values)
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: values

    integer,parameter :: u = 33
    character(3) :: ts
    integer :: i,ierror,j,k

    write(ts,'(i3.3)') tasks
    open(unit=u,file=name//ts//".txt",form="formatted",status="replace")
    if (size(values,2)<local_ewn) then
       do j = lbound(values,3),ubound(values,3)
          do i = lbound(values,2),ubound(values,2)
             write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
          end do
          write(u,'()')
       end do
    else
       do j = lbound(values,3),ubound(values,3)
          do i = lbound(values,2),ubound(values,2)
             write(u,'(2i6,100g15.5e3)') j,i,values(:,i,j)
          end do
          write(u,'()')
       end do
    end if
    close(u)
  end subroutine

  function distributed_put_var_integer_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_integer_2d,ncid,varid
    integer,dimension(:) :: start
    integer,dimension(:,:) :: values

    ! begin
    if (main_task) distributed_put_var_integer_2d = nf90_put_var(ncid,varid,values,start)
    call broadcast(distributed_put_var_integer_2d)
  end function

  function distributed_put_var_real4_1d(ncid,varid,values)
    implicit none
    integer :: distributed_put_var_real4_1d,ncid,varid
    real(4),dimension(:) :: values

    ! begin
    if (main_task) distributed_put_var_real4_1d = nf90_put_var(ncid,varid,values)
    call broadcast(distributed_put_var_real4_1d)
  end function

  function distributed_put_var_real4_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real4_2d,ncid,varid
    integer,dimension(:) :: start
    real(4),dimension(:,:) :: values

    ! begin
    if (main_task) distributed_put_var_real4_2d = nf90_put_var(ncid,varid,values,start)
    call broadcast(distributed_put_var_real4_2d)
  end function

  function distributed_put_var_real8_2d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real8_2d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:) :: values

    ! begin
    if (main_task) distributed_put_var_real8_2d = nf90_put_var(ncid,varid,values,start)
    call broadcast(distributed_put_var_real8_2d)
  end function

  function distributed_put_var_real8_3d(ncid,varid,values,start)
    implicit none
    integer :: distributed_put_var_real8_3d,ncid,varid
    integer,dimension(:) :: start
    real(8),dimension(:,:,:) :: values

    ! begin
    if (main_task) distributed_put_var_real8_3d = nf90_put_var(ncid,varid,values,start)
    call broadcast(distributed_put_var_real8_3d)
  end function

  subroutine distributed_scatter_var_integer_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    integer,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    integer,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(:,:) = global_values(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine

  subroutine distributed_scatter_var_logical_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    logical,dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    logical,dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(:,:) = global_values(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine

  subroutine distributed_scatter_var_real4_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    real(4),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(4),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(:,:) = global_values(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine

  subroutine distributed_scatter_var_real4_3d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    real(4),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(4),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(:,:,:) = global_values(:,:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine

  subroutine distributed_scatter_var_real8_2d(values, global_values)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    real(8),dimension(:,:),intent(inout) :: values  ! populated from values on main_task
    real(8),dimension(:,:),allocatable,intent(inout) :: global_values  ! only used on main_task

    values(:,:) = global_values(:,:)

    deallocate(global_values)
    ! automatic deallocation
  end subroutine

  subroutine distributed_scatter_var_real8_3d(values, global_values, deallocflag)
    ! JEFF Scatter a variable on the main_task node back to the distributed
    ! values = local portion of distributed variable
    ! global_values = reference to allocateable array into which the main_task holds the variable.
    ! global_values is deallocated at the end.
    implicit none
    real(8),dimension(:,:,:),intent(inout) :: values  ! populated from values on main_task
    real(8),dimension(:,:,:),allocatable,intent(inout) :: global_values  ! only used on main_task
    logical,optional :: deallocflag
    logical :: deallocmem

    if (present(deallocflag)) then
       deallocmem = deallocflag
    else
       deallocmem = .true.
    endif

    ! begin
    values(:,:,:) = global_values(:,:,:)

    if (deallocmem) deallocate(global_values)
    ! automatic deallocation
  end subroutine

  subroutine global_sum(x,y)
    implicit none
    real(8) :: x,y
  end subroutine

  subroutine not_parallel(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "WARNING: not parallel in ",file," at line ",line
  end subroutine

  subroutine parallel_barrier
    implicit none
  end subroutine

  function parallel_boundary(ew,ewn,ns,nsn)
    implicit none
    logical :: parallel_boundary
    integer :: ew,ewn,ns,nsn
    ! begin
    parallel_boundary = (ew==1.or.ew==ewn.or.ns==1.or.ns==nsn)
  end function

  function parallel_close(ncid)
    implicit none
    integer :: ncid,parallel_close
    ! begin
    if (main_task) parallel_close = nf90_close(ncid)
    call broadcast(parallel_close)
  end function

  function parallel_create(path,cmode,ncid)
    implicit none
    integer :: cmode,ncid,parallel_create
    character(len=*) :: path
    ! begin
    if (main_task) parallel_create = nf90_create(path,cmode,ncid)
    call broadcast(parallel_create)
    call broadcast(ncid)
  end function

  function parallel_def_dim(ncid,name,len,dimid)
    use netcdf
    implicit none
    integer :: dimid,len,ncid,parallel_def_dim
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_dim = nf90_def_dim(ncid,name,len,dimid)
    call broadcast(parallel_def_dim)
    call broadcast(dimid)
  end function

  function parallel_def_var_dimids(ncid,name,xtype,dimids,varid)
    implicit none
    integer :: ncid,parallel_def_var_dimids,varid,xtype
    integer,dimension(:) :: dimids
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_var_dimids = &
         nf90_def_var(ncid,name,xtype,dimids,varid)
    call broadcast(parallel_def_var_dimids)
    call broadcast(varid)
  end function

  function parallel_def_var_nodimids(ncid,name,xtype,varid)
    implicit none
    integer :: ncid,parallel_def_var_nodimids,varid,xtype
    character(len=*) :: name
    ! begin
    if (main_task) parallel_def_var_nodimids = &
         nf90_def_var(ncid,name,xtype,varid)
    call broadcast(parallel_def_var_nodimids)
    call broadcast(varid)
  end function

  function parallel_enddef(ncid)
    implicit none
    integer :: ncid,parallel_enddef
    ! begin
    if (main_task) parallel_enddef = nf90_enddef(ncid)
    call broadcast(parallel_enddef)
  end function

  subroutine parallel_finalise
    implicit none
  end subroutine

  function parallel_get_att_character(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_character,varid
    character(len=*) :: name,values
    ! begin
    if (main_task) parallel_get_att_character = &
         nf90_get_att(ncid,varid,name,values)
    call broadcast(parallel_get_att_character)
    call broadcast(values)
  end function

  function parallel_get_att_real4(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4,varid
    character(len=*) :: name
    real(4) :: values
    ! begin
    if (main_task) parallel_get_att_real4 = &
         nf90_get_att(ncid,varid,name,values)
  end function

  function parallel_get_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real4_1d,varid
    character(len=*) :: name
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real4_1d = &
         nf90_get_att(ncid,varid,name,values)
  end function

  function parallel_get_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8,varid
    character(len=*) :: name
    real(8) :: values
    ! begin
    if (main_task) parallel_get_att_real8 = &
         nf90_get_att(ncid,varid,name,values)
  end function

  function parallel_get_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_get_att_real8_1d,varid
    character(len=*) :: name
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_att_real8_1d = &
         nf90_get_att(ncid,varid,name,values)
  end function

  function parallel_get_var_integer_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_integer_1d,varid
    integer,dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_integer_1d = &
         nf90_get_var(ncid,varid,values)
  end function

  function parallel_get_var_real4_1d(ncid,varid,values)
    implicit none
    integer :: ncid,parallel_get_var_real4_1d,varid
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_get_var_real4_1d = &
         nf90_get_var(ncid,varid,values)
  end function

  function parallel_globalID(locns, locew, upstride)
    ! Returns a unique ID for a given row and column reference that is identical across all processors.
    ! For instance if Proc 2: (17,16) is the same global cell as Proc 3: (17,1), then the globalID will be the same for both.
    ! These IDs are spaced upstride apart.  upstride = number of vertical layers.  Typically (upn) + number of ghost layers (2 = top and bottom)
    integer,intent(IN) :: locns, locew, upstride
    integer :: parallel_globalID
    ! locns is local NS (row) grid index
    ! locew is local EW (col) grid index
    integer :: global_row, global_col, global_ID
    character(len=40) :: local_coord

    global_row = (locns - uhalo) + this_rank/ProcsEW * own_nsn
    	! Integer division required for this_rank/ProcsEW
    global_col = (locew - lhalo) + mod(this_rank, ProcsEW) * own_ewn
        ! There are ProcsEW processors per row.

    global_ID = ((global_row - 1) * global_ewn + (global_col - 1)) * upstride + 1

    ! Testing Code
    ! write(local_coord, "A13,I10.1,A2,I10.1,A1") " (NS, EW) = (", locns, ", ", locew, ")"
	! write(*,*) "Processor reference ", this_rank, local_coord, " globalID = ", global_ID

	!return value
	parallel_globalID = global_ID
  end function

  subroutine parallel_halo_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_real8_2d(a)
    implicit none
    real(8),dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_real8_3d(a)
    implicit none
    real(8),dimension(:,:,:) :: a
  end subroutine

  subroutine parallel_halo_verify_integer_2d(a)
    implicit none
    integer,dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_logical_2d(a)
    implicit none
    logical,dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_real4_2d(a)
    implicit none
    real(4),dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_temperature(a)
    !JEFF This routine is for updating the halo for the variable model%temper%temp.
    ! This variable is two larger in each dimension, because of the current advection code.
    ! Per Bill L, we will remove this difference when we update the remapping code.
    implicit none
    real(8),dimension(:,:,:) :: a

    if (tasks >= 2) then
		! Do nothing
	endif
  end subroutine

  subroutine parallel_halo_verify_real8_2d(a)
    implicit none
    real(8),dimension(:,:) :: a
  end subroutine

  subroutine parallel_halo_verify_real8_3d(a)
    implicit none
    real(8),dimension(:,:,:) :: a
  end subroutine

  subroutine parallel_initialise
    implicit none
  end subroutine

  subroutine parallel_print_integer_2d(name,a)
    implicit none
    character(*) :: name
    integer,dimension(:,:) :: a
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(a,2)
       do i = 1,size(a,1)
          write(u,*) j,i,a(i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine

  function parallel_inq_attname(ncid,varid,attnum,name)
    implicit none
    integer :: attnum,ncid,parallel_inq_attname,varid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_attname = &
         nf90_inq_attname(ncid,varid,attnum,name)
    call broadcast(parallel_inq_attname)
    call broadcast(name)
  end function

  function parallel_inq_dimid(ncid,name,dimid)
    implicit none
    integer :: dimid,ncid,parallel_inq_dimid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_dimid = nf90_inq_dimid(ncid,name,dimid)
    call broadcast(parallel_inq_dimid)
    call broadcast(dimid)
  end function

  function parallel_inq_varid(ncid,name,varid)
    implicit none
    integer :: ncid,parallel_inq_varid,varid
    character(len=*) :: name
    ! begin
    if (main_task) parallel_inq_varid = nf90_inq_varid(ncid,name,varid)
    call broadcast(parallel_inq_varid)
    call broadcast(varid)
  end function

  function parallel_inquire(ncid,nvariables)
    implicit none
    integer :: ncid,parallel_inquire,nvariables
    ! begin
    if (main_task) parallel_inquire = nf90_inquire(ncid,nvariables=nvariables)
    call broadcast(parallel_inquire)
    call broadcast(nvariables)
  end function

  function parallel_inquire_dimension(ncid,dimid,name,len)
    implicit none
    integer :: dimid,ncid,parallel_inquire_dimension
    integer,optional :: len
    character(len=*),optional :: name

    integer :: l

    ! begin

    if (present(name)) then
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,name,len=l)
       call broadcast(name)
    else
       if (main_task) parallel_inquire_dimension = &
            nf90_inquire_dimension(ncid,dimid,len=l)
    end if
    call broadcast(parallel_inquire_dimension)
    if (present(len)) then
       call broadcast(l)
       len = l
    end if
  end function

  function parallel_inquire_variable(ncid,varid,name,ndims,dimids,natts)
    implicit none
    integer :: ncid,parallel_inquire_variable,varid
    integer,optional :: ndims,natts
    character(len=*),optional :: name
    integer,dimension(:),optional :: dimids

    integer :: nd,na
    ! begin
    if (present(name)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,name=name)
       call broadcast(parallel_inquire_variable)
       call broadcast(name)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (present(dimids)) then
       if (main_task) parallel_inquire_variable = &
            nf90_inquire_variable(ncid,varid,dimids=dimids)
       call broadcast(parallel_inquire_variable)
       call broadcast(dimids)
       if (parallel_inquire_variable/=nf90_noerr) return
    end if
    if (main_task) parallel_inquire_variable = &
         nf90_inquire_variable(ncid,varid,ndims=nd,natts=na)
    call broadcast(parallel_inquire_variable)
    if (present(ndims)) then
       call broadcast(nd)
       ndims = nd
    end if
    if (present(natts)) then
       call broadcast(na)
       natts = na
    end if
  end function

  function parallel_open(path,mode,ncid)
    implicit none
    integer :: mode,ncid,parallel_open
    character(len=*) :: path
    ! begin
    if (main_task) parallel_open = nf90_open(path,mode,ncid)
    call broadcast(parallel_open)
  end function

  subroutine parallel_print_real8_2d(name,a)
    implicit none
    character(*) :: name
    real(8),dimension(:,:) :: a
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(a,2)
       do i = 1,size(a,1)
          write(u,*) j,i,a(i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine

  subroutine parallel_print_real8_3d(name,a)
    implicit none
    character(*) :: name
    real(8),dimension(:,:,:) :: a
    
    integer,parameter :: u = 33
    integer :: i,j
    ! begin
    open(unit=u,file=name//".txt",form="formatted",status="replace")
    do j = 1,size(a,3)
       do i = 1,size(a,2)
          write(u,'(2i6,100g15.5e3)') j,i,a(:,i,j)
       end do
       write(u,'()')
    end do
    close(u)
  end subroutine

  function parallel_put_att_character(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_character,varid
    character(len=*) :: name,values
    ! begin
    if (main_task) parallel_put_att_character = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_character)
  end function

  function parallel_put_att_real4(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4,varid
    character(len=*) :: name
    real(4) :: values
    ! begin
    if (main_task) parallel_put_att_real4 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4)
  end function

  function parallel_put_att_real4_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real4_1d,varid
    character(len=*) :: name
    real(4),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real4_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real4_1d)
  end function

  function parallel_put_att_real8(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8,varid
    character(len=*) :: name
    real(8) :: values
    ! begin
    if (main_task) parallel_put_att_real8 = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8)
  end function

  function parallel_put_att_real8_1d(ncid,varid,name,values)
    implicit none
    integer :: ncid,parallel_put_att_real8_1d,varid
    character(len=*) :: name
    real(8),dimension(:) :: values
    ! begin
    if (main_task) parallel_put_att_real8_1d = nf90_put_att(ncid,varid,name,values)
    call broadcast(parallel_put_att_real8_1d)
  end function

  function parallel_put_var_real4(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real4,varid
    integer,dimension(:) :: start
    real(4) :: values
    ! begin
    if (main_task) parallel_put_var_real4 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real4)
  end function

  function parallel_put_var_real8(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8,varid
    integer,dimension(:) :: start
    real(8) :: values
    ! begin
    if (main_task) parallel_put_var_real8 = &
         nf90_put_var(ncid,varid,values,start)
    call broadcast(parallel_put_var_real8)
  end function

  function parallel_put_var_real8_1d(ncid,varid,values,start)
    implicit none
    integer :: ncid,parallel_put_var_real8_1d,varid
    integer,dimension(:),optional :: start
    real(8),dimension(:) :: values
    ! begin
    if (main_task) then
       if (present(start)) then
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values,start)
       else
          parallel_put_var_real8_1d = nf90_put_var(ncid,varid,values)
       end if
    end if
    call broadcast(parallel_put_var_real8_1d)
  end function

  function parallel_redef(ncid)
    implicit none
    integer :: ncid,parallel_redef
    ! begin
    if (main_task) parallel_redef = nf90_redef(ncid)
    call broadcast(parallel_redef)
  end function

  function parallel_reduce_sum(x)
    ! Sum x across all of the nodes.
    ! In parallel_single mode just return x.
    implicit none
    real(8) :: x, parallel_reduce_sum

    parallel_reduce_sum = x
    return
  end function

  function parallel_reduce_max(x)
    ! Max x across all of the nodes.
    ! In parallel_single mode just return x.
    implicit none
    real(8) :: x, parallel_reduce_max

    parallel_reduce_max = x
    return
  end function

  subroutine parallel_show_minmax(label,values)
    implicit none
    character(*) :: label
    real(8),dimension(:,:,:) :: values
    ! begin
    print *,label,minval(values),maxval(values)
  end subroutine

  subroutine parallel_stop(file,line)
    implicit none
    integer :: line
    character(len=*) :: file
    ! begin
    write(0,*) "STOP in ",file," at line ",line
    ! stop
    write(0,*) "RUNNING in parallel_single mode, so STOP IGNORED."
  end subroutine

  function parallel_sync(ncid)
    implicit none
    integer :: ncid,parallel_sync
    ! begin
    if (main_task) parallel_sync = nf90_sync(ncid)
    call broadcast(parallel_sync)
  end function

  subroutine parallel_temp_halo(a)
    implicit none
    real(8),dimension(:,:,:) :: a
  end subroutine

  subroutine parallel_velo_halo(a)
    implicit none
    real(8),dimension(:,:) :: a
  end subroutine

    subroutine staggered_parallel_halo(a)
    implicit none
    real(8),dimension(:,:,:) :: a

    !JEFF Implements a one-sided (staggered) halo update spreading information from the NW to the SE.
    !Remember that the grid is laid out from the NE, then the upper right corner is assigned to task_id = 1.
    !It's western nbhr is task_id = 2, proceeding rowwise and starting from the eastern edge.
    ! What is required for correct calculations is vel(:,ew-1:ew,ns-1:ns), so correct info is needed from NW.
    ! (WAIT, is ns=1 northernmost or southernmost point?  If southernmost, then I need to reverse.)
    ! Communications DO NOT wrap around at the edges.

    if (tasks >= 2) then
	! Do nothing
    endif
  end subroutine

end module parallel
