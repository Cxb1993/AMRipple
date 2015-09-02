program AMRipple  

  use boxlib
  use bl_types
  use parallel
  use make_new_grids_module
  use regrid_module
  use define_bc_module
  use bc_module
  use bndry_reg_module
  use ml_layout_module
  use utils_mod, only : pwrite
  use init_state_module
  use calc_deltat_module
  use check_state_module
  use params_module
  use write_plotfile_module

  implicit none


  integer, parameter :: str_length = 80
  integer, parameter :: stdout = 6
  integer    :: un
  logical    :: fexist
  logical    :: haltNexit
  logical    :: new_grid
  integer    :: n,ierr,nl
  integer    :: nc,nlevs
  real(dp_t) :: start_time,run_time,run_time_IOproc
  real(dp_t) :: cluster_min_eff
  integer    :: max_levs,ndims,nsteps,plot_int
  integer    :: n_celli,n_cellj,n_cellk
  integer    :: max_grid_size,amr_buf_width
  integer    :: cluster_minwidth,cluster_blocking_factor,regrid_int
  integer    :: istep
  character(str_length) :: pws
  real(dp_t) :: T,deltaT

  ! size ndims
  integer, allocatable :: lo(:),hi(:)
  logical, allocatable :: is_periodic(:)
  real(dp_t), allocatable :: dom_lo(:),dom_hi(:)
  ! real(dp_t), allocatable :: xlo,xhi,ylo,yhi,zlo,zhi
  ! real(dp_t), allocatable :: bc_xlo,bc_xhi,bc_ylo,bc_yhi,bc_zlo,bc_zhi
  real(dp_t) :: xlo,xhi,ylo,yhi,zlo,zhi
  real(dp_t) :: bc_xlo,bc_xhi,bc_ylo,bc_yhi,bc_zlo,bc_zhi

  ! size (ndims,2)
  integer       , allocatable :: phys_bc(:,:)

  ! size (max_levs,ndims)
  real(dp_t), allocatable :: dx(:,:)


  type(box)         :: bx
  type(ml_boxarray) :: mbla
  type(ml_layout)   :: mla
  type(layout),    allocatable :: la_array(:)
  type(bndry_reg), allocatable :: bndry_flx(:)
  type(bc_tower)    :: the_bc_tower
  type(multifab),  allocatable :: state(:)   !holds current state
  type(multifab),  allocatable :: state1(:)  !holds state at n-1

  namelist /probin/ max_levs , ndims, nsteps, plot_int, n_celli, n_cellj, &
  n_cellk, max_grid_size, amr_buf_width, cluster_minwidth, &
  cluster_blocking_factor, cluster_min_eff, regrid_int, xlo, xhi, &
  ylo, yhi, zlo, zhi, dens_ratio, max_deltaT,CFL

  call boxlib_initialize()

  ! start the clock
  start_time = parallel_wtime()

  write(*,*)'params set'

  inquire(file='input.dat',exist=fexist)
  if(.not.fexist)then
    call bl_error('input file input.dat does not exist')
  endif
  un = unit_new()
  open(unit=un, file = 'input.dat', status = 'old', action = 'read')
  read(unit=un, nml = probin)
  close(unit=un)


  bc_xlo = PERIODIC !SLIP_WALL
  bc_xhi = PERIODIC !SLIP_WALL
  bc_ylo = SLIP_WALL
  bc_yhi = SLIP_WALL
  bc_zlo = SLIP_WALL
  bc_zhi = SLIP_WALL

  ! 2D:  f, rho, p, U, V
  ! 3D:  f, rho, p, U, V, W
  nc = 3 + ndims !number of state variables

  ! set the indices for the state variables
  if(ndims==2)then
    f_id    = 1
    rho_id  = 2
    pres_id = 3
    uvel_id = 4
    vvel_id = 5
    wvel_id = -1 !hopefully will set off runtime error
  elseif(ndims==3)then
    f_id    = 1
    rho_id  = 2
    pres_id = 3
    uvel_id = 4
    vvel_id = 5
    wvel_id = 6 
  else
    call bl_error('WTF?, ndims is not 2 or 3')
  endif

  ! sanity check on input params
  haltNexit = .false.
  if(max_grid_size > max(n_celli,n_cellj,n_cellk))then
    call pwrite('max grid size exceeds max of n_celli,j,k')
    haltNexit = .true.
  endif
  if(mod(n_celli,max_grid_size) /= 0)then
    call pwrite('n_celli must be a multiple of max_grid_size')
    haltNexit = .true.
  endif
  if(mod(n_cellj,max_grid_size) /= 0)then
    call pwrite('n_cellj must be a multiple of max_grid_size')
    haltNexit = .true.
  endif
  if(ndims==3)then
    if(mod(n_cellk,max_grid_size) /= 0)then
      call pwrite('n_cellk must be a multiple of max_grid_size')
      haltNexit = .true.
    endif
  endif
  if(haltNexit)then
    call pwrite('exiting....')
    call MPI_Finalize(ierr); stop
  endif


  ! allocate variables requiring ndims
  allocate(lo(ndims),hi(ndims))
  allocate(dom_lo(ndims),dom_hi(ndims))
  allocate(is_periodic(ndims))
  allocate(phys_bc(ndims,2))

  ! variables that depend on max_levs
  allocate(dx(max_levs,ndims)); dx = 0.0
  allocate(la_array(max_levs))
  allocate(state(max_levs))
  allocate(state1(max_levs))

  ! setup refinement and clustering parameters
  call cluster_set_minwidth(cluster_minwidth)
  call cluster_set_blocking_factor(cluster_blocking_factor)
  call cluster_set_min_eff(cluster_min_eff)
  call amr_ref_ratio_init(max_levs,ndims,2)

  ! tell mbla about max_levs and dimensionality of problem
  call ml_boxarray_build_n(mbla,max_levs,ndims)
  ! tell mbla about the ref_ratio between levels
  call ml_boxarray_set_ref_ratio(mbla)

  ! set physical problem domain dimsensions
  if(ndims==2)then
    dom_lo(1:2) = [xlo,ylo]
    dom_hi(1:2) = [xhi,yhi]
  else
    dom_lo(1:3) = [xlo,ylo,zlo]
    dom_hi(1:3) = [xhi,yhi,zhi]
  endif

  ! set grid spacing at each level for each direction => dx(max_levs,ndims)
  ! must set amr_ref_ratio_init and ml_boxarray_set_ref_ratio first
  do n=1,max_levs
    if(n==1)then
      dx(n,1) = (dom_hi(1)-dom_lo(1)) / n_celli
      dx(n,2) = (dom_hi(2)-dom_lo(2)) / n_cellj
      if(ndims==3)then
        dx(n,3) = (dom_hi(3)-dom_lo(3)) / n_cellk
      endif
    else
      dx(n,1) = dx(n-1,1) / mbla%rr(n-1,1)
      dx(n,2) = dx(n-1,2) / mbla%rr(n-1,2)
      if(ndims==3)then
        dx(n,3) = dx(n-1,3) / mbla%rr(n-1,3)
      endif
    endif
    ! if (parallel_IOProcessor())then
    !   write(*,*)
    !   if(ndims==3)then
    !     write(*,*)'level: ',n,' resolution: ',dx(n,1),dx(n,2),dx(n,3)
    !   else
    !     write(*,*)'level: ',n,' resolution: ',dx(n,1),dx(n,2)
    !   endif
    ! endif
  enddo

  !-----------------------------------------------------------------------------
  ! setup the boundary conditions
  !-----------------------------------------------------------------------------
  phys_bc(1,1) = bc_xlo
  phys_bc(1,2) = bc_xhi
  phys_bc(2,1) = bc_ylo
  phys_bc(2,2) = bc_yhi
  if (ndims == 3) then
    phys_bc(3,1) = bc_zlo
    phys_bc(3,2) = bc_zhi
  end if

  ! set periodicity
  ! build an array indicating periodicity in each direction
  is_periodic(:) = .false.
  do n=1,ndims
    if (phys_bc(n,1) == -1 .and. phys_bc(n,2) /= -1) then
      call bl_error("Invalid BC's - both lo and hi need to be periodic")
    end if
    if (phys_bc(n,2) == -1 .and. phys_bc(n,1) /= -1) then
      call bl_error("Invalid BC's - both lo and hi need to be periodic")
    end if
    if (phys_bc(n,1) == -1 .and. phys_bc(n,2) == -1) then
      is_periodic(n) = .true.
    end if
  end do

  ! put all the domain boundary conditions into phys_bc
  call bc_tower_init(the_bc_tower,max_levs,ndims,phys_bc)

  ! create a box (type(box) :: bx) from (0,0) to (n_cell-1,n_cell-1)
  if(ndims==3)then
    bx = make_box([0,0,0],[n_celli-1,n_cellj-1,n_cellk-1])
  else
    bx = make_box([0,0],[n_celli-1,n_cellj-1])
  endif

  ! tell multi-level box array mbla about the problem domain at each level
  mbla%pd(1) = bx
  do n=2,max_levs
    mbla%pd(n) = refine(mbla%pd(n-1),mbla%rr((n-1),:))
  enddo

  ! initialize the boxarray at level 1 to be one single box
  call boxarray_build_bx(mbla%bas(1),bx)

  ! overwrite the boxarray at level 1 to respect max_grid_size
  call boxarray_maxsize(mbla%bas(1),max_grid_size)

  ! build the level 1 layout
  call layout_build_ba(la_array(1),mbla%bas(1),mbla%pd(1),is_periodic)

  ! build the level 1 multifab with nc and components two ghost cells
  call multifab_build(state(1),la_array(1),nc,2)

  ! define level 1 of the_bc_tower
  call bc_tower_level_build(the_bc_tower,1,la_array(1),nc)

  ! todo, we need to work on bcs, the routines in the boxlib are not general.


  !-----------------------------------------------------------------------------
  ! initialize the solution on the AMR mesh
  !-----------------------------------------------------------------------------

  ! initialize the state variables on the base (coarse) level grid
  call init_state_on_level(state(1),dx(1,:),dom_lo,the_bc_tower%bc_tower_array(1))
  !call check_state_on_level(state(1))

  nl = 1
  new_grid = .true.

  do while ( (nl .lt. max_levs) .and. (new_grid) )

    ! determine whether we need finer grids based on tagging criteria
    ! if so, return new_grid=T and the la_array(nl+1)
    ! here, dx is a scalar (where we would like the option of it being dependent
    ! on direction).  However, it is not ever used in the call.
    ! make_new_grids has an optional argument that allows you to include a multifab
    ! for tagging cells.  However, if you follow this argument [aux_tag_mf] from
    ! make_new_grids to make_boxes to tag_boxes, you see that it errors out with
    ! aux_tag_mf passed to tag_boxes without implementation.  This is indication to
    ! me that we should perhaps be either modifying the source or copying certain routines
    ! like tag_boxes and having a local copy.  This also makes sense because tag_boxes
    ! is a bit arbitrary.  It uses the first field of the multfab passed in to tag boxes
    ! and those are tagged based on whether that value is greater than some arbitrary
    ! number.  I think we need a code-specific of tag_boxes.

    call make_new_grids(new_grid,la_array(nl),la_array(nl+1),state(nl),dx(nl,1), &
    amr_buf_width,mbla%rr(nl,1),nl,max_grid_size)


    if (new_grid) then

      ! tell mbla about the finer level boxarray
      call copy(mbla%bas(nl+1),get_boxarray(la_array(nl+1)))

      ! Build the level nl+1 data
      call multifab_build(state(nl+1),la_array(nl+1),state(nl)%nc,2)

      ! define level nl+1 of the_bc_tower
      call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1),nc)

      ! initialize phi on level nl+1
      write(*,*)'initing the state on level: ',nl+1
      call init_state_on_level(state(nl+1),dx(nl+1,:),dom_lo,the_bc_tower%bc_tower_array(nl+1))

      ! increment current level counter
      nl = nl+1

    endif

  end do
  !call bl_error('done initial refine')
  ! the current number of levels in the simulation is nlevs, not necessarily max_levs
  nlevs = nl

  ! destroy the state vector - we are going to build it again using the new multilevel
  ! layout after we have tested and reconfigured the grids due to proper nesting
  do n=1,nlevs
    call multifab_destroy(state(n))
  end do

  if (nlevs .ge. 3) then
    ! check for proper nesting
    call enforce_proper_nesting(mbla,la_array,max_grid_size)
  end if

  do n = 1,nlevs
    call destroy(la_array(n))
  end do

  call ml_layout_restricted_build(mla,mbla,nlevs,is_periodic)

  call destroy(mbla)

  ! this makes sure the boundary conditions are properly defined everywhere
  do n = 1,nlevs
    call bc_tower_level_build(the_bc_tower,n,mla%la(n),nc)
  end do

  ! We will use two ghost cells for now
  do n=1,nlevs
    call multifab_build(state(n),mla%la(n),nc,2)   !current
    call multifab_build(state1(n),mla%la(n),nc,2)  !last time step
  end do

  call init_state(mla,state,dx,dom_lo,the_bc_tower)


  !=============================================================================
  ! dump an initial condition file for visualization of AMR / ICs
  !=============================================================================
  call write_plotfile(mla,state,0,dx(1,:),ZERO,dom_lo,dom_hi)


  !=============================================================================
  ! main loop
  !=============================================================================
  T = 0
  do istep=1,nsteps

    !----------------------------------------------------------
    ! compute a time step using convection only 
    ! (u,v,w) is zero on the first iteration so we will have to
    ! threshold dt to maxdt set in input.dat
    !----------------------------------------------------------
    call calc_deltat(mla,state,dx,deltaT) 
    if (parallel_IOProcessor())then
      write(*,*)'time step : ',deltaT 
    endif
    call bl_error('geoff stopping here')

    !----------------------------------------------------------
    ! advection/diffusion/soure term
    ! right now we are including only the source term
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! compute divergence for RHS 
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! solve Poisson for new pressure field 
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! pressure correction on velocity 
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! vof equation 
    !----------------------------------------------------------
    !call advance_vof( )

    !----------------------------------------------------------
    ! update density 
    !----------------------------------------------------------

    !----------------------------------------------------------
    ! dump output 
    !----------------------------------------------------------
    if (plot_int .gt. 0) then
      if (mod(istep,plot_int) .eq. 0 .or. istep .eq. nsteps) then
        call write_plotfile(mla,state,istep,dx(1,:),deltaT,dom_lo,dom_hi)
      end if
    end if

    !----------------------------------------------------------
    ! Report progress
    !----------------------------------------------------------
    if (parallel_IOProcessor())then
      write(*,*)'completed : ',istep,' of ',nsteps
    endif

    call sleep(1)
   
    ! increment time
    T = T + deltaT

  end do

  !=============================================================================
  ! end main loop
  !=============================================================================


  ! check for memory that should have been deallocated
  if ( parallel_IOProcessor() ) then
    write(*,*) 'MEMORY STATS AT END OF PROGRAM'
  end if
  call print(multifab_mem_stats(),    "    multifab")
  call print(fab_mem_stats(),         "         fab")
  call print(boxarray_mem_stats(),    "    boxarray")
  call print(layout_mem_stats(),      "      layout")
  call print(boxassoc_mem_stats(),    "    boxassoc")
  call print(fgassoc_mem_stats(),     "     fgassoc")
  call print(syncassoc_mem_stats(),   "   syncassoc")
  call print(copyassoc_mem_stats(),   "   copyassoc")
  call print(fluxassoc_mem_stats(),   "   fluxassoc")


  ! end the clock
  run_time = parallel_wtime() - start_time

  ! collect run_time from each processor and store the maximum
  call parallel_reduce(run_time_IOproc, run_time, MPI_MAX, &
  proc = parallel_IOProcessorNode())

  ! report total program run time
  if ( parallel_IOProcessor() ) then
    write(*,'(A20,F10.2)'),"Total Run time (s) =",run_time_IOproc
  end if

  !finalize the code and exit
  call boxlib_finalize()

  

end program AMRipple  
