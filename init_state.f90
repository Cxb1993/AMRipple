module init_state_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use ml_restrict_fill_module
  use params_module

  implicit none

  private

  public :: init_state_on_level, init_state

contains

  subroutine init_state_on_level(state,dx,dom_lo,the_bc_level)

    type(multifab) , intent(inout) :: state
    real(kind=dp_t), intent(in   ) :: dx(state%dim)
    real(kind=dp_t), intent(in   ) :: dom_lo(state%dim)
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer i,ng,dm,nc
    integer :: lo(state%dim), hi(state%dim)

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ng = state%ng  !number of ghost cells
    dm = state%dim !problem dimension (2d/3d)
    nc = state%nc  !number of components

    do i=1,nfabs(state)
       dp => dataptr(state,i)
       lo = lwb(get_box(state,i))
       hi = upb(get_box(state,i))
       select case(dm)
       case (2)
          call init_state_2d(dp(:,:,1,:), nc, ng, lo, hi, dom_lo, dx, i)
       case (3)
         call bl_error('need to setup init_state_3d')
          !call init_state_3d(dp(:,:,:,:), ng, lo, hi, dom_lo, dx)
       end select
    end do

    ! I assume this sets the ghost cells in the interior?
    call multifab_fill_boundary(state)

    !arguments appear to be
    !multfab, first_state, first_bc_state?, nstate, bc_info
    !the loop inside over the data runs from first_state to first_state+nstate-1
    !note, if that fourth argument, which I thought was nstate is greater than 2,
    !the code crashes.  Are we doing something wrong setting up the_bc_level?
    call multifab_physbc(state,1,1,nc,the_bc_level)

  end subroutine init_state_on_level

  subroutine init_state(mla,state,dx,dom_lo,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: state(:)
    real(kind=dp_t), intent(in   ) :: dx(mla%nlevel,mla%dim)
    real(kind=dp_t), intent(in   ) :: dom_lo(mla%dim)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n, nc

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ng = state(1)%ng
    dm = mla%dim
    nlevs = mla%nlevel
    nc = state(1)%nc  !number of components

    do n=1,nlevs

       do i=1,nfabs(state(n))
          dp => dataptr(state(n),i)
          lo = lwb(get_box(state(n),i))
          hi = upb(get_box(state(n),i))
          select case(dm)
          case (2)
             call init_state_2d(dp(:,:,1,:), nc, ng, lo, hi, dom_lo, dx(n,:), i)
          case (3)
            call bl_error('need to setup init_state_3d')
             !call init_state_3d(dp(:,:,:,1), ng, lo, hi, dom_lo, dx(n,:))
          end select
       end do

    end do

    ! restrict the multi-level data, and
    ! fill all boundaries: same-level, coarse-fine, periodic, and domain boundaries
    call ml_restrict_and_fill(nlevs, state, mla%mba%rr, the_bc_tower%bc_tower_array)

  end subroutine init_state

  subroutine init_state_2d(state, ns, ng, lo, hi, dom_lo, dx, fabid)

    integer          :: lo(2), hi(2), ns, ng, fabid
    double precision :: state(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,ns)
    !double precision :: state(lo(1)-ng:,lo(2)-ng:,:)
    double precision :: dom_lo(2)
    double precision :: dx(2)

    ! local varables
    integer          :: i,j,n
    double precision :: x,y,ygauss

    state = 0.0  !U,V,W,P
    !do j=lo(2),hi(2)
    do j=lo(2)-ng,hi(2)+ng
       y = dom_lo(2) + (dble(j)+0.5d0) * dx(2)
       !do i=lo(1),hi(1)
       do i=lo(1)-ng,hi(1)+ng
          x = dom_lo(1) + (dble(i)+0.5d0) * dx(1)
          ygauss = .4 + .2*exp(-(x-0.5)*(x-0.5)/(2*.2*.2))
          if(y < ygauss)then
            state(i,j,f_id) = 1.0  !volume fraction f
          endif
          state(i,j,rho_id) = state(i,j,1)*dens_ratio + (1-state(i,j,1)) !rho
          state(i,j,pres_id) = 0.0 !P
          state(i,j,uvel_id) = 1.0 !0.0 !U
          state(i,j,vvel_id) = 0.0 !V
       end do
    end do


    end subroutine init_state_2d


  !   subroutine init_state_3d(state, ng, lo, hi, dom_lo, dx)
  !
  !   integer          :: lo(3), hi(3), ng
  !   double precision :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  !   double precision :: dom_lo(3)
  !   double precision :: dx
  !
  !   ! local varables
  !   integer          :: i,j,k
  !   double precision :: x,y,z,r1,r2,r3
  !
  !   !$omp parallel do private(i,j,k,x,y,z,r2)
  !   do k=lo(3),hi(3)
  !      z = dom_lo(3) + (dble(k)+0.5d0) * dx
  !      do j=lo(2),hi(2)
  !         y = dom_lo(2) + (dble(j)+0.5d0) * dx
  !         do i=lo(1),hi(1)
  !            x = dom_lo(1) + (dble(i)+0.5d0) * dx
  !
  !            r1 = ((x-0.5d0)**2 + (y-0.5d0)**2 + (z-0.0d0)**2) / 0.01d0
  !            r2 = ((x-0.0d0)**2 + (y-0.0d0)**2 + (z-0.0d0)**2) / 0.01d0
  !            r3 = ((x+0.5d0)**2 + (y+0.5d0)**2 + (z-0.0d0)**2) / 0.01d0
  !
  !            state(i,j,k) = 1.d0 + exp(-r1) + exp(-r2) + exp(-r3)
  !
  !         end do
  !      end do
  !   end do
  !   !$omp end parallel do
  !
  ! end subroutine init_state_3d

end module init_state_module
