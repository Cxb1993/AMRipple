module check_state_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use ml_restrict_fill_module

  implicit none

  private

  public :: check_state_on_level

contains

  subroutine check_state_on_level(state)

    type(multifab) , intent(in) :: state


    ! local
    integer i,ng,dm,ns
    integer :: lo(state%dim), hi(state%dim)

    real(kind=dp_t), pointer :: dp(:,:,:,:)

    ng = state%ng  !number of ghost cells
    dm = state%dim !problem dimension (2d/3d)
    ns = state%nc  !number of components

    do i=1,nfabs(state)
      dp => dataptr(state,i)
      lo = lwb(get_box(state,i))
      hi = upb(get_box(state,i))
      select case(dm)
      case (2)
        call check_state_2d(dp(:,:,1,:), ns, ng, lo, hi, i)
      case (3)
        write(*,*)'haven not set 3d check yet'
        stop
        !call init_state_3d(dp(:,:,:,:), ng, lo, hi, dom_lo, dx)
      end select
    end do


  end subroutine check_state_on_level

  ! subroutine init_state(mla,state,dx,dom_lo,the_bc_tower)
  !
  !   type(ml_layout), intent(in   ) :: mla
  !   type(multifab) , intent(inout) :: state(:)
  !   real(kind=dp_t), intent(in   ) :: dx(:)
  !   real(kind=dp_t), intent(in   ) :: dom_lo(mla%dim)
  !   type(bc_tower) , intent(in   ) :: the_bc_tower
  !
  !   ! local variables
  !   integer :: lo(mla%dim), hi(mla%dim)
  !   integer :: nlevs, dm, ng, i, n
  !
  !   real(kind=dp_t), pointer :: dp(:,:,:,:)
  !
  !   ng = state(1)%ng
  !   dm = mla%dim
  !   nlevs = mla%nlevel
  !
  !   do n=1,nlevs
  !
  !      do i=1,nfabs(state(n))
  !         dp => dataptr(state(n),i)
  !         lo = lwb(get_box(state(n),i))
  !         hi = upb(get_box(state(n),i))
  !         select case(dm)
  !         case (2)
  !            call init_state_2d(dp(:,:,1,1), ng, lo, hi, dom_lo, dx(n))
  !         case (3)
  !            call init_state_3d(dp(:,:,:,1), ng, lo, hi, dom_lo, dx(n))
  !         end select
  !      end do
  !
  !   end do
  !
  !   ! restrict the multi-level data, and
  !   ! fill all boundaries: same-level, coarse-fine, periodic, and domain boundaries
  !   call ml_restrict_and_fill(nlevs, state, mla%mba%rr, the_bc_tower%bc_tower_array)
  !
  ! end subroutine init_state

  subroutine check_state_2d(state, ns, ng, lo, hi, fabid)

    integer          :: lo(2), hi(2), ns, ng, fabid
    double precision :: state(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,ns)

    ! local varables
    integer          :: i,j,n,ngl

    ngl = 0!ng
    do n=1,ns
      do j=lo(2)-ngl,hi(2)+ngl
        !y = dom_lo(2) + (dble(j)+0.5d0) * dx(2)
        do i=lo(1)-ngl,hi(1)+ngl
          !x = dom_lo(1) + (dble(i)+0.5d0) * dx(1)
          write(*,*) 'check out', fabid,n,i,j,state(i,j,n)
        end do
      end do
    end do

  end subroutine check_state_2d

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

end module check_state_module
