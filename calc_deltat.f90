module calc_deltat_module

  use multifab_module
  use ml_layout_module
  use params_module

  implicit none

  private

  public :: calc_deltat 

contains

  !if we do not do time refinement we just need one global deltaT over
  !all levels.  If we do time refinement we will need a calc_deltat that
  !operates on a single level

  subroutine calc_deltat(mla,state,dx,deltaT) 

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: state(:)
    real(kind=dp_t), intent(in   ) :: dx(mla%nlevel,mla%dim)
    real(kind=dp_t), intent(out  ) :: deltaT

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: nlevs, dm, ng, i, n, nc

    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t) :: deltaT_fab,deltaT_loc

    ng = state(1)%ng    !number of ghosts
    dm = mla%dim        !dimension
    nlevs = mla%nlevel  !number of levels
    nc = state(1)%nc    !number of components

    !initialize deltaT to user-defined maximum value  
    deltaT = max_deltaT

    do n=1,nlevs

       do i=1,nfabs(state(n))
          dp => dataptr(state(n),i)
          lo = lwb(get_box(state(n),i))
          hi = upb(get_box(state(n),i))
          select case(dm)
          case (2)
             call calc_deltat_2d(dp(:,:,1,:), nc, ng, lo, hi, dx(n,:),deltaT_fab)
             deltaT = min(deltaT,deltaT_fab)
          case (3)
            call bl_error('need to setup calc_deltat_3d')
          end select
       end do

    end do

    !all reduce the min
    deltaT_loc = CFL*deltaT
    call parallel_reduce(deltaT, deltaT_loc, MPI_MIN)
    


  end subroutine calc_deltat 

  subroutine calc_deltat_2d(state, nc, ng, lo, hi, dx, deltaT_fab)

    integer    :: lo(2), hi(2), nc, ng, fabid
    real(dp_t) :: state(lo(1)-ng:hi(1)+ng,lo(2)-ng:hi(2)+ng,nc)
    real(dp_t) :: dx(2)
    real(dp_t), intent(out) :: deltaT_fab

    ! local varables
    real(dp_t) :: maxu,maxv

    !try using fortran intrinsic on array slice?
    maxu = maxval(abs(state(lo(1):hi(1),lo(2):hi(2),uvel_id)))
    maxv = maxval(abs(state(lo(1):hi(1),lo(2):hi(2),vvel_id)))

    deltaT_fab = min(dx(1)/(maxu+dp_tiny),dx(2)/(maxv+dp_tiny))


  end subroutine calc_deltat_2d


  !   subroutine init_state_3d(state, ng, lo, hi, dom_lo, dx)
  !
  !   integer          :: lo(3), hi(3), ng
  !   real(dp_t) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
  !   real(dp_t) :: dom_lo(3)
  !   real(dp_t) :: dx
  !
  !   ! local varables
  !   integer          :: i,j,k
  !   real(dp_t) :: x,y,z,r1,r2,r3
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

end module calc_deltat_module
