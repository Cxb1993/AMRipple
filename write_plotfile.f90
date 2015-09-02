module write_plotfile_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use params_module

  implicit none

contains

  subroutine write_plotfile(mla,state,istep,dx,time,prob_lo,prob_hi)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: state(:)
    integer        , intent(in   ) :: istep
    real(dp_t)     , intent(in   ) :: dx(:),time
    real(dp_t)     , intent(in   ) :: prob_lo(mla%dim), prob_hi(mla%dim)

    ! local variables
    character(len=8)  :: plotfile_name
    character(len=20) :: variable_names(state(1)%nc)

    type(multifab), allocatable :: plotdata(:)

    ! dimensioned as an array with size dm for fabio_ml_multifab_write_d
    real(dp_t) :: dx_vec(mla%dim)

    integer n, nlevs, nc

    nlevs = mla%nlevel
    nc    = state(1)%nc

    dx_vec = dx

    variable_names(f_id) = "f"
    variable_names(rho_id) = "rho"
    variable_names(pres_id) = "p"
    variable_names(uvel_id) = "U"
    variable_names(vvel_id) = "V"
    if(mla%dim == 3)then
      variable_names(wvel_id) = "W"
    endif


    allocate(plotdata(nlevs))

    do n=1,nlevs
      ! build plotdata with 1 component and 0 ghost cells
      call multifab_build(plotdata(n),mla%la(n),nc,0)

      ! copy the state into plotdata
      call multifab_copy_c(plotdata(n),1,state(n),1,nc)
    end do

    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='("plt",i5.5)') istep

    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), plotfile_name, variable_names, &
    mla%mba%pd(1), prob_lo, prob_hi, time, dx_vec)

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
      call multifab_destroy(plotdata(n))
    end do

    !sync up the processors
    call parallel_barrier()

  end subroutine write_plotfile

end module write_plotfile_module
