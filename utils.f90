module utils_mod


  use parallel


  interface pwrite
    module procedure pwrite_1
    module procedure pwrite_2
  end interface pwrite


contains
  subroutine pwrite_1(message)
    character(len= *), intent(in) :: message

    if(parallel_IOProcessor())then
      write(*,*) message
    endif

  end subroutine pwrite_1

  subroutine pwrite_2(iout,message)
    integer, intent(in) :: iout
    character(len= *), intent(in) :: message

    if(parallel_IOProcessor())then
      write(iout,*) message
    endif

  end subroutine pwrite_2

end module utils_mod
