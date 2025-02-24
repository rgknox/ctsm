module shr_sys_mod

  ! Majorly stripped down version of shr_sys_mod to
  ! provide a shutdown method for unit tests
  
  public :: shr_sys_abort

contains

  subroutine shr_sys_abort
    call exit(0)
  end subroutine shr_sys_abort
  
end module shr_sys_mod
