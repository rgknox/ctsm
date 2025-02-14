module MLCanopyVarPar

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing multilayer canopy model parameters
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none
  save

  public :: MLCanopySetVarPar
  
  ! Parameters for multilayer canopy

  integer, parameter :: nlevmlcan = 100     ! Number of layers in multilayer canopy model
  integer, parameter :: nleaf = 2           ! Number of leaf types (sunlit and shaded)
  integer, parameter :: isun = 1            ! Sunlit leaf index
  integer, parameter :: isha = 2            ! Shaded leaf index

  ! Parameter constants passed in from the host
  real(r8) :: spval = 1.e36_r8              ! special value for real data
  integer  :: ispval = -9999                ! special value for int data
  integer  :: nlevgrnd                      ! Number of ground layers
  integer  :: numrad                        ! Number of SW radiation bands (vis/nir)
  integer  :: ivis                          ! Index of visible radiation band
  integer  :: inir                          ! Index of near-infrared radiation band
  integer  :: iulog                         ! File unit number for text logging

contains
  
  subroutine MLCanopySetVarPar(varname,rval,ival,cval)
    
    integer, optional, intent(in)         :: ival
    real(r8), optional, intent(in)        :: rval
    character(len=*),optional, intent(in) :: cval
    character(len=*),intent(in)           :: varname
    
    select case(trim(varname))
    case('spval') spval       = rval
    case('ispval') ispval     = ival
    case('nlevgrnd') nlevgrnd = ival
    case('numrad') numrad     = ival
    case('ivis') ivis         = ival
    case('inir') inir         = ival
    case('iulog') iulog       = ival
    end select
       
  end subroutine MLCanopySetVarPar
  
end module MLCanopyVarPar
