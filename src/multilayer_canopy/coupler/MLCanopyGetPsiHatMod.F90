module MLCanopyGetPsiHatMod

  public LookupPsihatINI

  contains
  
  !-----------------------------------------------------------------------
  subroutine LookupPsihatINI(rslfile,masterproc)
    !
    ! !DESCRIPTION:
    ! Initialize the look-up tables needed to calculate the RSL psihat functions.
    ! Remember that in a netcdf file the dimensions appear in the opposite order:
    ! netcdf: psigridM_nc(nL,nZ) -> Fortran: psigridM(nZ,nL)
    ! netcdf: psigridH_nc(nL,nZ) -> Fortran: psigridH(nZ,nL)
    !
    ! !USES:
    use fileutils, only : getfil
    use ncdio_pio, only : ncd_io, ncd_pio_closefile, ncd_pio_openfile, file_desc_t
    use ncdio_pio, only : ncd_inqdid, ncd_inqdlen
    use MLCanopyVarCon, only : nZ, nL, dtLgridM, zdtgridM, psigridM, dtLgridH, zdtgridH, psigridH
    !
    implicit none

    ! !ARGUMENTS:
    character(len=256), intent(in) :: rslfile
    logical,intent(in) :: masterproc
    
    !
    !LOCAL VARIABLES
    character(len=256) :: locfn        ! Local file name
    type(file_desc_t) :: ncid          ! pio netCDF file id
    integer :: dimid                   ! netCDF dimension id
    logical :: readv                   ! read variable in or not

    real(r8) :: zdtgridM_nc(nZ)        ! netcdf data: Grid of zdt on which psihat is given for momentum
    real(r8) :: dtLgridM_nc(nL)        ! netcdf data: Grid of dtL on which psihat is given for momentum
    real(r8) :: psigridM_nc(nL,nZ)     ! netcdf data: Grid of psihat values for momentum
    real(r8) :: zdtgridH_nc(nZ)        ! netcdf data: Grid of zdt on which psihat is given for heat
    real(r8) :: dtLgridH_nc(nL)        ! netcdf data: Grid of dtL on which psihat is given for heat
    real(r8) :: psigridH_nc(nL,nZ)     ! netcdf data: Grid of psihat values for heat
    integer :: nZ_nc, nL_nc            ! netcdf data: dimensions
    integer :: ii, jj                  ! Looping indices
    !---------------------------------------------------------------------

    if (masterproc) then
       write(iulog,*) 'Attempting to read RSL look-up table .....'
    end if

    ! Get netcdf file

    call getfil (rslfile, locfn, 0)

    ! Open netcdf file

    call ncd_pio_openfile (ncid, trim(locfn), 0)

    ! Check dimensions

    call ncd_inqdid (ncid, 'nZ', dimid)
    call ncd_inqdlen (ncid, dimid, nZ_nc)

    if (nZ_nc /= nZ) then
       call endrun (msg=' ERROR: LookupPsihatINI: nZ does not equal expected value')
    end if

    call ncd_inqdid (ncid, 'nL', dimid)
    call ncd_inqdlen (ncid, dimid, nL_nc)

    if (nL_nc /= nL) then
       call endrun (msg=' ERROR: LookupPsihatINI: nL does not equal expected value')
    end if

    ! Read variables

    call ncd_io('dtLgridM', dtLgridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading dtLgridM')

    call ncd_io('zdtgridM', zdtgridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading zdtgridM')

    call ncd_io('psigridM', psigridM_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading psigridM')

    call ncd_io('dtLgridH', dtLgridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading dtLgridH')

    call ncd_io('zdtgridH', zdtgridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading zdtgridH')

    call ncd_io('psigridH', psigridH_nc, 'read', ncid, readvar=readv, posNOTonfile=.true.)
    if (.not. readv) call endrun (msg=' ERROR: LookupPsihatINI: error reading psigridH')

    ! Close netcdf file

    call ncd_pio_closefile(ncid)

    if (masterproc) then
       write(iulog,*) 'Successfully read RSL look-up table'
    end if

    ! Copy netcdf variables

    do jj = 1, nL
       dtLgridM(1,jj) = dtLgridM_nc(jj)
       dtLgridH(1,jj) = dtLgridH_nc(jj)
    end do

    do ii = 1, nZ
       zdtgridM(ii,1) = zdtgridM_nc(ii)
       zdtgridH(ii,1) = zdtgridH_nc(ii)
    end do

    do ii = 1, nZ
       do jj = 1, nL
          psigridM(ii,jj) = psigridM_nc(jj,ii)
          psigridH(ii,jj) = psigridH_nc(jj,ii)
       end do
    end do

    return

  end subroutine LookupPsihatINI
  
end module MLCanopyGetPsiHatMod
