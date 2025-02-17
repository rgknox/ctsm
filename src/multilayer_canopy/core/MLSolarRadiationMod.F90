module MLSolarRadiationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate solar radiation transfer through canopy
  !
  ! Important Note*  This module assumes a canopy scattering environment
  !                  is the same PFT for all elements in that patch
  !
  !
  ! !USES:
  use MLCanopyVarCtl, only : endrun
  use clm_varctl, only : iulog
  use shr_kind_mod, only : r8 => shr_kind_r8
  use MLCanopyFluxesType, only : mlcanopy_type
  !
  ! !PUBLIC TYPES:
  implicit none
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SolarRadiation        ! Main driver for radiative transfer
  public :: AllocateMLSolarParams
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: Norman               ! Norman radiative transfer
  private :: TwoStream            ! Two-stream approximation

  !-----------------------------------------------------------------------

  type, public :: rad_params_type
     
     ! From the parameter file
     real(r8), allocatable :: rhol(:,:)         ! leaf material reflectance:   (pft x band)
     real(r8), allocatable :: rhos(:,:)         ! stem material reflectance:   (pft x band)
     real(r8), allocatable :: taul(:,:)         ! leaf material transmittance: (pft x band)
     real(r8), allocatable :: taus(:,:)         ! stem material transmittance: (pft x band)
     real(r8), allocatable :: xl(:)             ! leaf/stem orientation (pft)
     real(r8), allocatable :: clump_fac(:)      ! clumping index 0-1, when
                                                ! leaves stick together (pft)
  end type rad_params_type

  type(rad_params_type),public :: rad_params
  
contains

  !-----------------------------------------------------------------------

  subroutine AllocateMLSolarParams(n_pft)
      
    integer,intent(in) :: n_pft
    
    ! Include the zeroth pft index for air
    
    allocate(rad_params%rhol(0:n_pft,numrad))
    allocate(rad_params%rhos(0:n_pft,numrad))
    allocate(rad_params%taul(0:n_pft,numrad))
    allocate(rad_params%taus(0:n_pft,numrad))
    allocate(rad_params%xl(0:n_pft))
    allocate(rad_params%clump_fac(0:n_pft))

  end subroutine AllocateMLSolarParams
  
  !-----------------------------------------------------------------------
  
  subroutine SolarRadiation (num_filter, filter, mlcanopy_inst, pft)
    !
    ! !DESCRIPTION:
    ! Solar radiation transfer through canopy
    !
    ! !USES:
    use MLCanopyVarCon, only : pi => rpi
    use MLCanopyVarPar, only : numrad, ivis
    use MLCanopyVarCon, only : chil_max, chil_min, kb_max, J_to_umol
    use MLCanopyVarCtl, only : light_type, leaf_optics_type
    use MLCanopyVarPar, only : nlevmlcan, isun, isha
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filter                     ! Number of patches in filter
    integer, intent(in) :: filter(:)                      ! Patch filter
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    integer, intent(in) :: pft(:)                         ! plant functional type associated with the filter index (patch)
    
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                       ! Filter index
    integer  :: p                                        ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                        ! Column index for CLM g/l/c/p hierarchy
    integer  :: ic                                       ! Aboveground layer index
    integer  :: ib                                       ! Waveband index
    integer  :: j                                        ! Sky angle index
    real(r8) :: angle                                    ! Sky angle (5, 15, 25, 35, 45, 55, 65, 75, 85 degrees)
    real(r8) :: gdirj                                    ! Relative projected area of leaf elements in the direction of sky angle
    real(r8) :: wl, ws                                   ! Leaf and stem fraction of canopy layer

    real(r8) :: chil(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan) ! Departure of leaf angle from spherical orientation (-0.4 <= xl <= 0.6)
    real(r8) :: phi1(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan) ! Term in Ross-Goudriaan function for gdir
    real(r8) :: phi2(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan) ! Term in Ross-Goudriaan function for gdir
    real(r8) :: gdir(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan) ! Relative projected area of leaf in the direction of solar beam

    real(r8) :: clump_fac_ic(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan)    ! Foliage clumping index (-)

    real(r8) :: rho(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem reflectance
    real(r8) :: tau(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem transmittance
    real(r8) :: omega(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)  ! Leaf/stem scattering coefficient

    ! For two-stream radiation
    real(r8) :: asu                                                  ! Single scattering albedo
    real(r8) :: tmp0,tmp1,tmp2                                       ! Intermediate variables
    real(r8) :: avmu(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan)            ! Average inverse diffuse optical depth per unit leaf area
    real(r8) :: betad(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)  ! Upscatter parameter for diffuse radiation
    real(r8) :: betab(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)  ! Upscatter parameter for direct beam radiation
    !---------------------------------------------------------------------

    associate ( &
                                                        ! *** Input ***
    solar_zen   => mlcanopy_inst%solar_zen_forcing , &  ! Solar zenith angle (radians)
    ncan        => mlcanopy_inst%ncan_canopy       , &  ! Number of aboveground layers
    ntop        => mlcanopy_inst%ntop_canopy       , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy       , &  ! Index for bottom leaf layer
    dlai        => mlcanopy_inst%dlai_profile      , &  ! Canopy layer leaf area index (m2/m2)
    dsai        => mlcanopy_inst%dsai_profile      , &  ! Canopy layer stem area index (m2/m2)
    dpai        => mlcanopy_inst%dpai_profile      , &  ! Canopy layer plant area index (m2/m2)
                                                        ! *** Output ***
    fracsun     => mlcanopy_inst%fracsun_profile   , &  ! Canopy layer sunlit fraction (-)
    kb          => mlcanopy_inst%kb_profile        , &  ! Direct beam extinction coefficient (-)
    tb          => mlcanopy_inst%tb_profile        , &  ! Canopy layer transmittance of direct beam radiation (-)
    td          => mlcanopy_inst%td_profile        , &  ! Canopy layer transmittance of diffuse radiation (-)
    tbi         => mlcanopy_inst%tbi_profile       , &  ! Cumulative transmittance of direct beam onto canopy layer (-)
    apar        => mlcanopy_inst%apar_leaf         , &  ! Leaf absorbed PAR (umol photon/m2 leaf/s)
                                                        ! *** Output from radiation routines ***
    swveg       => mlcanopy_inst%swveg_canopy      , &  ! Absorbed solar radiation: vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun_canopy   , &  ! Absorbed solar radiation: sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha_canopy   , &  ! Absorbed solar radiation: shaded canopy (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy     , &  ! Albedo above canopy (-)
    swsoi       => mlcanopy_inst%swsoi_soil        , &  ! Absorbed solar radiation: ground (W/m2)
    swleaf      => mlcanopy_inst%swleaf_leaf       , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    swupw       => mlcanopy_inst%swupw_profile     , &  ! Upward diffuse solar flux above canopy layer (W/m2)
    swdwn       => mlcanopy_inst%swdwn_profile     , &  ! Downward diffuse solar flux above canopy layer (W/m2)
    swbeam      => mlcanopy_inst%swbeam_profile      &  ! Direct beam solar flux above canopy layer (W/m2)
    )

    ! Calculate canopy layer optical properties

    do fp = 1, num_filter
       p = filter(fp)

       ! Zero out variables for all layers

       chil(p,1:ncan(p)) = 0._r8
       phi1(p,1:ncan(p)) = 0._r8
       phi2(p,1:ncan(p)) = 0._r8
       gdir(p,1:ncan(p)) = 0._r8
       kb(p,1:ncan(p)) = 0._r8
       fracsun(p,1:ncan(p)) = 0._r8
       tb(p,1:ncan(p)) = 0._r8
       td(p,1:ncan(p)) = 0._r8
       tbi(p,0:ncan(p)) = 0._r8
       avmu(p,1:ncan(p)) = 0._r8
       clump_fac_ic(p,1:ncan(p)) = 0._r8

       rho(p,1:ncan(p),1:numrad) = 0._r8
       tau(p,1:ncan(p),1:numrad) = 0._r8
       omega(p,1:ncan(p),1:numrad) = 0._r8
       betab(p,1:ncan(p),1:numrad) = 0._r8
       betad(p,1:ncan(p),1:numrad) = 0._r8

       do ic = ntop(p), nbot(p) , -1

          ! Weight reflectance and transmittance by lai and sai and calculate
          ! leaf scattering coefficient

          wl = dlai(p,ic) / dpai(p,ic) ; ws = dsai(p,ic) / dpai(p,ic)
          do ib = 1, numrad
             select case (leaf_optics_type)
             case (0)
                 rho(p,ic,ib) = max(rad_params%rhol(pft(fp),ib)*wl + rad_params%rhos(pft(fp),ib)*ws, 1.e-06_r8)
                 tau(p,ic,ib) = max(rad_params%taul(pft(fp),ib)*wl + rad_params%taus(pft(fp),ib)*ws, 1.e-06_r8)
              case (1)
                call endrun (msg=' ERROR: SolarRadiation: need to specify vertical profile for rho & tau')
             end select
             omega(p,ic,ib) = rho(p,ic,ib) + tau(p,ic,ib)
          end do

          ! Leaf angle distribution

          select case (leaf_optics_type)
          case (0)
             chil(p,ic) = rad_params%xl(pft(fp))
          case (1)
             call endrun (msg=' ERROR: SolarRadiation: need to specify vertical profile for chil')
          end select
          chil(p,ic) = min(max(chil(p,ic), chil_min), chil_max)
          if (abs(chil(p,ic)) <= 0.01_r8) chil(p,ic) = 0.01_r8

          ! Terms in Ross-Goudriaan function for gdir

          phi1(p,ic) = 0.5_r8 - 0.633_r8 * chil(p,ic) - 0.330_r8 * chil(p,ic) * chil(p,ic)
          phi2(p,ic) = 0.877_r8 * (1._r8 - 2._r8 * phi1(p,ic))

          ! Relative projected area of leaf in the direction of solar beam

          gdir(p,ic) = phi1(p,ic) + phi2(p,ic) * cos(solar_zen(p))

          ! Direct beam extinction coefficient

          kb(p,ic) = gdir(p,ic) / cos(solar_zen(p))
          kb(p,ic) = min(kb(p,ic), kb_max)

          ! Foliage clumping factor

          select case (leaf_optics_type)
          case (0)
             clump_fac_ic(p,ic) = rad_params%clump_fac(pft(fp))
          case (1)
             call endrun (msg=' ERROR: SolarRadiation: need to specify vertical profile for clump_fac')
          end select

          ! Direct beam transmittance (tb) through a single layer

          tb(p,ic) = exp(-kb(p,ic) * dpai(p,ic) * clump_fac_ic(p,ic))

          ! Diffuse transmittance through a single layer (also needed for longwave
          ! radiation). Estimated for nine sky angles in increments of 10 degrees.

          td(p,ic) = 0._r8
          do j = 1, 9
             angle = (5._r8 + (j - 1) * 10._r8) * pi / 180._r8
             gdirj = phi1(p,ic) + phi2(p,ic) * cos(angle)
             td(p,ic) = td(p,ic) &
                      + exp(-gdirj / cos(angle) * dpai(p,ic) * clump_fac_ic(p,ic)) * sin(angle) * cos(angle)
          end do
          td(p,ic) = td(p,ic) * 2._r8 * (10._r8 * pi / 180._r8)

          ! Transmittance (tbi) of unscattered direct beam onto layer i

          if (ic == ntop(p)) then
             tbi(p,ic) = 1._r8
          else
             tbi(p,ic) = tbi(p,ic+1) * exp(-kb(p,ic+1) * dpai(p,ic+1) * clump_fac_ic(p,ic+1))
          end if

          ! Sunlit fraction of layer. Make sure fracsun > 0 and < 1.

          fracsun(p,ic) = tbi(p,ic) / (kb(p,ic) * dpai(p,ic)) &
                        * (1._r8 - exp(-kb(p,ic) * clump_fac_ic(p,ic) * dpai(p,ic)))

          if (fracsun(p,ic) <= 0._r8) then
             call endrun (msg=' ERROR: SolarRadiation: fracsun is too small')
          end if

          if ((1._r8 - fracsun(p,ic)) <= 0._r8) then
             call endrun (msg=' ERROR: SolarRadiation: fracsha is too small')
          end if

          !-----------------------------------------------------
          ! Special parameters for two-stream radiative transfer
          !-----------------------------------------------------

          ! avmu - average inverse diffuse optical depth per unit leaf area

          avmu(p,ic) = (1._r8 - phi1(p,ic)/phi2(p,ic) * log((phi1(p,ic)+phi2(p,ic))/phi1(p,ic))) / phi2(p,ic)

          do ib = 1, numrad

             ! betad - upscatter parameter for diffuse radiation

             betad(p,ic,ib) = 0.5_r8 / omega(p,ic,ib) * ( rho(p,ic,ib) + tau(p,ic,ib) &
                            + (rho(p,ic,ib)-tau(p,ic,ib)) * ((1._r8+chil(p,ic))/2._r8)**2 )

             ! betab - upscatter parameter for direct beam radiation

             tmp0 = gdir(p,ic) + phi2(p,ic) * cos(solar_zen(p))
             tmp1 = phi1(p,ic) * cos(solar_zen(p))
             tmp2 = 1._r8 - tmp1/tmp0 * log((tmp1+tmp0)/tmp1)
             asu = 0.5_r8 * omega(p,ic,ib) * gdir(p,ic) / tmp0 * tmp2
             betab(p,ic,ib) = (1._r8 + avmu(p,ic)*kb(p,ic)) / (omega(p,ic,ib)*avmu(p,ic)*kb(p,ic)) * asu
          end do

       end do

       ! Direct beam transmittance onto ground

       tbi(p,0) = tbi(p,nbot(p)) * exp(-kb(p,nbot(p)) * dpai(p,nbot(p)) * clump_fac_ic(p,nbot(p)))

    end do

    ! Calculate radiative transfer through canopy

    select case (light_type)
    case (1)
       call Norman (num_filter, filter, rho, tau, omega, mlcanopy_inst)
    case (2)
       call TwoStream (num_filter, filter, omega, avmu, betad, betab, clump_fac_ic, mlcanopy_inst)
    case default
       call endrun (msg=' ERROR: SolarRadiation: light_type not valid')
    end select

    ! APAR per unit sunlit and shaded leaf area (W/m2 -> umol/m2/s)

    do fp = 1, num_filter
       p = filter(fp)
       do ic = 1, ncan(p)
          apar(p,ic,isun) = swleaf(p,ic,isun,ivis) * J_to_umol
          apar(p,ic,isha) = swleaf(p,ic,isha,ivis) * J_to_umol
       end do
    end do

    end associate
  end subroutine SolarRadiation

  !-----------------------------------------------------------------------
  subroutine Norman (num_filter, filter, rho, tau, omega, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using Norman (1979)
    !
    ! !USES:
    use MLCanopyVarPar, only : numrad
    use MLCanopyVarPar, only : nlevmlcan, isun, isha
    use MLMathToolsMod, only : tridiag
    !
    ! !ARGUMENTS:
    implicit none
    integer , intent(in) :: num_filter                                           ! Number of patches in filter
    integer , intent(in) :: filter(:)                                            ! Patch filter
    real(r8), intent(in) :: rho(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem reflectance
    real(r8), intent(in) :: tau(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)    ! Leaf/stem transmittance
    real(r8), intent(in) :: omega(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,1:numrad)  ! Leaf/stem scattering coefficient
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                               ! Filter index
    integer  :: p                                                ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ic                                               ! Aboveground layer index
    integer  :: icm1                                             ! Layer below ic (ic-1)
    integer  :: ib                                               ! Waveband index
    real(r8) :: suminc                                           ! Incident radiation (W/m2)
    real(r8) :: sumref                                           ! Reflected radiation (W/m2)
    real(r8) :: sumabs                                           ! Absorbed radiation (W/m2)
    real(r8) :: err                                              ! Error check (W/m2)
    real(r8) :: refld                                            ! Term for diffuse radiation reflected by layer
    real(r8) :: trand                                            ! Term for diffuse radiation transmitted by layer
    integer  :: m                                                ! Index to the tridiagonal matrix
    real(r8) :: aic, bic                                         ! Intermediate terms for tridiagonal matrix
    real(r8) :: eic, fic                                         ! Intermediate terms for tridiagonal matrix
    integer, parameter :: neq = (nlevmlcan+1)*2                  ! Number of tridiagonal equations to solve
    real(r8) :: atri(neq), btri(neq)                             ! Entries in tridiagonal matrix
    real(r8) :: ctri(neq), dtri(neq)                             ! Entries in tridiagonal matrix
    real(r8) :: utri(neq)                                        ! Tridiagonal solution
    real(r8) :: swabsb                                           ! Absorbed direct beam solar flux for canopy layer (W/m2 ground)
    real(r8) :: swabsd                                           ! Absorbed diffuse solar flux for canopy layer (W/m2 ground)
    real(r8) :: swsun                                            ! Absorbed solar radiation, sunlit fraction of layer (W/m2 ground)
    real(r8) :: swsha                                            ! Absorbed solar radiation, shaded fraction of layer (W/m2 ground)
    !---------------------------------------------------------------------

    associate ( &
                                                      ! *** Input ***
    swskyb      => mlcanopy_inst%swskyb_forcing  , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing  , &  ! Atmospheric diffuse solar radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy     , &  ! Number of aboveground layers
    ntop        => mlcanopy_inst%ntop_canopy     , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy     , &  ! Index for bottom leaf layer
    albsoib     => mlcanopy_inst%albsoib_soil    , &  ! Direct beam albedo of ground (-)
    albsoid     => mlcanopy_inst%albsoid_soil    , &  ! Diffuse albedo of ground (-)
    dpai        => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    fracsun     => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
    tb          => mlcanopy_inst%tb_profile      , &  ! Canopy layer transmittance of direct beam radiation (-)
    td          => mlcanopy_inst%td_profile      , &  ! Canopy layer transmittance of diffuse radiation (-)
    tbi         => mlcanopy_inst%tbi_profile     , &  ! Cumulative transmittance of direct beam onto canopy layer (-)
                                                      ! *** Output ***
    swveg       => mlcanopy_inst%swveg_canopy    , &  ! Absorbed solar radiation: vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun_canopy , &  ! Absorbed solar radiation: sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha_canopy , &  ! Absorbed solar radiation: shaded canopy (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy   , &  ! Albedo above canopy (-)
    swsoi       => mlcanopy_inst%swsoi_soil      , &  ! Absorbed solar radiation: ground (W/m2)
    swleaf      => mlcanopy_inst%swleaf_leaf     , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    swupw       => mlcanopy_inst%swupw_profile   , &  ! Upward diffuse solar flux above canopy layer (W/m2)
    swdwn       => mlcanopy_inst%swdwn_profile   , &  ! Downward diffuse solar flux above canopy layer (W/m2)
    swbeam      => mlcanopy_inst%swbeam_profile    &  ! Direct beam solar flux above canopy layer (W/m2)
    )

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          ! Zero out radiative fluxes for all layers

          swbeam(p,0,ib) = 0._r8
          swupw(p,0,ib) = 0._r8
          swdwn(p,0,ib) = 0._r8

          do ic = 1, ncan(p)
             swbeam(p,ic,ib) = 0._r8
             swupw(p,ic,ib) = 0._r8
             swdwn(p,ic,ib) = 0._r8
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

          ! Set up and solve tridiagonal system of equations for upward and
          ! downward fluxes. There are two equations for each layer and the
          ! soil. These equations are referenced by "m". The first equation
          ! is the upward flux and the second equation is the downward flux.

          m = 0

          ! Soil: upward flux

          m = m + 1 
          atri(m) = 0._r8 
          btri(m) = 1._r8 
          ctri(m) = -albsoid(p,ib) 
          dtri(m) = swskyb(p,ib) * tbi(p,0) * albsoib(p,ib) 

          ! Soil: downward flux

          refld = (1._r8 - td(p,nbot(p))) * rho(p,nbot(p),ib) 
          trand = (1._r8 - td(p,nbot(p))) * tau(p,nbot(p),ib) + td(p,nbot(p)) 
          aic = refld - trand * trand / refld 
          bic = trand / refld 

          m = m + 1 
          atri(m) = -aic 
          btri(m) = 1._r8 
          ctri(m) = -bic 
          dtri(m) = swskyb(p,ib) * tbi(p,nbot(p)) * (1._r8 - tb(p,nbot(p))) * (tau(p,nbot(p),ib) - rho(p,nbot(p),ib) * bic)

          ! Leaf layers, excluding top layer

          do ic = nbot(p), ntop(p)-1

             ! Upward flux

             refld = (1._r8 - td(p,ic)) * rho(p,ic,ib)
             trand = (1._r8 - td(p,ic)) * tau(p,ic,ib) + td(p,ic)
             fic = refld - trand * trand / refld
             eic = trand / refld

             m = m + 1
             atri(m) = -eic
             btri(m) = 1._r8
             ctri(m) = -fic
             dtri(m) = swskyb(p,ib) * tbi(p,ic) * (1._r8 - tb(p,ic)) * (rho(p,ic,ib) - tau(p,ic,ib) * eic)

             ! Downward flux

             refld = (1._r8 - td(p,ic+1)) * rho(p,ic+1,ib)
             trand = (1._r8 - td(p,ic+1)) * tau(p,ic+1,ib) + td(p,ic+1)
             aic = refld - trand * trand / refld
             bic = trand / refld

             m = m + 1
             atri(m) = -aic
             btri(m) = 1._r8
             ctri(m) = -bic
             dtri(m) = swskyb(p,ib) * tbi(p,ic+1) * (1._r8 - tb(p,ic+1)) * (tau(p,ic+1,ib) - rho(p,ic+1,ib) * bic)

          end do

          ! Top canopy layer: upward flux

          ic = ntop(p)
          refld = (1._r8 - td(p,ic)) * rho(p,ic,ib)
          trand = (1._r8 - td(p,ic)) * tau(p,ic,ib) + td(p,ic)
          fic = refld - trand * trand / refld
          eic = trand / refld

          m = m + 1
          atri(m) = -eic
          btri(m) = 1._r8
          ctri(m) = -fic
          dtri(m) = swskyb(p,ib) * tbi(p,ic) * (1._r8 - tb(p,ic)) * (rho(p,ic,ib) - tau(p,ic,ib) * eic)

          ! Top canopy layer: downward flux

          m = m + 1
          atri(m) = 0._r8
          btri(m) = 1._r8
          ctri(m) = 0._r8
          dtri(m) = swskyd(p,ib)

          ! Solve tridiagonal system of equations for upward and downward fluxes

          call tridiag (atri, btri, ctri, dtri, utri, m)

          ! Now copy the solution (utri) to the upward (swupw) and downward (swdwn)
          ! fluxes for each layer
          ! swupw =  Upward diffuse flux above layer
          ! swdwn =  Downward diffuse flux onto layer

          m = 0

          ! Soil fluxes

          m = m + 1
          swupw(p,0,ib) = utri(m)
          m = m + 1
          swdwn(p,0,ib) = utri(m)

          ! Leaf layer fluxes

          do ic = nbot(p), ntop(p)
             m = m + 1
             swupw(p,ic,ib) = utri(m)
             m = m + 1
             swdwn(p,ic,ib) = utri(m)
          end do

       end do             ! end patch loop
    end do                ! end waveband loop

    ! Compute fluxes

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          ! Solar radiation absorbed by ground (soil)

          swbeam(p,0,ib) = tbi(p,0) * swskyb(p,ib)
          swabsb = swbeam(p,0,ib) * (1._r8 - albsoib(p,ib))
          swabsd = swdwn(p,0,ib) * (1._r8 - albsoid(p,ib))
          swsoi(p,ib) = swabsb + swabsd

          ! Leaf layer fluxes

          swveg(p,ib) = 0._r8
          swvegsun(p,ib) = 0._r8
          swvegsha(p,ib) = 0._r8

          do ic = nbot(p), ntop(p)

             ! Downward direct beam incident on layer and absorbed direct
             ! beam and diffuse for layer. Note the special case for first
             ! leaf layer, where the upward flux from below is from the ground.
             ! The ground is ic=0, but nbot-1 will not equal 0 if there are lower
             ! canopy layers without leaves. This coding accommodates open
             ! layers beneath the bottom of the canopy.

             swbeam(p,ic,ib) = tbi(p,ic) * swskyb(p,ib)
             swabsb = swbeam(p,ic,ib) * (1._r8 - tb(p,ic)) * (1._r8 - omega(p,ic,ib))
             if (ic == nbot(p)) then
                icm1 = 0
             else
                icm1 = ic - 1
             end if
             swabsd = (swdwn(p,ic,ib) + swupw(p,icm1,ib)) * (1._r8 - td(p,ic)) * (1._r8 - omega(p,ic,ib))

             ! Absorbed radiation for shaded and sunlit portions of layer

             swsha = swabsd * (1._r8 - fracsun(p,ic))
             swsun = swabsd * fracsun(p,ic) + swabsb 

             ! Per unit sunlit and shaded leaf area

             swleaf(p,ic,isun,ib) = swsun / (fracsun(p,ic) * dpai(p,ic))
             swleaf(p,ic,isha,ib) = swsha / ((1._r8 - fracsun(p,ic)) * dpai(p,ic))

             ! Sum solar radiation absorbed by vegetation and sunlit/shaded leaves

             swveg(p,ib) = swveg(p,ib) + (swabsb + swabsd)
             swvegsun(p,ib) = swvegsun(p,ib) + swsun
             swvegsha(p,ib) = swvegsha(p,ib) + swsha

          end do  

          ! Albedo

          suminc = swskyb(p,ib) + swskyd(p,ib)
          if (suminc > 0._r8) then
             albcan(p,ib) = swupw(p,ntop(p),ib) / suminc
          else
             albcan(p,ib) = 0._r8
          end if

          ! Conservation check for total radiation balance: absorbed = incoming - outgoing

          sumref = albcan(p,ib) * (swskyb(p,ib) + swskyd(p,ib))
          sumabs = suminc - sumref
          err = sumabs - (swveg(p,ib) + swsoi(p,ib))
          if (abs(err) > 1.e-03_r8) then
             call endrun (msg='ERROR: Norman: total solar conservation error')
          end if

          ! Sunlit and shaded absorption

          err = (swvegsun(p,ib) + swvegsha(p,ib)) - swveg(p,ib) 
          if (abs(err) > 1.e-03_r8) then
             call endrun (msg='ERROR: Norman: sunlit/shade solar conservation error')
          end if

       end do             ! end patch loop
    end do                ! end waveband loop

    end associate
  end subroutine Norman

  !-----------------------------------------------------------------------
  subroutine TwoStream (num_filter, filter, omega, avmu, betad, betab, clump_fac_ic, mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute solar radiation transfer through canopy using the two-stream
    ! approximation. This solution is integrated over each layer and so
    ! does not require thin canopy layers. Also, optical properties
    ! can vary with depth in the canopy.
    !
    ! !USES:
    use MLCanopyVarPar, only : numrad
    use MLCanopyVarPar, only : nlevmlcan, isun, isha
    !
    ! !ARGUMENTS:
    implicit none
    integer,  intent(in) :: num_filter                                        ! Number of patches in filter
    integer,  intent(in) :: filter(:)                                         ! Patch filter
    real(r8), intent(in) :: omega(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad) ! Leaf/stem scattering coefficient
    real(r8), intent(in) :: avmu(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan)         ! Average inverse diffuse optical depth per unit leaf area
    real(r8), intent(in) :: betad(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad) ! Upscatter parameter for diffuse radiation
    real(r8), intent(in) :: betab(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad) ! Upscatter parameter for direct beam radiation
    real(r8), intent(in) :: clump_fac_ic(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan) ! Foliage clumping index (-)
    type(mlcanopy_type), intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: fp                                                            ! Filter index
    integer  :: p                                                             ! Patch index for CLM g/l/c/p hierarchy
    integer  :: ib                                                            ! Waveband index
    integer  :: ic                                                            ! Aboveground layer index
    real(r8) :: b,c,d,h,u,v                                                   ! Intermediate two-stream parameters
    real(r8) :: g1,g2                                                         ! Intermediate two-stream parameters
    real(r8) :: s1,s2                                                         ! Intermediate two-stream parameters
    real(r8) :: num1,num2                                                     ! Intermediate two-stream parameters
    real(r8) :: den1,den2                                                     ! Intermediate two-stream parameters
    real(r8) :: n1b,n2b                                                       ! Two-stream parameters
    real(r8) :: n1d,n2d                                                       ! Two-stream parameters
    real(r8) :: a1b,a2b                                                       ! Two-stream parameters
    real(r8) :: a1d,a2d                                                       ! Two-stream parameters
    real(r8) :: dir,dif                                                       ! Direct beam and diffuse fluxes (W/m2)
    real(r8) :: sun,sha                                                       ! Sunlit and shaded fluxes (W/m2)
    real(r8) :: suminc                                                        ! Incident radiation (W/m2)
    real(r8) :: sumref                                                        ! Reflected radiation (W/m2)
    real(r8) :: sumabs                                                        ! Absorbed radiation (W/m2)
    real(r8), parameter :: unitb = 1._r8                                      ! Unit direct beam radiation (W/m2)
    real(r8), parameter :: unitd = 1._r8                                      ! Unit diffuse radiation (W/m2)

    real(r8) :: iupwb0(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)            ! Direct beam flux scattered upward (reflected) above canopy layer (W/m2)
    real(r8) :: iupwb(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)             ! Direct beam flux scattered upward at the canopy layer depth (W/m2)
    real(r8) :: idwnb(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)             ! Direct beam flux scattered downward below canopy layer (W/m2)
    real(r8) :: iabsb(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)             ! Direct beam flux absorbed by canopy layer (W/m2)
    real(r8) :: iabsbb(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)            ! Unscattered component of iabsb (W/m2)
    real(r8) :: iabsbs(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)            ! Scattered component of iabsb (W/m2)
    real(r8) :: iabsb_sun(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)         ! Sunlit component of iabsb (W/m2)
    real(r8) :: iabsb_sha(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)         ! Shaded component of iabsb (W/m2)

    real(r8) :: iupwd0(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)            ! Diffuse flux scattered upward (reflected) above canopy layer (W/m2)
    real(r8) :: iupwd(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)             ! Diffuse flux scattered upward at the canopy layer depth (W/m2)
    real(r8) :: idwnd(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)             ! Diffuse flux scattered downward below canopy layer (W/m2)
    real(r8) :: iabsd(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)             ! Diffuse flux absorbed by canopy layer (W/m2)
    real(r8) :: iabsd_sun(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)         ! Sunlit component of iabsd (W/m2)
    real(r8) :: iabsd_sha(mlcanopy_inst%begp:mlcanopy_inst%endp,1:nlevmlcan,numrad)         ! Shaded component of iabsd (W/m2)

    real(r8) :: albb_below(mlcanopy_inst%begp:mlcanopy_inst%endp,numrad)                    ! Direct beam albedo for canopy layer below current layer
    real(r8) :: albd_below(mlcanopy_inst%begp:mlcanopy_inst%endp,numrad)                    ! Diffuse albedo for canopy layer below current layer
    !---------------------------------------------------------------------

    associate ( &
                                                      ! *** Input ***
    swskyb      => mlcanopy_inst%swskyb_forcing  , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd      => mlcanopy_inst%swskyd_forcing  , &  ! Atmospheric diffuse solar radiation (W/m2)
    ncan        => mlcanopy_inst%ncan_canopy     , &  ! Number of aboveground layers
    ntop        => mlcanopy_inst%ntop_canopy     , &  ! Index for top leaf layer
    nbot        => mlcanopy_inst%nbot_canopy     , &  ! Index for bottom leaf layer
    albsoib     => mlcanopy_inst%albsoib_soil    , &  ! Direct beam albedo of ground (-)
    albsoid     => mlcanopy_inst%albsoid_soil    , &  ! Diffuse albedo of ground (-)
    dpai        => mlcanopy_inst%dpai_profile    , &  ! Canopy layer plant area index (m2/m2)
    fracsun     => mlcanopy_inst%fracsun_profile , &  ! Canopy layer sunlit fraction (-)
    kb          => mlcanopy_inst%kb_profile      , &  ! Direct beam extinction coefficient (-)
    tbi         => mlcanopy_inst%tbi_profile     , &  ! Cumulative transmittance of direct beam onto canopy layer (-)
                                                      ! *** Output ***
    swveg       => mlcanopy_inst%swveg_canopy    , &  ! Absorbed solar radiation: vegetation (W/m2)
    swvegsun    => mlcanopy_inst%swvegsun_canopy , &  ! Absorbed solar radiation: sunlit canopy (W/m2)
    swvegsha    => mlcanopy_inst%swvegsha_canopy , &  ! Absorbed solar radiation: shaded canopy (W/m2)
    albcan      => mlcanopy_inst%albcan_canopy   , &  ! Albedo above canopy (-)
    swsoi       => mlcanopy_inst%swsoi_soil      , &  ! Absorbed solar radiation: ground (W/m2)
    swleaf      => mlcanopy_inst%swleaf_leaf     , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    swupw       => mlcanopy_inst%swupw_profile   , &  ! Upward diffuse solar flux above canopy layer (W/m2)
    swdwn       => mlcanopy_inst%swdwn_profile   , &  ! Downward diffuse solar flux above canopy layer (W/m2)
    swbeam      => mlcanopy_inst%swbeam_profile    &  ! Direct beam solar flux above canopy layer (W/m2)
    )

    ! Zero out radiative fluxes for all layers

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          swbeam(p,0,ib) = 0._r8
          swupw(p,0,ib) = 0._r8
          swdwn(p,0,ib) = 0._r8

          do ic = 1, ncan(p)
             swbeam(p,ic,ib) = 0._r8
             swupw(p,ic,ib) = 0._r8
             swdwn(p,ic,ib) = 0._r8
             swleaf(p,ic,isun,ib) = 0._r8
             swleaf(p,ic,isha,ib) = 0._r8
          end do

       end do
    end do


    ! Calculate radiative transfer in two steps:
    ! (1) Working from bottom to top of canopy, calculate radiative fluxes for
    !     a unit of direct beam radiation and a unit of diffuse radiation at
    !     the top of the layer. Boundary conditions are the albedos (direct,
    !     diffuse) for the immediate layer below. For the bottom of the canopy
    !     (nbot), these are the soil albedos. For all other layers, these are
    !     the upward diffuse fluxes above the lower layer.
    ! (2) Then working from top to bottom of canopy, calculate the fluxes
    !     incident on a layer and the absorption by sunlit and shaded leaves.

    do ib = 1, numrad
       do fp = 1, num_filter
          p = filter(fp)

          ! Initialize albedos below current layer

          albb_below(p,ib) = albsoib(p,ib) ; albd_below(p,ib) = albsoid(p,ib)

          ! Layer fluxes, from bottom to top

          do ic = nbot(p), ntop(p)

             ! Common terms

             b = (1._r8 - (1._r8 - betad(p,ic,ib)) * omega(p,ic,ib)) / avmu(p,ic)
             c = betad(p,ic,ib) * omega(p,ic,ib) / avmu(p,ic)
             h = sqrt(b*b - c*c)
             u = (h - b - c) / (2._r8 * h)
             v = (h + b + c) / (2._r8 * h)
             d = omega(p,ic,ib) * kb(p,ic) * unitb / (h*h - kb(p,ic)*kb(p,ic))
             g1 = (betab(p,ic,ib) * kb(p,ic) - b * betab(p,ic,ib) - c * (1._r8 - betab(p,ic,ib))) * d
             g2 = ((1._r8 - betab(p,ic,ib)) * kb(p,ic) + c * betab(p,ic,ib) + b * (1._r8 - betab(p,ic,ib))) * d
             s1 = exp(-h * clump_fac_ic(p,ic) * dpai(p,ic))
             s2 = exp(-kb(p,ic) * clump_fac_ic(p,ic) * dpai(p,ic))

             ! Terms for direct beam radiation

             num1 = v * (g1 + g2 * albd_below(p,ib) + albb_below(p,ib) * unitb) * s2
             num2 = g2 * (u + v * albd_below(p,ib)) * s1
             den1 = v * (v + u * albd_below(p,ib)) / s1
             den2 = u * (u + v * albd_below(p,ib)) * s1
             n2b = (num1 - num2) / (den1 - den2)
             n1b = (g2 - n2b * u) / v

             a1b = -g1 *      (1._r8 - s2*s2) / (2._r8 * kb(p,ic)) &
                 +  n1b * u * (1._r8 - s2*s1) / (kb(p,ic) + h) + n2b * v * (1._r8 - s2/s1) / (kb(p,ic) - h)
             a2b =  g2 *      (1._r8 - s2*s2) / (2._r8 * kb(p,ic)) &
                 -  n1b * v * (1._r8 - s2*s1) / (kb(p,ic) + h) - n2b * u * (1._r8 - s2/s1) / (kb(p,ic) - h)
             a1b = a1b * tbi(p,ic)     ! To account for fracsun in multilayer canopy
             a2b = a2b * tbi(p,ic)     ! To account for fracsun in multilayer canopy

             ! Direct beam radiative fluxes

             iupwb0(p,ic,ib) = -g1 + n1b * u + n2b * v
             iupwb(p,ic,ib) = -g1 * s2 + n1b * u * s1 + n2b * v / s1
             idwnb(p,ic,ib) =  g2 * s2 - n1b * v * s1 - n2b * u / s1
             iabsb(p,ic,ib) = unitb * (1._r8 - s2) - iupwb0(p,ic,ib) + iupwb(p,ic,ib) - idwnb(p,ic,ib)
             iabsbb(p,ic,ib) = (1._r8 - omega(p,ic,ib)) * unitb * (1._r8 - s2)
             iabsbs(p,ic,ib) = omega(p,ic,ib) * unitb * (1._r8 - s2) - iupwb0(p,ic,ib) + iupwb(p,ic,ib) - idwnb(p,ic,ib)
             iabsb_sun(p,ic,ib) = (1._r8-omega(p,ic,ib)) * ((1._r8-s2) * unitb + clump_fac_ic(p,ic) / avmu(p,ic) * (a1b+a2b))
             iabsb_sha(p,ic,ib) = iabsb(p,ic,ib) - iabsb_sun(p,ic,ib)

             ! Compare my integration for sunlit leaves with that used in the code
             ! sun = (1._r8-omega(p,ic,ib)) * (1._r8-s2) * unitb
             ! sun = sun + 0.5_r8 * omega(p,ic,ib) * clump_fac_ic(p,ic) * (1._r8-s2*s2) * unitb * tbi(p,ic)
             ! sun = sun + 0.5_r8 * clump_fac_ic(p,ic) * g1 * (1._r8-s2*s2) * tbi(p,ic)
             ! sun = sun + 0.5_r8 * clump_fac_ic(p,ic) * g2 * (1._r8-s2*s2) * tbi(p,ic)
             ! sun = sun + (-n1b * u * (1._r8 - s2*s1) / (kb(p,ic) + h) + n2b * v * (1._r8 - s2/s1) / (kb(p,ic) - h)) * h * clump_fac_ic(p,ic) * tbi(p,ic)
             ! sun = sun - ( n1b * v * (1._r8 - s2*s1) / (kb(p,ic) + h) - n2b * u * (1._r8 - s2/s1) / (kb(p,ic) - h)) * h * clump_fac_ic(p,ic) * tbi(p,ic)
             ! if (abs(iabsb_sun(p,ic,ib)-sun) > 1.e-10_r8) then
             !    write (6,*) iabsb_sun(p,ic,ib), sun
             !    stop
             ! end if

             ! Terms for diffuse radiation

             num1 = unitd * (u + v * albd_below(p,ib)) * s1
             den1 = v * (v + u * albd_below(p,ib)) / s1
             den2 = u * (u + v * albd_below(p,ib)) * s1
             n2d = num1 / (den1 - den2)
             n1d = -(unitd + n2d * u) / v

             a1d =  n1d * u * (1._r8 - s2*s1) / (kb(p,ic) + h) + n2d * v * (1._r8 - s2/s1) / (kb(p,ic) - h)
             a2d = -n1d * v * (1._r8 - s2*s1) / (kb(p,ic) + h) - n2d * u * (1._r8 - s2/s1) / (kb(p,ic) - h)
             a1d = a1d * tbi(p,ic)     ! To account for fracsun in multilayer canopy
             a2d = a2d * tbi(p,ic)     ! To account for fracsun in multilayer canopy

             ! Diffuse radiative fluxes

             iupwd0(p,ic,ib) = n1d * u + n2d * v
             iupwd(p,ic,ib) =  n1d * u * s1 + n2d * v / s1
             idwnd(p,ic,ib) = -n1d * v * s1 - n2d * u / s1
             iabsd(p,ic,ib) = unitd - iupwd0(p,ic,ib) + iupwd(p,ic,ib) - idwnd(p,ic,ib)
             iabsd_sun(p,ic,ib) = (1._r8 - omega(p,ic,ib)) * clump_fac_ic(p,ic) / avmu(p,ic) * (a1d + a2d)
             iabsd_sha(p,ic,ib) = iabsd(p,ic,ib) - iabsd_sun(p,ic,ib)

             ! Compare my integration for sunlit leaves with that used in the code
             ! sun = -n1d * u * (1._r8 - s2*s1) / (kb(p,ic) + h) + n2d * v * (1._r8 - s2/s1) / (kb(p,ic) - h) &
             !       -n1d * v * (1._r8 - s2*s1) / (kb(p,ic) + h) + n2d * u * (1._r8 - s2/s1) / (kb(p,ic) - h)
             ! sun = sun * h * clump_fac_ic(p,ic) * tbi(p,ic)
             ! if (abs(iabsd_sun(p,ic,ib)-sun) > 1.e-10_r8) then
             !    write (6,*) iabsd_sun(p,ic,ib), sun
             !    stop
             ! end if

             ! Update albedos to be used for the next layer

             albb_below(p,ib) = iupwb0(p,ic,ib) ; albd_below(p,ib) = iupwd0(p,ic,ib)

          end do

          ! Now working from top to bottom of canopy, calculate the fluxes
          ! incident on a layer and the absorption by sunlit and shades leaves

          dir = swskyb(p,ib) ; dif = swskyd(p,ib) 
          do ic = ntop(p), nbot(p), -1

             ! Downward and upward fluxes above layer

             swbeam(p,ic,ib) = dir
             swdwn(p,ic,ib) = dif
             swupw(p,ic,ib) = iupwd0(p,ic,ib) * dif + iupwb0(p,ic,ib) * dir

             ! Absorption by canopy layer (W/m2 leaf)

             sun = (iabsb_sun(p,ic,ib) * dir + iabsd_sun(p,ic,ib) * dif) / (fracsun(p,ic) * dpai(p,ic))
             sha = (iabsb_sha(p,ic,ib) * dir + iabsd_sha(p,ic,ib) * dif) / ((1._r8 - fracsun(p,ic)) * dpai(p,ic))
             swleaf(p,ic,isun,ib) = sun ; swleaf(p,ic,isha,ib) = sha

             ! Diffuse and direct beam radiation incident on top of lower layer

             dif = dir * idwnb(p,ic,ib) + dif * idwnd(p,ic,ib)
             dir = dir * exp(-kb(p,ic) * clump_fac_ic(p,ic) * dpai(p,ic))

          end do

          ! Downward and upward fluxes above ground (soil)

          swbeam(p,0,ib) = dir
          swdwn(p,0,ib) = dif
          swupw(p,0,ib) = albsoid(p,ib) * dif + albsoib(p,ib) * dir

          ! Solar radiation absorbed by ground (soil)

           swsoi(p,ib) = dir * (1._r8 - albsoib(p,ib)) + dif * (1._r8 - albsoid(p,ib))

          ! Canopy albedo

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumref = iupwb0(p,ntop(p),ib) * swskyb(p,ib) + iupwd0(p,ntop(p),ib) * swskyd(p,ib)
          if (suminc > 0._r8) then
             albcan(p,ib) = sumref / suminc
          else
             albcan(p,ib) = 0._r8
          end if

          ! Sum canopy absorption (W/m2 ground) using leaf fluxes per unit sunlit
          ! and shaded leaf area (W/m2 leaf)

          swveg(p,ib) = 0._r8
          swvegsun(p,ib) = 0._r8
          swvegsha(p,ib) = 0._r8

          do ic = nbot(p), ntop(p)
             sun = swleaf(p,ic,isun,ib) * fracsun(p,ic) * dpai(p,ic)
             sha = swleaf(p,ic,isha,ib) * (1._r8 - fracsun(p,ic)) * dpai(p,ic)
             swveg(p,ib) = swveg(p,ib) + (sun +  sha)
             swvegsun(p,ib) = swvegsun(p,ib) + sun
             swvegsha(p,ib) = swvegsha(p,ib) + sha
          end do

          ! Conservation check: total incident = total reflected + total absorbed

          suminc = swskyb(p,ib) + swskyd(p,ib)
          sumref = albcan(p,ib) * suminc
          sumabs = swveg(p,ib) + swsoi(p,ib)

          if (abs(suminc - (sumabs+sumref)) >= 1.e-06_r8) then
             call endrun (msg='ERROR: TwoStream: total solar radiation conservation error')
          end if

       end do      ! end grid point loop
    end do         ! end waveband loop

    end associate
  end subroutine TwoStream

end module MLSolarRadiationMod
