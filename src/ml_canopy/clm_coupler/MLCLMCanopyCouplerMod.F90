module MLCLMCanopyCouplerMod
  
  use abortutils, only : endrun
  use clm_varcon, only : ispval, spval
  use clm_varpar, only : nlevgrnd, numrad
  use decompMod , only : bounds_type
  use shr_kind_mod, only : r8 => shr_kind_r8use

  use clm_varcon, only : pi => rpi
  use clm_varpar, only : numrad, ivis
  use MLCanopyFluxesType, only : mlcanopy_type

  use PatchType, only : patch
  use pftconMod, only : pftcon
  
  use decompMod           , only : bounds_type
  use atm2lndType         , only : atm2lnd_type
  use CanopyStateType     , only : canopystate_type
  use ColumnType          , only : col
  use EnergyFluxType      , only : energyflux_type
  use FrictionVelocityMod , only : frictionvel_type
  use GridcellType        , only : grc
  use PatchType           , only : patch
  use pftconMod           , only : pftcon
  use SoilStateType       , only : soilstate_type
  use SolarAbsorbedType   , only : solarabs_type
  use SurfaceAlbedoType   , only : surfalb_type
  use TemperatureType     , only : temperature_type
  use WaterFluxType       , only : waterflux_type
  use WaterStateType      , only : waterstate_type
  use MLCanopyFluxesType  , only : mlcanopy_type
  
  implicit none
  save
  private

contains

  subroutine TransferMLCLMParams()

    use MLSolarRadiationMod, only : rad_params
    use MLSolarRadiationMod, only : AllocateRadParams
    
    integer :: ft,ib  ! loop indices

    call AllocateRadParams(numpft)

    do ft = 1,numpft
       do ib = 1,numrad

          rad_params%rhol(ib,ft) = EDPftvarcon_inst%rhol(ft,ib)
          rad_params%rhos(ib,ft) = EDPftvarcon_inst%rhos(ft,ib)
          rad_params%taul(ib,ft) = EDPftvarcon_inst%taul(ft,ib)
          rad_params%taus(ib,ft) = EDPftvarcon_inst%taus(ft,ib)

       end do
       rad_params%xl(ft) = EDPftvarcon_inst%xl(ft)
       rad_params%clumping_index(ft) = EDPftvarcon_inst%clumping_index(ft)
    end do
    
    call RadParamPrep()

    return
  end subroutine TransferCLMMLParams
  
  

  !-----------------------------------------------------------------------

  subroutine Init (mlcanopy,bounds)
    !
    ! !DESCRIPTION:
    !
    ! Initialization of the data type. Allocate data, setup variables
    ! for history output, and initialize values needed for a cold-start
    !
    type(mlcanopy_type) :: mlcanopy
    type(bounds_type), intent(in) :: bounds

    mlcanopy%begp = bounds%begp
    mlcanopy%endp = bounds%endp
    
    call mlcanopy%InitAllocate()
    call InitMLHistory(mlcanopy,bounds)
    call mlcanopy%InitCold()

  end subroutine Init
  
  !-----------------------------------------------------------------------

  subroutine InitMLHistory (mlcanopy, bounds)
    !
    ! !DESCRIPTION:
    ! Setup the fields that can be output on history files
    !
    ! !USES:
    use histFileMod, only: hist_addfld1d, hist_addfld2d
    !
    ! !ARGUMENTS:
    type(mlcanopy_type) :: mlcanopy
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    !---------------------------------------------------------------------

    begp = bounds%begp ; endp= bounds%endp

    mlcanopy%gppveg_canopy(begp:endp) = spval
    call hist_addfld1d (fname='GPP_ML', units='umol/m2s', &
         avgflag='A', long_name='Gross primary production', &
         ptr_patch=mlcanopy%gppveg_canopy, set_lake=spval, set_urb=spval)

    mlcanopy%lwp_mean_profile(begp:endp,1:nlevmlcan) = spval
    call hist_addfld2d (fname='LWP_ML', units='MPa', type2d='nlevmlcan', &
         avgflag='A', long_name='Weighted maean leaf water potential of canopy layer', &
         ptr_patch=mlcanopy%lwp_mean_profile, set_lake=spval, set_urb=spval)

  end subroutine InitCLMHistory

    !-----------------------------------------------------------------------
  subroutine Restart (mlcanopy, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file
    !
    ! !USES:
    use ncdio_pio, only : file_desc_t, ncd_defvar, ncd_io, ncd_double, ncd_int, ncd_inqvdlen
    use restUtilMod, only : restartvar
    !
    ! !ARGUMENTS:
    type(mlcanopy_type) :: mlcanopy
    type(bounds_type), intent(in)    :: bounds
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    logical :: readvar      ! determine if variable is on initial file
    !---------------------------------------------------------------------

    ! Example for 1-d patch variable

    call restartvar(ncid=ncid, flag=flag, varname='taf_ml', xtype=ncd_double,  &
       dim1name='pft', long_name='air temperature at canopy top', units='K', &
       interpinic_flag='interp', readvar=readvar, data=mlcanopy%taf_canopy)

    ! Example for 2-d patch variable

    call restartvar(ncid=ncid, flag=flag, varname='lwp_ml', xtype=ncd_double,  &
       dim1name='pft', dim2name='nlevmlcan', switchdim=.true., &
       long_name='leaf water potential of canopy layer', units='MPa', &
       interpinic_flag='interp', readvar=readvar, data=mlcanopy%lwp_mean_profile)

  end subroutine Restart

  !-----------------------------------------------------------------------

  subroutine MLCLMCanopyFluxes (bounds, num_exposedvegp, filter_exposedvegp, &
       atm2lnd_inst, canopystate_inst, soilstate_inst, temperature_inst, waterstate_inst, &
       waterflux_inst, energyflux_inst, frictionvel_inst, surfalb_inst, solarabs_inst, &
       mlcanopy_inst)
    !
    ! !DESCRIPTION:
    ! Compute fluxes for sunlit and shaded leaves at each level
    ! and for soil surface
    !
    ! !USES:
    use clm_time_manager, only : get_nstep, get_step_size, get_curr_calday
    use clm_varcon, only : grav, pi => rpi, spval
    use clm_varorb, only : eccen, obliqr, lambm0, mvelpp
    use clm_varpar, only : ivis, inir
    use shr_orb_mod, only : shr_orb_decl, shr_orb_cosz
    use spmdMod, only : masterproc
    use MLclm_varcon, only : mmh2o, mmdry, cpd, cpw, rgas, wind_forc_min, lapse_rate
    use MLclm_varctl, only : mlcan_to_clm, dtime_substep, ml_vert_init, fracdir
    use MLclm_varpar, only : isun, isha, nlevmlcan, nleaf
    use MLCanopyNitrogenProfileMod, only : CanopyNitrogenProfile
    use MLCanopyTurbulenceMod, only : CanopyTurbulence
    use MLCanopyWaterMod, only : CanopyInterception, CanopyEvaporation
    use MLinitVerticalMod, only : initVerticalProfiles, initVerticalStructure
    use MLLeafBoundaryLayerMod, only : LeafBoundaryLayer
    use MLLeafHeatCapacityMod, only : LeafHeatCapacity
    use MLLeafPhotosynthesisMod, only : LeafPhotosynthesis
    use MLLongwaveRadiationMod, only : LongwaveRadiation
    use MLPlantHydraulicsMod, only: SoilResistance, PlantResistance, LeafWaterPotential
    use MLSolarRadiationMod, only: SolarRadiation
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_exposedvegp           ! Number of non-snow-covered veg patches in CLM patch filter
    integer, intent(in) :: filter_exposedvegp(:)     ! CLM patch filter for non-snow-covered vegetation

    type(atm2lnd_type)     , intent(in)    :: atm2lnd_inst
    type(canopystate_type) , intent(inout) :: canopystate_inst
    type(soilstate_type)   , intent(inout) :: soilstate_inst
    type(temperature_type) , intent(inout) :: temperature_inst
    type(waterstate_type)  , intent(inout) :: waterstate_inst
    type(waterflux_type)   , intent(inout) :: waterflux_inst
    type(energyflux_type)  , intent(inout) :: energyflux_inst
    type(frictionvel_type) , intent(inout) :: frictionvel_inst
    type(surfalb_type)     , intent(inout) :: surfalb_inst
    type(solarabs_type)    , intent(inout) :: solarabs_inst
    type(mlcanopy_type)    , intent(inout) :: mlcanopy_inst
    !
    ! !LOCAL VARIABLES:
    integer  :: num_mlcan                               ! Number of vegetated patches for multilayer canopy
    integer  :: fp                                      ! Filter index
    integer  :: p                                       ! Patch index for CLM g/l/c/p hierarchy
    integer  :: c                                       ! Column index for CLM g/l/c/p hierarchy
    integer  :: g                                       ! Gridcell index for CLM g/l/c/p hierarchy
    integer  :: ic                                      ! Aboveground layer index
    integer  :: nstep                                   ! Current model time step number
    integer  :: num_sub_steps                           ! Number of sub-time steps
    integer  :: niter                                   ! Current sub-time step
    real(r8) :: dtime                                   ! Model time step (s)
    real(r8) :: caldaym1                                ! Calendar day for zenith angle (1.000 on 0Z January 1 of current year)
    real(r8) :: declinm1                                ! Solar declination angle for zenith angle (radians)
    real(r8) :: eccf                                    ! Earth orbit eccentricity factor
    real(r8) :: coszen                                  ! Cosine solar zenith angle
    real(r8) :: lat, lon                                ! Latitude and longitude (radians)
    real(r8) :: totpai                                  ! Canopy lai+sai for error check (m2/m2)
    real(r8) :: totrad                                  ! Direct beam + diffuse solar radiation (W/m2)

    ! These are used to accumulate flux variables over the model sub-time step.
    ! The last dimension is the number of variables

    real(r8) :: flux_accumulator(bounds%begp:bounds%endp,nvar1d)                          ! Single-level fluxes
    real(r8) :: flux_accumulator_profile(bounds%begp:bounds%endp,1:nlevmlcan,nvar2d)      ! Multi-level profile fluxes
    real(r8) :: flux_accumulator_leaf(bounds%begp:bounds%endp,1:nlevmlcan,1:nleaf,nvar3d) ! Multi-level leaf fluxes
    !---------------------------------------------------------------------

    ! Variables used in this subroutine. See README.txt for a complete list of
    ! required CLM variables and also required CLM files. See MLCanopyFluxesType.F90
    ! for a complete list of multilayer canopy variables.

    associate ( &
                                                                  ! *** CLM variables ***
    forc_u         => atm2lnd_inst%forc_u_grc                , &  ! INPUT: Atmospheric wind speed in east direction (m/s)
    forc_v         => atm2lnd_inst%forc_v_grc                , &  ! INPUT: Atmospheric wind speed in north direction (m/s)
    forc_pco2      => atm2lnd_inst%forc_pco2_grc             , &  ! INPUT: Atmospheric CO2 partial pressure (Pa)
    forc_po2       => atm2lnd_inst%forc_po2_grc              , &  ! INPUT: Atmospheric O2 partial pressure (Pa)
    forc_solad     => atm2lnd_inst%forc_solad_grc            , &  ! INPUT: Atmospheric direct beam radiation (W/m2)
    forc_solai     => atm2lnd_inst%forc_solai_grc            , &  ! INPUT: Atmospheric diffuse radiation (W/m2)
    forc_t         => atm2lnd_inst%forc_t_downscaled_col     , &  ! INPUT: Atmospheric temperature (K)
    forc_q         => atm2lnd_inst%forc_q_downscaled_col     , &  ! INPUT: Atmospheric specific humidity (kg/kg)
    forc_pbot      => atm2lnd_inst%forc_pbot_downscaled_col  , &  ! INPUT: Atmospheric pressure (Pa)
    forc_lwrad     => atm2lnd_inst%forc_lwrad_downscaled_col , &  ! INPUT: Atmospheric longwave radiation (W/m2)
    forc_rain      => atm2lnd_inst%forc_rain_downscaled_col  , &  ! INPUT: Rainfall rate (mm/s)
    forc_snow      => atm2lnd_inst%forc_snow_downscaled_col  , &  ! INPUT: Snowfall rate (mm/s)
    elai           => canopystate_inst%elai_patch            , &  ! INPUT: Leaf area index of canopy (m2/m2)
    esai           => canopystate_inst%esai_patch            , &  ! INPUT: Stem area index of canopy (m2/m2)
    snl            => col%snl                                , &  ! INPUT: Number of snow layers
    z              => col%z                                  , &  ! INPUT: Soil layer depth (m)
    zi             => col%zi                                 , &  ! INPUT: Soil layer depth at layer interface (m)
    soilresis      => soilstate_inst%soilresis_col           , &  ! INPUT: Soil evaporative resistance (s/m)
    thk            => soilstate_inst%thk_col                 , &  ! INPUT: Soil layer thermal conductivity (W/m/K)
    smp_l          => soilstate_inst%smp_l_col               , &  ! INPUT: Soil layer matric potential (mm)
    albgrd         => surfalb_inst%albgrd_col                , &  ! INPUT: Direct beam albedo of ground (soil)
    albgri         => surfalb_inst%albgri_col                , &  ! INPUT: Diffuse albedo of ground (soil)
    t_a10_patch    => temperature_inst%t_a10_patch           , &  ! INPUT: 10-day running mean of the 2-m temperature (K)
    t_soisno       => temperature_inst%t_soisno_col          , &  ! INPUT: Soil temperature (K)
    eflx_lh_tot    => energyflux_inst%eflx_lh_tot_patch      , &  ! OUTPUT: patch total latent heat flux (W/m2)
    eflx_sh_tot    => energyflux_inst%eflx_sh_tot_patch      , &  ! OUTPUT: patch total sensible heat flux (W/m2)
    eflx_lwrad_out => energyflux_inst%eflx_lwrad_out_patch   , &  ! OUTPUT: patch emitted infrared (longwave) radiation (W/m2)
    taux           => energyflux_inst%taux_patch             , &  ! OUTPUT: patch wind (shear) stress: e-w (kg/m/s2)
    tauy           => energyflux_inst%tauy_patch             , &  ! OUTPUT: patch wind (shear) stress: n-s (kg/m/s2)
    fv             => frictionvel_inst%fv_patch              , &  ! OUTPUT: patch friction velocity (m/s)
    u10_clm        => frictionvel_inst%u10_clm_patch         , &  ! OUTPUT: patch 10-m wind (m/s)
    fsa            => solarabs_inst%fsa_patch                , &  ! OUTPUT: patch solar radiation absorbed (total) (W/m2)
    albd           => surfalb_inst%albd_patch                , &  ! OUTPUT: patch surface albedo (direct)
    albi           => surfalb_inst%albi_patch                , &  ! OUTPUT: patch surface albedo (diffuse)
    t_ref2m        => temperature_inst%t_ref2m_patch         , &  ! OUTPUT: patch 2 m height surface air temperature (K)
    qflx_evap_tot  => waterflux_inst%qflx_evap_tot_patch     , &  ! OUTPUT: patch total evapotranspiration flux (kg H2O/m2/s)
    q_ref2m        => waterstate_inst%q_ref2m_patch          , &  ! OUTPUT: patch 2 m height surface specific humidity (kg/kg)

                                                                  ! *** Multilayer canopy variables ***
    zref           => mlcanopy_inst%zref_forcing             , &  ! Atmospheric reference height (m)
    tref           => mlcanopy_inst%tref_forcing             , &  ! Air temperature at reference height (K)
    thref          => mlcanopy_inst%thref_forcing            , &  ! Atmospheric potential temperature at reference height (K)
    thvref         => mlcanopy_inst%thvref_forcing           , &  ! Atmospheric virtual potential temperature at reference height (K)
    qref           => mlcanopy_inst%qref_forcing             , &  ! Specific humidity at reference height (kg/kg)
    eref           => mlcanopy_inst%eref_forcing             , &  ! Vapor pressure at reference height (Pa)
    uref           => mlcanopy_inst%uref_forcing             , &  ! Wind speed at reference height (m/s)
    pref           => mlcanopy_inst%pref_forcing             , &  ! Air pressure at reference height (Pa)
    co2ref         => mlcanopy_inst%co2ref_forcing           , &  ! Atmospheric CO2 at reference height (umol/mol)
    o2ref          => mlcanopy_inst%o2ref_forcing            , &  ! Atmospheric O2 at reference height (mmol/mol)
    rhoair         => mlcanopy_inst%rhoair_forcing           , &  ! Air density at reference height (kg/m3)
    rhomol         => mlcanopy_inst%rhomol_forcing           , &  ! Molar density at reference height (mol/m3)
    mmair          => mlcanopy_inst%mmair_forcing            , &  ! Molecular mass of air at reference height (kg/mol)
    cpair          => mlcanopy_inst%cpair_forcing            , &  ! Specific heat of air (constant pressure) at reference height (J/mol/K)
    solar_zen      => mlcanopy_inst%solar_zen_forcing        , &  ! Solar zenith angle (radians)
    swskyb         => mlcanopy_inst%swskyb_forcing           , &  ! Atmospheric direct beam solar radiation (W/m2)
    swskyd         => mlcanopy_inst%swskyd_forcing           , &  ! Atmospheric diffuse solar radiation (W/m2)
    lwsky          => mlcanopy_inst%lwsky_forcing            , &  ! Atmospheric longwave radiation (W/m2)
    qflx_rain      => mlcanopy_inst%qflx_rain_forcing        , &  ! Rainfall (mm H2O/s = kg H2O/m2/s)
    qflx_snow      => mlcanopy_inst%qflx_snow_forcing        , &  ! Snowfall (mm H2O/s = kg H2O/m2/s)
    tacclim        => mlcanopy_inst%tacclim_forcing          , &  ! Average air temperature for acclimation (K)
    ncan           => mlcanopy_inst%ncan_canopy              , &  ! Number of aboveground layers
    lai            => mlcanopy_inst%lai_canopy               , &  ! Leaf area index of canopy (m2/m2)
    sai            => mlcanopy_inst%sai_canopy               , &  ! Stem area index of canopy (m2/m2)
    swveg          => mlcanopy_inst%swveg_canopy             , &  ! Absorbed solar radiation: vegetation (W/m2)
    lwup           => mlcanopy_inst%lwup_canopy              , &  ! Upward longwave radiation above canopy (W/m2)
    shflx          => mlcanopy_inst%shflx_canopy             , &  ! Total sensible heat flux, including soil (W/m2)
    lhflx          => mlcanopy_inst%lhflx_canopy             , &  ! Total latent heat flux, including soil (W/m2)
    etflx          => mlcanopy_inst%etflx_canopy             , &  ! Total water vapor flux, including soil (mol H2O/m2/s)
    ustar          => mlcanopy_inst%ustar_canopy             , &  ! Friction velocity (m/s)
    albsoib        => mlcanopy_inst%albsoib_soil             , &  ! Direct beam albedo of ground (-)
    albsoid        => mlcanopy_inst%albsoid_soil             , &  ! Diffuse albedo of ground (-)
    swsoi          => mlcanopy_inst%swsoi_soil               , &  ! Absorbed solar radiation: ground (W/m2)
    lwsoi          => mlcanopy_inst%lwsoi_soil               , &  ! Absorbed longwave radiation: ground (W/m2)
    rnsoi          => mlcanopy_inst%rnsoi_soil               , &  ! Net radiation: ground (W/m2)
    tg             => mlcanopy_inst%tg_soil                  , &  ! Soil surface temperature (K)
    tg_bef         => mlcanopy_inst%tg_bef_soil              , &  ! Soil surface temperature for previous timestep (K)
    rhg            => mlcanopy_inst%rhg_soil                 , &  ! Relative humidity of airspace at soil surface (fraction)
    soilres        => mlcanopy_inst%soilres_soil             , &  ! Soil evaporative resistance (s/m)
    soil_t         => mlcanopy_inst%soil_t_soil              , &  ! Temperature of first snow/soil layer (K)
    soil_dz        => mlcanopy_inst%soil_dz_soil             , &  ! Depth to temperature of first snow/soil layer (m)
    soil_tk        => mlcanopy_inst%soil_tk_soil             , &  ! Thermal conductivity of first snow/soil layer (W/m/K)
    dlai           => mlcanopy_inst%dlai_profile             , &  ! Canopy layer leaf area index (m2/m2)
    dsai           => mlcanopy_inst%dsai_profile             , &  ! Canopy layer stem area index (m2/m2)
    dpai           => mlcanopy_inst%dpai_profile             , &  ! Canopy layer plant area index (m2/m2)
    dlai_frac      => mlcanopy_inst%dlai_frac_profile        , &  ! Canopy layer leaf area index (fraction of canopy total)
    dsai_frac      => mlcanopy_inst%dsai_frac_profile        , &  ! Canopy layer stem area index (fraction of canopy total)
    fracsun        => mlcanopy_inst%fracsun_profile          , &  ! Canopy layer sunlit fraction (-)
    tair           => mlcanopy_inst%tair_profile             , &  ! Canopy layer air temperature (K)
    eair           => mlcanopy_inst%eair_profile             , &  ! Canopy layer vapor pressure (Pa)
    cair           => mlcanopy_inst%cair_profile             , &  ! Canopy layer atmospheric CO2 (umol/mol)
    tair_bef       => mlcanopy_inst%tair_bef_profile         , &  ! Canopy layer air temperature for previous timestep (K)
    eair_bef       => mlcanopy_inst%eair_bef_profile         , &  ! Canopy layer vapor pressure for previous timestep (Pa)
    cair_bef       => mlcanopy_inst%cair_bef_profile         , &  ! Canopy layer atmospheric CO2 for previous timestep (umol/mol)
    swleaf         => mlcanopy_inst%swleaf_leaf              , &  ! Leaf absorbed solar radiation (W/m2 leaf)
    lwleaf         => mlcanopy_inst%lwleaf_leaf              , &  ! Leaf absorbed longwave radiation (W/m2 leaf)
    rnleaf         => mlcanopy_inst%rnleaf_leaf              , &  ! Leaf net radiation (W/m2 leaf)
    tleaf          => mlcanopy_inst%tleaf_leaf               , &  ! Leaf temperature (K)
    tleaf_bef      => mlcanopy_inst%tleaf_bef_leaf           , &  ! Leaf temperature for previous timestep (K)
    tleaf_hist     => mlcanopy_inst%tleaf_hist_leaf          , &  ! Leaf temperature (not sun/shade average) for history files (K)
    lwp            => mlcanopy_inst%lwp_leaf                 , &  ! Leaf water potential (MPa)
    lwp_hist       => mlcanopy_inst%lwp_hist_leaf              &  ! Leaf water potential (not sun/shade average) for history files (MPa)
    )

    ! Get current step counter (nstep) and step size (dtime)

    nstep = get_nstep()
    dtime = get_step_size()

    ! Set number of sub-time steps for flux calculations

    num_sub_steps = int(dtime / dtime_substep)

    ! Build filter of patches to process with multilayer canopy

    mlcanopy_inst%filtern = 0
    mlcanopy_inst%filter(:) = ispval
    
    do fp = 1, num_exposedvegp
       p = filter_exposedvegp(fp)
       g = patch%gridcell(fp)
!      if (grc%latdeg(g) .gt. -2.9_r8 .and. grc%latdeg(g) .lt. -2.7_r8) then
!      if (grc%londeg(g) .gt. 294.5_r8 .and. grc%londeg(g) .lt. 295.5_r8) then
       mlcanopy_inst%filtern = mlcanopy_inst%filtern + 1
       mlcanopy_inst%filter(mlcanopy_inst%active_n) = p
!      end if
!      end if
    end do

    ! Initialize canopy vertical structure and profiles. This is only done
    ! once (on first time step) because the forcing height (and therefore
    ! the volume of air in the surface layer) changes between time steps.
    ! The call for initialization is triggered by zref = spval.

    ml_vert_init = 0
    do fp = 1, ml%filtern
       p = mlcanopy_inst%filter(fp)
       if (zref(p) == spval) ml_vert_init = 1
    end do

    if (ml_vert_init == 1) then
       if (masterproc) then
          write (iulog,*) 'Attempting to initialize multilayer canopy vertical structure .....'
       end if
       
       call initVerticalStructure (bounds, mlcanopy_inst%filtern, mlcanopy%filter, &
       canopystate_inst, frictionvel_inst, mlcanopy_inst)

       call initVerticalProfiles (mlcanopy_inst%filtern, mlcanopy%filter, &
       atm2lnd_inst, mlcanopy_inst)

       if (masterproc) then
          write (iulog,*) 'Successfully initialized multilayer canopy vertical structure'
       end if
    end if

    ! Copy CLM variables to multilayer canopy variables. Note the
    ! distinction between grid cell (g), column (c), and patch (p)
    ! variables. All multilayer canopy variables are for patches.

    do fp = 1, mlcanopy_inst%filtern
       p = mlcanopy%filter(fp)
       c = patch%column(p)
       g = patch%gridcell(p)

       ! Atmospheric forcing: CLM grid cell (g) variables -> patch (p) variables
       ! RGK: THESE 5 are ML VARS
       uref(p) = max (wind_forc_min, sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))
       swskyb(p,ivis) = forc_solad(g,ivis) 
       swskyd(p,ivis) = forc_solai(g,ivis)
       swskyb(p,inir) = forc_solad(g,inir)
       swskyd(p,inir) = forc_solai(g,inir)

       ! Re-partition direct and diffuse radiation if desired

       if (fracdir >= 0._r8) then
          totrad = swskyb(p,ivis) + swskyd(p,ivis)
          swskyb(p,ivis) = totrad * fracdir
          swskyd(p,ivis) = totrad * (1._r8 - fracdir)
          totrad = swskyb(p,inir) + swskyd(p,inir)
          swskyb(p,inir) = totrad * fracdir
          swskyd(p,inir) = totrad * (1._r8 - fracdir)
       end if

       ! Atmospheric forcing: CLM column (c) variables -> patch (p) variables

       tref(p) = forc_t(c)
       qref(p) = forc_q(c)
       pref(p) = forc_pbot(c)
       lwsky(p) = forc_lwrad(c)
       qflx_rain(p) = forc_rain(c)
       qflx_snow(p) = forc_snow(c)

       ! CO2 and O2: CLM grid cell (g) -> patch (p). Note that the units
       ! conversion requires pbot.

       co2ref(p) = forc_pco2(g) / forc_pbot(c) * 1.e06_r8  ! Pa -> umol/mol
       o2ref(p)  = forc_po2(g) / forc_pbot(c) * 1.e03_r8   ! Pa -> mmol/mol

       ! Miscellaneous

       tacclim(p) = t_a10_patch(p)

       ! Ground (soil) albedos: CLM column (c) -> patch (p)

       albsoib(p,ivis) = albgrd(c,ivis) ; albsoib(p,inir) = albgrd(c,inir)
       albsoid(p,ivis) = albgri(c,ivis) ; albsoid(p,inir) = albgri(c,inir)

       ! Soil evaporative resistance: CLM column (c) -> patch (p)

       soilres(p) = soilresis(c)

       ! Properties of first soil layer: CLM column (c) -> patch (p)

       soil_t(p) = t_soisno(c,snl(c)+1)          ! Temperature of first snow/soil layer (K)
       soil_dz(p) = (z(c,snl(c)+1)-zi(c,snl(c))) ! Depth to temperature of first snow/soil layer (m)
       soil_tk(p) = thk(c,snl(c)+1)              ! Thermal conductivity of first snow/soil layer (W/m/K)

    end do

    ! Solar zenith angle. Need to subtract one time step (-dtime) because
    ! zenith angle is calculated for the beginning of the time step. So use
    ! calendar day at beginning of the time step (nstep-1).

    caldaym1 = get_curr_calday(offset=-int(dtime))
    call shr_orb_decl (caldaym1, eccen, mvelpp, lambm0, obliqr, declinm1, eccf)

    do fp = 1, mlcanopy_inst%filtern
       p = mlcanopy%filter(fp)
       c = patch%column(p)
       g = patch%gridcell(p)
       lat = grc%latdeg(g) * pi / 180._r8
       lon = grc%londeg(g) * pi / 180._r8
       coszen = shr_orb_cosz (caldaym1, lat, lon, declinm1)
       solar_zen(p) = acos(max(0.01_r8,coszen))

       ! Compare coszen to that expected from CLM

       if (abs(coszen-surfalb_inst%coszen_col(c)) .gt. 1.e-03_r8) then
          write (iulog,*) nstep, coszen, surfalb_inst%coszen_col(c)
          call endrun (msg=' ERROR: MLCanopyFluxes: coszen error')
       end if
    end do

    ! Derived atmospheric input

    do fp = 1, mlcanopy_inst%filtern
       p = mlcanopy%filter(fp)
       eref(p) = qref(p) * pref(p) / (mmh2o / mmdry + (1._r8 - mmh2o / mmdry) * qref(p))
       rhomol(p) = pref(p) / (rgas * tref(p))
       rhoair(p) = rhomol(p) * mmdry * (1._r8 - (1._r8 - mmh2o/mmdry) * eref(p) / pref(p))
       mmair(p) = rhoair(p) / rhomol(p)
       cpair(p) = cpd * (1._r8 + (cpw/cpd - 1._r8) * qref(p)) * mmair(p)
       thref(p) = tref(p) + lapse_rate * zref(p)
       thvref(p) = thref(p) * (1._r8 + 0.61_r8 * qref(p))
    end do

    ! Update leaf and stem area profile for current values

    do fp = 1, mlcanopy_inst%filtern
       p = mlcanopy%filter(fp)

       ! Get values for current time step from CLM

       lai(p) = elai(p)
       sai(p) = esai(p)

       ! Vertical profiles

       do ic = 1, ncan(p)
          dlai(p,ic) = dlai_frac(p,ic) * lai(p)
          dsai(p,ic) = dsai_frac(p,ic) * sai(p)
          dpai(p,ic) = dlai(p,ic) + dsai(p,ic)
       end do

       totpai = sum(dpai(p,1:ncan(p)))
       if (abs(totpai - (lai(p)+sai(p))) > 1.e-06_r8) then
          call endrun (msg=' ERROR: MLCanopyFluxes: plant area index not updated correctly')
       end if

    end do

    ! Solar radiation transfer through the canopy

    call SolarRadiation (mlcanopy_inst,patch%itype(bounds%begp:bounds%endp))

    ! Plant hydraulics

    call SoilResistance (mlcanopy_inst%filtern, mlcanopy%filter, &
    soilstate_inst, waterstate_inst, mlcanopy_inst)

    call PlantResistance (mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

    ! Canopy profile of photosynthetic capacity

    call CanopyNitrogenProfile (mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

    ! Use sub-stepping to calculate fluxes over the full time step

    do niter = 1, num_sub_steps

       ! Save values for previous timestep

       do fp = 1, mlcanopy_inst%filtern
          p = mlcanopy%filter(fp)
          tg_bef(p) = tg(p)
          do ic = 1, ncan(p)
             tleaf_bef(p,ic,isun) = tleaf(p,ic,isun)
             tleaf_bef(p,ic,isha) = tleaf(p,ic,isha)
             tair_bef(p,ic) = tair(p,ic)
             eair_bef(p,ic) = eair(p,ic)
             cair_bef(p,ic) = cair(p,ic)
          end do
       end do

       ! Canopy interception

       call CanopyInterception (mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

       ! Longwave radiation transfer through the canopy

       call LongwaveRadiation (bounds, mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

       ! Net radiation at each layer and at ground

       do fp = 1, mlcanopy_inst%filtern
          p = mlcanopy%filter(fp)
          do ic = 1, ncan(p)
             rnleaf(p,ic,isun) = swleaf(p,ic,isun,ivis) + swleaf(p,ic,isun,inir) + lwleaf(p,ic,isun)
             rnleaf(p,ic,isha) = swleaf(p,ic,isha,ivis) + swleaf(p,ic,isha,inir) + lwleaf(p,ic,isha)
          end do
          rnsoi(p) = swsoi(p,ivis) + swsoi(p,inir) + lwsoi(p)
       end do

       ! Leaf heat capacity

       call LeafHeatCapacity (mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

       ! Leaf boundary layer conductance

       call LeafBoundaryLayer (mlcanopy_inst%filtern, mlcanopy%filter, isun, mlcanopy_inst)
       call LeafBoundaryLayer (mlcanopy_inst%filtern, mlcanopy%filter, isha, mlcanopy_inst)

       ! Photosynthesis and stomatal conductance

       call LeafPhotosynthesis (mlcanopy_inst%filtern, mlcanopy%filter, isun, mlcanopy_inst)
       call LeafPhotosynthesis (mlcanopy_inst%filtern, mlcanopy%filter, isha, mlcanopy_inst)

       ! Relative humidity in soil airspace

       do fp = 1, mlcanopy_inst%filtern
          p = mlcanopy%filter(fp)
          c = patch%column(p)
          rhg(p) = exp(grav * mmh2o * smp_l(c,1)*1.e-03_r8 / (rgas * t_soisno(c,1)))
       end do

       ! Canopy turbulence, scalar source/sink fluxes for leaves and soil, and
       ! scalar profiles above and within the canopy

       call CanopyTurbulence (niter, mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

       ! Update leaf water potential for the current transpiration rate

       call LeafWaterPotential (mlcanopy_inst%filtern, mlcanopy%filter, isun, mlcanopy_inst)
       call LeafWaterPotential (mlcanopy_inst%filtern, mlcanopy%filter, isha, mlcanopy_inst)

       ! Update canopy intercepted water for evaporation and dew

       call CanopyEvaporation (mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

       ! Fluxes need to be accumulated over all sub-time steps. Other
       ! variables are instantaneous for the final sub-time step.

       call SubTimeStepFluxIntegration (niter, num_sub_steps, mlcanopy_inst%filtern, mlcanopy%filter, &
       flux_accumulator, flux_accumulator_profile, flux_accumulator_leaf, mlcanopy_inst)

    end do    ! End sub-stepping loop

    ! Sum leaf and soil fluxes and other canopy diagnostics

    call CanopyFluxesDiagnostics (mlcanopy_inst%filtern, mlcanopy%filter, mlcanopy_inst)

    ! Leaf temperature and leaf water potential are prognostic variables
    ! for sunlit and shaded leaves. But sun/shade fractions change over
    ! time. Merge temperature and leaf water potential for sunlit and
    ! shaded leaves to layer-average value, which is used at next time step.
    ! Use this method (rather than retaining sun/shade states) because some
    ! shade leaf becomes sun leaf (and vice versa) between time steps as fsun
    ! changes.

    do fp = 1, mlcanopy_inst%filtern
       p = mlcanopy%filter(fp)
       do ic = 1, ncan(p)

          ! First save sun/shade leaves for model output

          tleaf_hist(p,ic,isun) = tleaf(p,ic,isun)
          tleaf_hist(p,ic,isha) = tleaf(p,ic,isha)
          lwp_hist(p,ic,isun) = lwp(p,ic,isun)
          lwp_hist(p,ic,isha) = lwp(p,ic,isha)

          ! Now merge sun/shade leaves

          if (dpai(p,ic) > 0._r8) then
             tleaf(p,ic,isun) = tleaf(p,ic,isun) * fracsun(p,ic) + tleaf(p,ic,isha) * (1._r8 - fracsun(p,ic))
             tleaf(p,ic,isha) = tleaf(p,ic,isun)
             lwp(p,ic,isun) = lwp(p,ic,isun) * fracsun(p,ic) + lwp(p,ic,isha) * (1._r8 - fracsun(p,ic))
             lwp(p,ic,isha) = lwp(p,ic,isun)
          end if

       end do
    end do

    ! Copy multilayer canopy variables to CLM variables. These are
    ! passed from CLM to CAM. The variables passed to CAM are found
    ! in the CLM routine: main/lnd2atmType.F90. These are CLM grid
    ! cell variables. The mapping from CLM patch to CLM gridcell is
    ! done in: main/lnd2atmMod.F90

    if (mlcan_to_clm == 1) then
       do fp = 1, mlcanopy_inst%filtern
          p = mlcanopy%filter(fp)
          albd(p,ivis) = 0._r8 ; albd(p,inir) = 0._r8
          albi(p,ivis) = 0._r8 ; albi(p,inir) = 0._r8
          taux(p) = 0._r8
          tauy(p) = 0._r8
          eflx_lh_tot(p) = lhflx(p)
          eflx_sh_tot(p) = shflx(p)
          eflx_lwrad_out(p) = lwup(p)
          qflx_evap_tot(p) = etflx(p) * mmh2o
          fv(p) = ustar(p)
          u10_clm(p) = 0._r8
          t_ref2m(p) = 0._r8
          q_ref2m(p) = 0._r8
          fsa(p) = swveg(p,ivis) + swveg(p,inir) + swsoi(p,ivis) + swsoi(p,inir)
       end do
    end if

    end associate
  end subroutine MLCLMCanopyFluxes
  
end module MLCLMCanopyCouplerMod
