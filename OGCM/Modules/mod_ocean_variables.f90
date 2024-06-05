!input physical parameters of the run
module key_switches
implicit none

integer ksw_ts,  ksw_uv,              &     !TS (left) and UV (right) equation solving (0-no, 1-yes)                                                                                                      
        ksw_age, ksw_pt,              &     !Ideal age (left) and PassTracer (right) equation solving (0-no, 1-yes)                                                                                                                       
        ksw_lat, ksw_lat4,            &     !Lateral 2nd (left, 0-no, 1-sigma, 2-z, 3-rho (2&3 for TS)) and 4th (right, 0-no, 1-yes, UV only) order mix
        ksw_vert,                     &     !vertical mix (0-constant coeff, 1-PP, 2-MY)                                                                                      
        ksw_dens,                     &     !pressure gradient computation (0-no (constant density), 1-local potential density)
        ksw_ith, ksw_itr, ksw_idn,    &     !sea ice thermodynamics, transport and dynamics using (0-no, 1-yes)                                                                                                                      
        ksw_ssbc,                     &     !Type of surface boundary conditions (1-surface T&S and wind stress are prescribed; 2-T&S fluxes and wind stress are prescribed; 3-T&S fluxes and wind stress are simulated
        ksw_wflux,                    &     !normalize global mean salt balance (0-no, 1-normalize water/salt flux (depending on variable volume))                                                               
        ksw_lbc_ts, ksw_lbc_uv, ksw_lbc_ssh !open boundary conditions for T&S, U&V and SSH using (0-no, 1-yes)                                                                                                         

real(4) sst_relax,    sss_relax,      &     !Surface nudging rate for temperature and salinity [m/s]            
        ldiff_ts_ref, ldiff_ts_upw, ldiff_ts_smg,   &     !lateral diffusivity ratio (undim) for temperature: background, upwind and smagorinsky                               
        lvisc_ref,    lvisc_upw,    lvisc_smg,      &     !lateral vicosity(2nd order) ratio (undim): background, upwind and smagorinsky
        lvisc_4_ref,  lvisc_4_upw,  lvisc_4_smg,    &     !lateral vicosity(4th order) ratio (undim): background, upwind and smagorinsky
        tsfrac_lat,                   &     !fraction of salinity lateral diffusivity due to one for temperature
        vdiff_ts_min, vdiff_ts_max,   &     !vertical background and top diffusivity for T [m**2/s]
        vvisc_min,    vvisc_max,      &     !vertical background and top viscosity for UV [m**2/s]        
        tsfrac_vert,                  &     !fraction of salinity lateral diffusion due to one for temperature
        rfric_min, rfric_max                !minimum and maximum rayleigh dissipation rate (1/s)

endmodule key_switches

!-------------module for description common ogcm variables and task control parameters---------
module ocean_variables
use mpi_parallel_tools
implicit none

save

!barotropic dynamics arrays
real(8),allocatable:: ssh(:,:),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
                    ssh_n(:,:),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
                    ubrtr(:,:),     &  !barotropic velocity      zonal[m/s] at current time step (internal mode)
                    vbrtr(:,:),     &  !barotropic velocity meridional[m/s] at current time step (internal mode)
                   RHSx2d(:,:),     &  !x-component of external force(barotropic)
                   RHSy2d(:,:)         !y-component of external force(barotropic)

real(4), allocatable:: sls(:,:)

!3d dynamics arrays
real(4),allocatable::  uu(:,:,:),   &  !     zonal velocity [m/s]
                       vv(:,:,:),   &  !meridional velocity [m/s]
                       uh(:,:,:),   &  !     zonal velocity by depth [m*m/s]
                       vh(:,:,:),   &  !meridional velocity by depth [m*m/s]
                       ww(:,:,:),   &  !  vertical velocity in sigma-coord [m/s]
                       tt(:,:,:),   &  !potential temperature[°C]
                       ss(:,:,:),   &  !salinity [psu]
                  den_sgt(:,:,:),   &  !in situ density sigma-t (p= RefDen*gravity*Z ) (kg/m^3)
                  den_sg0(:,:,:),   &  !potential density sigma-0 (p=0) (kg/m^3)
                  den_sgm(:,:,:),   &  !potential density sigma-zmean (for p=hmp2 dbar, 2000.0 for Global Ocean), see basinpar.fi
                  den_ini(:,:,:),   &  !Initial potential density (for steric level computation)
                  yng_slp(:,:,:)       !Joung module for sea level pressure (Pa)

real(4),allocatable::  q2(:,:,:),   &  !turbulent kinetic energy (m/s)^2
                      q2l(:,:,:)       !turbulent kinetic energy by turbulent length scale (m/s)^2*m

real(4), allocatable::  rhs_q2(:,:,:),  &
                       rhs_q2l(:,:,:),  &
                       coef_q2(:,:,:),  &
                      coef_q2l(:,:,:)

real(4),allocatable:: stress_t(:,:,:),    &     !Horizontal tension tensor component
                      stress_s(:,:,:)           !Horizontal shearing tensor component

!surface and bottom boundary condition types    ! igrz[T,S]= 1 : f = f0 
integer,allocatable:: igrzts(:,:)               !          = 2 : df/dz = f0

!Passive tracer arrays
real(4),allocatable:: pass_tracer(:,:,:),   &     !Passive tracer
                      pt_forc_adv(:,:),     &     !Passive tracer surface advective forcing
                      pt_forc_dif(:,:)            !Passive tracer surface diffusive forcing

integer,allocatable:: igrzpt(:,:)                 !Passive tracer boundary condition type

real(4),allocatable:: age(:,:,:),   &     !Water ideal age [years]
               age_forc_adv(:,:),   &     !Ideal age surface advective forcing
               age_forc_dif(:,:)          !Ideal age surface diffusive forcing

integer,allocatable:: igrzage(:,:)                !Passive tracer boundary condition type


! coefficients of viscosity  and diffusion
! lateral viscosity and diffusion coefficients for t, s, u, v:
! horizontal:
real(4),allocatable::  amts_ref(:,:,:),     &   !T lateral diffusion in T-points [m^2/s]
                       amuv_ref(:,:,:),     &   !U and V lateral viscosity in T-points [m^2/s]
                       amuv4_ref(:,:,:)         !U and V 4-th order lateral viscosity in T-points [m^4/s]^(1/2)

real(4),allocatable::  amts(:,:,:),     &   !T lateral diffusion in T-points [m^2/s]
                       amuv(:,:,:),     &   !U and V lateral viscosity in T-points [m^2/s]
                       amuv4(:,:,:)         !U and V 4-th order lateral viscosity in T-points [m^4/s]^(1/2)

! vertical viscous and diffusion functions
real(4),allocatable:: anzt(:,:,:),     &     !T vertical diffusion [m^2/s]
                      anzu(:,:,:),     &     !U and V vertical viscosity [m^2/s]
                      vbf2(:,:,:),     &     !Brunt–Vaisala  frequency square [1/s^2]
                      shf2(:,:,:)            !Vertical shear frequency square [1/s^2]

real(4),allocatable:: mld_dens(:,:)          !mixed layer depth [m] based on potential density

! sea surface boundary condition
real(4), allocatable:: tflux_adv(:,:),   &       !total advective surface temperature flux [°C*m/s]
                       tflux_dif(:,:),   &       !total diffusive surface temperature flux [°C*m/s]
                      swflux(:,:),       &       !total surface shortwave radiation flux [°C*m/s]
                       sflux_adv(:,:),   &       !total advective surface salinity flux [psu*m/s]
                       sflux_dif(:,:),   &       !total diffusive surface salinity flux [psu*m/s]
                    vol_flux(:,:),       &       !total volume flux [m/s] (precip+evap+ice) [m/s]
                    bot_fric(:,:),       &       !bottom fricrion bulk coefficient [m/s]
                       cd_ao(:,:),       &       !drag coefficient air-ocean [kg/m^2/s]
                       cd_ai(:,:),       &       !drag coefficient air-ice [kg/m^2/s]
                       cd_io(:,:),       &       !drag coefficient ice-ocean [kg/m^2/s]
                  divswrad(:,:,:),       &       !shortwave radiation divergence coefficients
                        dkft(:,:),       &       !relaxation coefficient for SST, [m/s]
                        dkfs(:,:),       &       !relaxation coefficient for SSS, [m/s]
                    sensheat(:,:),       &       !sensible heat flux
                     latheat(:,:),       &       !latent heat flux
                      lw_bal(:,:),       &       !longwave radiation balance
                  sw_bal_atm(:,:),       &       !shortwave radiation balance from atmosphere
                   sw_bal_oc(:,:),       &       !shortwave radiation balance to ocean
                  hf_tot_atm(:,:),       &       !total heat flux from air
                   hf_tot_oc(:,:),       &       !total heat flux into ocean
                   sf_tot_oc(:,:),       &       !total salt flux into ocean
                  wf_tot_atm(:,:),       &       !total water flux from air
                   wf_tot_oc(:,:)                !total water flux into ocean

!Atmospheric arrays for bulk-formulae
real(4),allocatable:: tatm(:,:),   &    !Air temperature, [°C]
                      qatm(:,:),   &    !Air humidity, [kg/kg]
                      rain(:,:),   &    !rain, [kg/m^2/s]
                      snow(:,:),   &    !snow, [kg/m^2/s]
                       lwr(:,:),   &    !Downward  longwave radiation, [W/m^2]
                       swr(:,:),   &    !Downward shortwave radiation, [W/m^2]
                      slpr(:,:),   &    !Sea level pressure, [Pa]
                      uwnd(:,:),   &    !Zonal      wind speed, [m/s]
                      vwnd(:,:),   &    !Meridional wind speed, [m/s]
                      taux(:,:),   &    !Zonal      wind stress, [Pa]
                      tauy(:,:)         !Meridional wind stress, [Pa]

real(4), allocatable::  r_diss(:,:),  & !Rayleigh friction scale (1/s)
                    AlbOpWater(:,:)     !Open water albedo

! array description for ice.
real(4), allocatable:: hice(:,:,:),    &       !Ice volume, [m]
                       aice(:,:,:),    &       !Ice compactness, [0-1]
                      hsnow(:,:,:),    &       !Snow volume, [m]
                       sice(:,:,:),    &       !Ice salinity, [PSU*m]
                       mistot(:,:),    &       !total ice and snow mass
                       aistot(:,:),    &       !Total ice compactness, [0-1]
                        hitot(:,:),    &       !Total ice volume, m
                        hstot(:,:),    &       !Total snow volume, m
                         pice(:,:),    &       !Ice pressure, [Pa*m]
                        sig_d(:,:),    &       !Stress tensor divergence component
                        sig_t(:,:),    &       !Stress tensor tension component
                        sig_s(:,:),    &       !Stress tensor shear component
                        eps_d(:,:),    &       !Strain rate tensor divergence component
                        eps_t(:,:),    &       !Strain rate tensor tension component
                        eps_s(:,:),    &       !Strain rate tensor shear component
                         uice(:,:),    &       !Ice zonal velocity, [m/s]
                         vice(:,:)             !Ice meridional velocity, [m/s]

real(4), allocatable::  hf_sugar(:,:)          !heat  flux due to new ice formation

real(4), allocatable:: idyn_coefx(:,:),      &
                       idyn_coefy(:,:),      &
                   idyn_coefx_tot(:,:),      &
                   idyn_coefy_tot(:,:),      &
                        idyn_rhsx(:,:),      &
                        idyn_rhsy(:,:),      &
                    idyn_rhsx_tot(:,:),      &
                    idyn_rhsy_tot(:,:),      &
                          mistotu(:,:),      &
                          mistotv(:,:),      &           
                            piceh(:,:)

real(4), allocatable::  taux_atm(:,:),       & !Zonal      stress from air, [Pa]
                        tauy_atm(:,:),       & !Meridional stress from air, [Pa]
                        taux_oc(:,:),        & !Zonal      stress to water, [Pa]
                        tauy_oc(:,:),        & !Meridional stress to water, [Pa]
                        tx_adv(:,:),         & !Total advective surface zonal momentum flux [m/s*m/s]
                        tx_dif(:,:),         & !Total diffusive surface zonal momentum flux [m/s*m/s]
                        ty_adv(:,:),         & !Total advective surface meridional momentum flux [m/s*m/s]
                        ty_dif(:,:),         & !Total diffusive surface meridional momentum flux [m/s*m/s]
                        tx_coef(:,:),        & !Total diffusive surface zonal momentum cofficient [m/s]
                        ty_coef(:,:)           !Total diffusive surface meridional momentum cofficient [m/s]

real(4), allocatable::  wind_ao(:,:),   &    !Air-ocean speed shear, [m/s]
                        wind_ai(:,:),   &    !Air-ice speed shear, [m/s]
                        ioshear(:,:),   &    !ice-ocean speed shear, [m/s]
                        windx(:,:),     &    !Zonal      wind speed in u-points, [m/s]
                        windy(:,:)           !Meridional wind speed in v-points, [m/s]

real(4),allocatable:: tt_calc(:,:,:),     &
                      ss_calc(:,:,:),     &
                      uu_calc(:,:,:),     &
                      vv_calc(:,:,:),     &
                  ssdata_calc(:,:,:),     &
                  ssflux_calc(:,:,:),     &
                   vflux_calc(:,:,:),     &
                     hhq_calc(:,:),       &
                     hhu_calc(:,:),       &
                     hhv_calc(:,:),       &
                     txo_calc(:,:),       &
                     tyo_calc(:,:),       &
                    uwnd_calc(:,:),       &
                    vwnd_calc(:,:)

real(4),allocatable:: pt_calc(:,:,:), age_calc(:,:,:)
real(4),allocatable:: den_sgt_calc(:,:,:), den_sg0_calc(:,:,:), den_sgm_calc(:,:,:)

real(4), allocatable::  aice_calc(:,:),   &
                        hice_calc(:,:),   &
                       hsnow_calc(:,:),   &
                        uice_calc(:,:),   &
                        vice_calc(:,:)

real(8), allocatable:: ssh_calc(:,:),     &
                   ssh_hhq_calc(:,:)

real(4), allocatable:: sls_calc(:,:),     &
                  mld_dens_calc(:,:)

real(4), allocatable:: aux_array3d_tgr1(:,:,:),      &
                       aux_array3d_tgr2(:,:,:),      &
                       aux_array3d_tgr3(:,:,:),      &
                       aux_array3d_tgr4(:,:,:),      &
                       aux_array3d_tgr5(:,:,:),      &
                                    z3d(:,:,:),      &
                             den_isopyc(:,:,:),      &
                       aux_array3d_wgr1(:,:,:),      &
                       aux_array3d_wgr2(:,:,:),      &
                       aux_array3d_wgr3(:,:,:),      &
                       aux_array3d_wgr4(:,:,:),      &
                    aux_array3d_zgr_loc(:,:,:),      &
                   aux_array3d_zgr_glob(:,:,:)

real(4), allocatable:: aux_array_itr1(:,:,:),   &
                       aux_array_itr2(:,:,:),   &
                       aux_array_itr3(:,:,:),   &
                       aux_array_itr4(:,:,:)

real(4), allocatable:: ssd(:,:,:),        &
                       ssf(:,:,:),        &
                       vfl(:,:,:)

real(4), allocatable:: aux_array2d_01(:,:),     &
                       aux_array2d_02(:,:)

integer meancalc             !calculator for time mean output

endmodule ocean_variables
!-------------end module for description common ogcm variables and task control parameters---------
