module init_arrays_routes
implicit none

contains

!array boundary definition for non-mpi arrays
!=======================================================================
subroutine non_mpi_array_boundary_definition
use basin_grid

       nx_start=mmm
       nx_end  =mm
       ny_start =nnn
       ny_end  =nn
 
       bnd_x1=nx_start-2
       bnd_x2=nx_end  +2
       bnd_y1=ny_start-2
       bnd_y2=ny_end  +2

endsubroutine non_mpi_array_boundary_definition

!allocation of arrays
!=======================================================================
   subroutine model_grid_allocate
   use basin_grid

      allocate( lu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of t-grid
               lu1(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of t-grid (1 everywhere)
               luu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of h-grid (0 on boundary)
               luh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of h-grid (1 on boundary)
               lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of u-grid (0 on boundary)
               lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of v-grid (0 on boundary)
               llu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of u-grid (0 on boundary)
               llv(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )                !mask of v-grid (0 on boundary)

               lu=0.0; lu1=0.0; luu=0.0; luh=0.0; lcu=0.0; lcv=0.0; llu=0.0; llv=0.0 

      allocate (lbasins(nx,ny))       !integer masks of regional basins
    
               lbasins=0
      
      allocate(      hhh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points)
                   hhh_n(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points) at previous step
                hhq_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at rest state
                hhu_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at rest state
                hhv_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points) at rest state
                     hhq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points)
                   hhq_n(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at previous step
                     hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points)
                   hhu_n(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at previous step
                     hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points)
                   hhv_n(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points) at previous step
                     rlh_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !coriolis-1 parameter on edge (t-centers) points
                     rlh_c(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !coriolis-2 parameter on edge (t-centers) points
                    z(nz), zw(nz+1),    &  !vertical sigma-levels (t-points and w-points)
                hzt(nz+1), dz(nz),      &  !steps between t-levels and w-levels
             dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between   t-points (in radians or meters)
             dx (bnd_x1:bnd_x2,bnd_y1:bnd_y2),dy (bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between u,v-points (in radians or meters)
             dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between   h-points (in radians or meters)
             dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between v,u-points (in radians or meters)
                                                               xt(nx),yt(ny),        &  !horizontal t-grid            x- and y-coordinates (in degrees)
                                                               xu(nx),yv(ny)    )       !horizontal u-grid and v-grid x- and y-coordinates (in degrees)

               hhh=0.0; hhh_n=0.0; hhq_rest=0.0; hhu_rest=0.0; hhv_rest=0.0
               hhq=0.0; hhq_n=0.0
               hhu=0.0; hhu_n=0.0;  hhv=0.0; hhv_n=0.0 
               rlh_s=0.0; rlh_c=0.0
               z=0.0; zw=0.0; hzt=0.0; dz=0.0
               dxt=0.0; dyt=0.0; dx=0.0; dy=0.0
               dxh=0.0; dyh=0.0; dxb=0.0; dyb=0.0
               xt=0.0d0; yt=0.0d0; xu=0.0d0; yv=0.0d0

      allocate( sqt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !grid area in T-points
                squ(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !grid area in U-points
                sqv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !grid area in V-points
                sqh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)     )    !grid area in H-points
      
                sqt=0.0; squ=0.0; sqv=0.0; sqh=0.0
              
      allocate( geo_lon_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of T-points
                geo_lat_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of T-points
                geo_lon_u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of U-points
                geo_lat_u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of U-points
                geo_lon_v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of V-points
                geo_lat_v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of V-points
                geo_lon_h(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of H-points
                geo_lat_h(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of H-points
                rotvec_coeff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)    )

                geo_lon_t=0.0d0; geo_lat_t=0.0d0; geo_lon_u=0.0d0; geo_lat_u=0.0d0
                geo_lon_v=0.0d0; geo_lat_v=0.0d0; geo_lon_h=0.0d0; geo_lat_h=0.0d0
                rotvec_coeff=0.0d0

  endsubroutine model_grid_allocate

!deallocation of arrays
!=======================================================================
  subroutine model_grid_deallocate
  use basin_grid

           deallocate(rotvec_coeff,geo_lat_h,geo_lon_h,geo_lat_v,geo_lon_v,geo_lat_u,geo_lon_u,geo_lat_t,geo_lon_t)
           deallocate(sqh,sqv,squ,sqt)
           deallocate(yv,xu,yt,xt,dyb,dxb,dyh,dxh,dy,dx,dyt,dxt,dz,hzt,zw,z,rlh_c,rlh_s)
           deallocate(hhv_n,hhv,hhu_n,hhu,hhq_n,hhq,hhv_rest,hhu_rest,hhq_rest,hhh_n,hhh)
           deallocate(lbasins)
           deallocate(llv,llu,lcv,lcu,luh,luu,lu1,lu)
  
  endsubroutine model_grid_deallocate

!allocation of arrays   
!=======================================================================
   subroutine ocean_variables_allocate
   use ocean_variables
   use basin_grid
   include 'locout.fi'
   include 'globout.fi'

      allocate( ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
              ssh_n(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
              ubrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !barotropic velocity      zonal[m/s] at current  time step [m] (internal mode)
              vbrtr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !barotropic velocity meridional[m/s] at current  time step [m] (internal mode)
             RHSx2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !x-component of external force(barotropic)
             RHSy2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      )  !y-component of external force(barotropic)
     
               ssh =0.0d0; ssh_n =0.0d0
               ubrtr=0.0d0;  vbrtr=0.0d0
               RHSx2d=0.0d0; RHSy2d=0.0d0
      
      allocate(sls(bnd_x1:bnd_x2,bnd_y1:bnd_y2))
               sls = 0.0

      allocate(uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !     zonal velocity [m/s]
               vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !meridional velocity [m/s]
               uh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !     zonal velocity by depth [m*m/s]
               vh(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !meridional velocity by depth [m*m/s]
               ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &  !  vertical velocity in sigma-coord [m/s]
               tt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !potential temperature[°C]
               ss(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !salinity [psu-35ppt]
          den_sgt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !in situ density sigma-t (p= RefDen*gravity*Z ) (kg/m^3)
          den_sg0(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !potential density sigma-0 (p=0) (kg/m^3)
          den_sgm(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !potential density sigma-zmean (for p=hmp2 dbar, 2000.0 for Global Ocean), see basinpar.fi 
          den_ini(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !Initial potential density (for steric level computation)
          yng_slp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz))       !Joung module for sea level pressure (Pa)            

               uu=0.0; vv=0.0; uh=0.0; vh=0.0; ww=0.0
               tt=0.0; ss=0.0
               den_sgt=0.0; den_sg0=0.0; den_sgm=0.0; den_ini=0.0; yng_slp=0.0

     allocate( q2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &  !turbulent kinetic energy (m/s)^2
              q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)    )  !turbulent kinetic energy by turbulent length scale (m/s)^2*m
     
               q2=0.0; q2l=0.0

     allocate( rhs_q2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
              rhs_q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
              coef_q2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
             coef_q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1) )
     
               rhs_q2 = 0.0; rhs_q2l = 0.0
               coef_q2 = 0.0; coef_q2l = 0.0

     allocate( stress_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &      !Horizontal tension tensor component
               stress_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) )         !Horizontal shearing tensor component      
     
               stress_t=0.0; stress_s=0.0     
     
     allocate (igrzts(bnd_x1:bnd_x2,bnd_y1:bnd_y2)    )      ! igrz[T,S]= 1 : f = f0
                                                             !          = 2 : df/dz = f0 
               igrzts=0

     allocate( pass_tracer(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),  &   !Passive tracer
               pt_forc_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &   !Passive tracer surface advective forcing
               pt_forc_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &   !Passive tracer surface diffusive forcing
                    igrzpt(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )        !Types of boundary condition

               pass_tracer=0.0; pt_forc_adv=0.0;  pt_forc_dif=0.0; igrzpt=0

     allocate( age(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),        &     !Age
      age_forc_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),           &     !Age surface advective forcing
      age_forc_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),           &     !Age surface diffusive forcing
           igrzage(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )                !Types of boundary condition
               
               age=0.0; age_forc_adv=0.0; age_forc_dif=0.0; igrzage=0 
                                   
     allocate( amts_ref(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &   !T lateral diffusion in T-points [m^2/s]               
               amuv_ref(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &   !U and V 2th order lateral diffusion in T-points[m^2/s]
              amuv4_ref(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)  )       !U and V 4th order lateral viscosity in T-points[m^4/s]^(1/2) 

               amts_ref=0.0; amuv_ref=0.0; amuv4_ref=0.0

     allocate( amts(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &   !T lateral diffusion in T-points [m^2/s]               
               amuv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &   !U and V 2th order lateral diffusion in T-points[m^2/s]
              amuv4(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)  )       !U and V 4th order lateral viscosity in T-points[m^4/s]^(1/2) 

               amts=0.0; amuv=0.0; amuv4=0.0
                                                        
     allocate( anzt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &     !T vertical diffusion [m^2/s]
               anzu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &     !U and V vertical viscosity [m^2/s]
               vbf2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
               shf2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1) )
               
               anzt=0.0; anzu=0.0; vbf2 = 0.0; shf2 = 0.0 

     allocate ( mld_dens(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )          !mixed layer depth [m] based on potential density
               
               mld_dens=0.0

     allocate( tflux_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &       !total advective surface temperature flux [°C*m/s]
               tflux_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &       !total diffusive surface temperature flux [°C*m/s]
              swflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total surface shortwave radiation flux [°C*m/s]
               sflux_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &       !total advective surface salinity flux [psu*m/s]
               sflux_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &       !total diffusive surface salinity flux [psu*m/s]
            vol_flux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total volume flux (precip+evap+ice+rivers) [m/s]
            bot_fric(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !bottom friction rate [m/s]
               cd_ao(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !drag coefficient air-ocean [kg/m^2/s]
               cd_ai(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !drag coefficient air-ice [kg/m^2/s]
               cd_io(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !drag coefficient ice-ocean [kg/m^2/s]              
            divswrad(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &       !shortwave radiation divergence coefficients
                dkft(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !relaxation coefficient for SST, [m/s]
                dkfs(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !relaxation coefficient for SSS, [m/s]
            sensheat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !sensible heat flux
             latheat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !latent heat flux
              lw_bal(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !longwave radiation balance
          sw_bal_atm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !shortwave radiation balance
           sw_bal_oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !shortwave radiation balance
          hf_tot_atm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total heat flux
           hf_tot_oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total heat flux
           sf_tot_oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total salt flux
          wf_tot_atm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total heat flux from air
           wf_tot_oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )              !total water flux into ocean

               tflux_adv=0.0; tflux_dif=0.0; swflux=0.0
               sflux_adv=0.0; sflux_dif=0.0; vol_flux=0.0
               bot_fric=0.0; cd_ao=0.0; cd_ai=0.0; cd_io=0.0
               divswrad=0.0; dkft=0.0; dkfs=0.0
               sensheat=0.0; latheat=0.0; lw_bal=0.0; sw_bal_atm=0.0; sw_bal_oc=0.0
               hf_tot_atm=0.0; hf_tot_oc=0.0; sf_tot_oc=0.0; wf_tot_atm=0.0; wf_tot_oc=0.0

      allocate( tatm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Air temperature, [°C]
                qatm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Air humidity, [kg/kg]
                rain(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !rain, [m/s]
                snow(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !snow, [m/s]
                 lwr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Downward  longwave radiation, [W/m^2]
                 swr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Downward shortwave radiation, [W/m^2]
                slpr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Sea level pressure, [Pa]
                uwnd(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Zonal      wind speed, [m/s]
                vwnd(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Meridional wind speed, [m/s]
                taux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Zonal      wind stress, [Pa]
                tauy(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )       !Meridional wind stress, [Pa]
               
                tatm=0.0; qatm=0.0; rain=0.0; snow=0.0
                 lwr=0.0;  swr=0.0; slpr=0.0; uwnd=0.0; vwnd=0.0
                taux=0.0; tauy=0.0


      allocate( r_diss(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
            AlbOpWater(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )
               
               r_diss = 0.0; AlbOpWater = 0.0
       
      allocate( hice(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &       !Ice volume, [m]
                aice(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &       !Ice compactness, [0-1]
               hsnow(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &       !Snow volume, [m]
                sice(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &       !Ice salinity, [PSU*m]
              mistot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !total ice and snow mass
              aistot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Total ice compactness, [0-1]
               hitot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Total ice volume, m
               hstot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Total snow volume, m
                pice(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Ice pressure, [Pa*m]
               sig_d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Stress tensor divergence component
               sig_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Stress tensor tension component
               sig_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Stress tensor shear component
               eps_d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Strain rate tensor divergence component
               eps_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Strain rate tensor tension component
               eps_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Strain rate tensor shear component
                uice(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &       !Ice zonal velocity, [m/s]
                vice(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  )                !Ice meridional velocity, [m/s]

               hice=0.0; aice=0.0; hsnow=0.0; sice=0.0
               mistot=0.0; aistot=0.0; hitot=0.0; hstot=0.0; pice=0.0
               sig_d=0.0; sig_t=0.0; sig_s=0.0
               eps_d=0.0; eps_t=0.0; eps_s=0.0
               uice=0.0; vice=0.0
      
      allocate( hf_sugar(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )             !heat  flux due to new ice formation
               
               hf_sugar=0.0
      
      allocate(idyn_coefx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
               idyn_coefy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
           idyn_coefx_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
           idyn_coefy_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                idyn_rhsx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                idyn_rhsy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
            idyn_rhsx_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
            idyn_rhsy_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  mistotu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
                  mistotv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &           
                    piceh(bnd_x1:bnd_x2,bnd_y1:bnd_y2))
                
                idyn_coefx=0.0; idyn_coefy=0.0; idyn_coefx_tot=0.0; idyn_coefy_tot=0.0
                idyn_rhsx=0.0; idyn_rhsy=0.0; idyn_rhsx_tot=0.0; idyn_rhsy_tot=0.0
                mistotu=0.0; mistotv=0.0; piceh=0.0

      allocate( taux_atm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &   !Zonal      stress from air, [Pa]
                tauy_atm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &   !Meridional stress from air, [Pa]
                taux_oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &   !Zonal      stress to water, [Pa]
                tauy_oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )       !Meridional stress to water, [Pa]

                taux_atm=0.0; tauy_atm=0.0; taux_oc=0.0; tauy_oc=0.0
           
      allocate(  tx_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                 tx_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                 ty_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                 ty_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                 tx_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
                 ty_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )
               
               tx_adv=0.0; tx_dif=0.0
               ty_adv=0.0; ty_dif=0.0
               tx_coef=0.0; ty_coef=0.0

      allocate( wind_ao(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &   !Air-ocean speed shear, [m/s]
                wind_ai(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &   !Air-ice speed shear, [m/s]
                ioshear(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &   !ice-ocean speed shear, [m/s]
                  windx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &   !Zonal      wind speed in u-points, [m/s]
                windy(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )      !Meridional wind speed in v-points, [m/s]
                
               wind_ao=0.0; wind_ai=0.0; ioshear=0.0; windx=0.0; windy=0.0

      allocate( tt_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                ss_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                uu_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                vv_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
            ssdata_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,11),     &
            ssflux_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,10),     &
             vflux_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,5),      &
               hhq_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               hhu_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               hhv_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               ssh_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
           ssh_hhq_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               sls_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
          mld_dens_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               txo_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               tyo_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
              uwnd_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
              vwnd_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2))
     
               tt_calc=0.0; ss_calc=0.0; uu_calc=0.0; vv_calc=0.0
               ssdata_calc=0.0;  ssflux_calc=0.0;  vflux_calc=0.0  
               hhq_calc=0.0;  hhu_calc=0.0;  hhv_calc=0.0
               ssh_calc=0.0d0; ssh_hhq_calc=0.0d0; sls_calc=0.0
               mld_dens_calc=0.0
               txo_calc=0.0;  tyo_calc=0.0; uwnd_calc=0.0; vwnd_calc=0.0 

      allocate( pt_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) )
               
               pt_calc=0.0

      allocate (age_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz))
               
               age_calc=0.0

      allocate ( den_sgt_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                 den_sg0_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                 den_sgm_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) )

               den_sgt_calc=0.0; den_sg0_calc=0.0; den_sgm_calc=0.0
                   
      allocate( aice_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                hice_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               hsnow_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                uice_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                vice_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )

               aice_calc=0.0; hice_calc=0.0; hsnow_calc=0.0; uice_calc=0.0; vice_calc=0.0
     
               meancalc=0

      allocate( aux_array3d_tgr1(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                aux_array3d_tgr2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                aux_array3d_tgr3(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                aux_array3d_tgr4(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                aux_array3d_tgr5(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                             z3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                      den_isopyc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                aux_array3d_wgr1(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
                aux_array3d_wgr2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
                aux_array3d_wgr3(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
                aux_array3d_wgr4(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
                aux_array3d_zgr_loc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev_loc),    &
                aux_array3d_zgr_glob(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev_glob) )

                aux_array3d_tgr1 = 0.0
                aux_array3d_tgr2 = 0.0
                aux_array3d_tgr3 = 0.0
                aux_array3d_tgr4 = 0.0
                aux_array3d_tgr5 = 0.0
                             z3d = 0.0
                      den_isopyc = 0.0
                aux_array3d_wgr1 = 0.0
                aux_array3d_wgr2 = 0.0
                aux_array3d_wgr3 = 0.0
                aux_array3d_wgr4 = 0.0
                aux_array3d_zgr_loc  = 0.0
                aux_array3d_zgr_glob = 0.0

      allocate( aux_array_itr1(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),     &
                aux_array_itr2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),     &
                aux_array_itr3(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),     &
                aux_array_itr4(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad) )

                aux_array_itr1 = 0.0
                aux_array_itr2 = 0.0
                aux_array_itr3 = 0.0
                aux_array_itr4 = 0.0
      
      allocate( ssd(bnd_x1:bnd_x2,bnd_y1:bnd_y2,11),        &
                ssf(bnd_x1:bnd_x2,bnd_y1:bnd_y2,10),        &
                vfl(bnd_x1:bnd_x2,bnd_y1:bnd_y2,5) )
                
                ssd=0.0; ssf=0.0; vfl=0.0

      allocate(aux_array2d_01(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
               aux_array2d_02(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )
               
               aux_array2d_01 = 0.0;  aux_array2d_02 = 0.0

   endsubroutine ocean_variables_allocate

!deallocation of arrays
!=======================================================================
   subroutine ocean_variables_deallocate
   use ocean_variables

      deallocate(aux_array2d_02,aux_array2d_01)
      
      deallocate(vfl,ssf,ssd)
      
      deallocate(aux_array_itr4,aux_array_itr3,aux_array_itr2,     &
                 aux_array_itr1)
                
      deallocate(aux_array3d_zgr_glob,aux_array3d_zgr_loc,        &
                          aux_array3d_wgr4,aux_array3d_wgr3,      &
               aux_array3d_wgr2,aux_array3d_wgr1,den_isopyc,      &
                      z3d,aux_array3d_tgr5,aux_array3d_tgr4,      &
               aux_array3d_tgr3,aux_array3d_tgr2,aux_array3d_tgr1)

      deallocate(vice_calc, uice_calc, hsnow_calc, hice_calc, aice_calc)
      deallocate(den_sgm_calc, den_sg0_calc, den_sgt_calc)
      deallocate(age_calc)
      deallocate(pt_calc)

      deallocate( vwnd_calc, uwnd_calc, tyo_calc, txo_calc, mld_dens_calc,   &
             sls_calc, ssh_hhq_calc, ssh_calc, hhv_calc, hhu_calc, hhq_calc, &
             vflux_calc,ssflux_calc,ssdata_calc,vv_calc,uu_calc,ss_calc,tt_calc)

      deallocate(windy, windx, ioshear, wind_ai, wind_ao)
      deallocate(ty_coef,tx_coef,ty_dif,ty_adv,tx_dif,tx_adv)
      deallocate(tauy_oc,taux_oc,tauy_atm,taux_atm)

      deallocate(piceh,mistotv,mistotu,idyn_rhsy_tot,idyn_rhsx_tot,      &
                 idyn_rhsy,idyn_rhsx,idyn_coefy_tot,idyn_coefx_tot,      &
                 idyn_coefy,idyn_coefx)

      deallocate(hf_sugar)

      deallocate(vice, uice, eps_s, eps_t, eps_d, sig_s, sig_t, sig_d, &
                 pice, hstot, hitot, mistot,sice,hsnow,aice,hice) 

      deallocate(AlbOpWater, r_diss)
      deallocate(tauy,taux,vwnd,uwnd,slpr,swr,lwr,snow,rain,qatm,tatm)

      deallocate(wf_tot_oc,wf_tot_atm,sf_tot_oc,hf_tot_oc,hf_tot_atm,           &
                 sw_bal_oc,sw_bal_atm,lw_bal,latheat,sensheat,                  &
                 dkfs,dkft,divswrad,cd_io,cd_ai,cd_ao,bot_fric,vol_flux,        &
                 sflux_dif,sflux_adv,swflux,tflux_dif,tflux_adv)             
      deallocate(mld_dens)
      deallocate(shf2,vbf2,anzu,anzt,amuv4,amuv,amts)
      deallocate(amuv4_ref,amuv_ref,amts_ref)
      deallocate(igrzage,age_forc_dif,age_forc_adv,age)   
      deallocate(igrzpt,pt_forc_dif,pt_forc_adv,pass_tracer)   
      deallocate(igrzts)
      deallocate(stress_s, stress_t)

      deallocate(coef_q2l, coef_q2, rhs_q2l, rhs_q2)   
      deallocate(q2l,q2)

      deallocate(yng_slp,den_ini,den_sgm,den_sg0,den_sgt,ss,tt,ww,vh,uh,vv,uu)
      deallocate(sls)
      deallocate(RHSy2d,RHSx2d,vbrtr,ubrtr,ssh_n,ssh)
   
   endsubroutine ocean_variables_deallocate
   
!Allocation of arrays
!=======================================================================
  subroutine oceanbc_arrays_allocate
  use ocean_bc
  use basin_grid
  
   allocate ( sst_obs(nx, ny),     &       !Observed SST [�C]
              sss_obs(nx, ny),     &       !Observed SSS [psu-35ppt])
               runoff(nx, ny),     &       !River runoff [kg/m^2/s]
         runoff_solid(nx, ny)  )           !Solid runoff [kg/m^2/s]

               sst_obs=0.0; sss_obs=0.0
               runoff=0.0; runoff_solid=0.0
         
   allocate( lqpx(numb_of_lqp_max),        &
             lqpy(numb_of_lqp_max),        &
             tlqbw(numb_of_lqp_max,nz),    &
             slqbw(numb_of_lqp_max,nz),    &
             ulqbw(numb_of_lqp_max,nz),    &
             vlqbw(numb_of_lqp_max,nz),    &
           sshlqbw(numb_of_lqp_max) )
               
               lqpx=0; lqpy=0; tlqbw=0.0; slqbw=0.0
               ulqbw=0.0; vlqbw=0.0; sshlqbw=0.0

   allocate( agelqbw(numb_of_lqp_max,nz) )
               
               agelqbw=0.0
         
   allocate( index_of_lb(numb_of_lqp_max) )
   
               index_of_lb=' '

  endsubroutine oceanbc_arrays_allocate

!=======================================================================
  subroutine oceanbc_arrays_deallocate
  use ocean_bc

   deallocate(index_of_lb)
   deallocate(agelqbw)
   deallocate(sshlqbw,vlqbw,ulqbw,slqbw,tlqbw,lqpy,lqpx)
   deallocate(runoff_solid, runoff,sss_obs,sst_obs)

  endsubroutine oceanbc_arrays_deallocate

!=======================================================================
subroutine atm_arrays_allocate
use atm_forcing

 allocate( xa(nxa),ya(nya), atm_mask(nxa,nya) )
               
               xa=0.0d0; ya=0.0d0; atm_mask=0

 allocate( a_hflux(nxa,nya),       &   !heat balance [w/m**2] 
           a_swrad(nxa,nya),       &   !sw radiation balance[w/m**2]
           a_wflux(nxa,nya),       &   !precipitation-evaporation[m/s]
           a_stress_x(nxa,nya),    &   !zonal wind stress[pA=n/m**2] 
           a_stress_y(nxa,nya),    &   !meridional wind stress[pA=n/m**2]
           a_slpr(nxa,nya),        &   !pressure at sea surface
           a_lwr(nxa,nya),         &   !dw-lw-rad[w/m**2]
           a_swr(nxa,nya),         &   !dw-sw-rad[w/m**2]
           a_rain(nxa,nya),        &   !precipit[m/s]
           a_snow(nxa,nya),        &   !precipit[m/s]
           a_tatm(nxa,nya),        &   !temp of atmosphere[�c]
           a_qatm(nxa,nya),        &   !humidity [g/kg]
           a_uwnd(nxa,nya),        &   !u-wind speed[m/s]
           a_vwnd(nxa,nya)  )          !v-wind speed[m/s]

               a_hflux=0.0; a_swrad=0.0; a_wflux=0.0 
               a_stress_x=0.0; a_stress_y=0.0
               a_slpr=0.0; a_lwr=0.0; a_swr=0.0
               a_rain=0.0; a_snow=0.0
               a_tatm=0.0; a_qatm=0.0 
               a_uwnd=0.0; a_vwnd=0.0

               ind_change_heat =0
               ind_change_water=0
               ind_change_stress=0
               ind_change_rad  =0
               ind_change_prec =0
               ind_change_tatm =0
               ind_change_qatm =0
               ind_change_wind =0
               ind_change_slpr =0
               
               num_rec_heat =0
               num_rec_water=0
               num_rec_stress=0
               num_rec_rad  =0
               num_rec_prec =0
               num_rec_tatm =0
               num_rec_qatm =0
               num_rec_wind =0
               num_rec_slpr =0

endsubroutine atm_arrays_allocate

!=======================================================================
subroutine atm_arrays_deallocate
use atm_forcing

 deallocate(a_vwnd,a_uwnd,a_qatm,a_tatm,a_snow,a_rain,a_swr,a_lwr,     &
             a_slpr,a_stress_y,a_stress_x,a_wflux,a_swrad,a_hflux)
 deallocate(atm_mask,ya,xa)

endsubroutine atm_arrays_deallocate

!=====================================================================================
 subroutine atm2oc_allocate
 use atm2oc_interpol

      allocate( wght_mtrx_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &
                  i_input_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &
                  j_input_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)    )
               wght_mtrx_a2o=0.0d0; i_input_a2o=0; j_input_a2o=0

 endsubroutine atm2oc_allocate

!=======================================================================
 subroutine atm2oc_deallocate
 use atm2oc_interpol

      deallocate(j_input_a2o, i_input_a2o, wght_mtrx_a2o)
 
 endsubroutine atm2oc_deallocate

 endmodule init_arrays_routes