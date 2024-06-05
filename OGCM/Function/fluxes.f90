module flux_routes
use constants
use basin_grid
use ocean_variables
use ocean_bc
use key_switches
use ocalg_routes
use seaice_routes
implicit none

include 'atmforcing.fi'

contains
!========================================================================
subroutine sea_surface_fluxes(tau,nicestep)

real(4), intent(in):: tau
integer, intent(in):: nicestep

integer m,n,k,grad
integer ierr

real(4) tf, wf, hf, sf, sh, lh, ev, lw ,sw, tx, ty, qu, qv, aopw, tau2, drag
real(4) wf_calc, wf_ave, sf_calc, sf_ave, m_calc, saice, frac, vflux, sf1,sf2
real(4) ut,vt,uit,vit
real(4) buf_real4

m_calc=0.0
wf_calc=0.0
sf_calc=0.0
wf_ave=0.0
sf_ave=0.0
   
 if(ksw_wflux>0) then
   !$omp parallel do private(m,n) reduction(+:m_calc) 
   do n=ny_start,ny_end
     do m=nx_start,nx_end
      if(lu(m,n)>0.5) then
       m_calc = m_calc + sqt(m,n)
      endif    
     enddo
   enddo
   !$omp end parallel do
   buf_real4 = m_calc
   call mpi_allreduce(buf_real4,m_calc,1,mpi_real4,mpi_sum,cart_comm,ierr)
 endif

!If 1st boundary condition (SST=SST_obs, SSS=SSS_obs, Tx and Ty are set)
if(ksw_ssbc==1) then

!$omp parallel do private(m,n)
 do n=ny_start,ny_end
   do m=nx_start,nx_end 
     if(lu(m,n)>0.5) then
        vol_flux(m,n) = 0.0
       tflux_adv(m,n) = 0.0
       tflux_dif(m,n) = sst_obs(m,n)
       sflux_adv(m,n) = 0.0
       sflux_dif(m,n) = sss_obs(m,n)
          swflux(m,n) = 0.0
     endif
   enddo
 enddo
!$omp end parallel do

endif

!If 2nd boundary condition (Qt, Qs, Tx and Ty are set)
if(ksw_ssbc==2) then
  wf_tot_atm = wf_tot_oc
endif

!If 1st or 2nd boundary condition (Tx and Ty are set)
if(ksw_ssbc<3) then

!$omp parallel do private(m,n)
  do n=ny_start,ny_end
    do m=nx_start,nx_end 

      if (lcu(m,n)>0.5) then
       taux_atm(m,n)=(taux(m,n)*dy(m,n)+taux(m+1,n)*dy(m+1,n))/(2.0*dyh(m,n))
       taux_oc(m,n) = taux_atm(m,n)
        tx_dif(m,n) = taux_oc(m,n)
        tx_coef(m,n)= 0.0
      endif

      if (lcv(m,n)>0.5) then
       tauy_atm(m,n)=(tauy(m,n)*dx(m,n)+tauy(m,n+1)*dx(m,n+1))/(2.0*dxh(m,n))
       tauy_oc(m,n) = tauy_atm(m,n)
        ty_dif(m,n) = tauy_oc(m,n)
        ty_coef(m,n)= 0.0
      endif
       
    enddo
  enddo
!$omp end parallel do

endif

if(ksw_ssbc==3) then

!Computing wind velocity in C-grid velocity points
!$omp parallel do private(m,n,ut,vt,uit,vit)
  do n=ny_start, ny_end
    do m=nx_start, nx_end
     
     if(lcu(m,n)>0.5) then
      windx(m,n)= (uwnd(m,n)*dy(m,n)+uwnd(m+1,n)*dy(m+1,n))/(2.0*dyh(m,n))
     endif
     
     if(lcv(m,n)>0.5) then
      windy(m,n)= (vwnd(m,n)*dx(m,n)+vwnd(m,n+1)*dx(m,n+1))/(2.0*dxh(m,n))
     endif

     if(lu(m,n)>0.5) then
      ut = (uu(m  ,n,1)*dyh(m  ,n)*hhu(m  ,n)    &
          + uu(m-1,n,1)*dyh(m-1,n)*hhu(m-1,n))/(2.0*hhq(m,n)*dy(m,n))       !u on t-grid
      vt = (vv(m,n  ,1)*hhv(m,n  )*dxh(m,n  )    &
          + vv(m,n-1,1)*hhv(m,n-1)*dxh(m,n-1))/(2.0*hhq(m,n)*dx(m,n))       !v on t-grid   
      
      uit = (uice(m  ,n)*dyh(m  ,n)              &
           + uice(m-1,n)*dyh(m-1,n) )/(2.0*dy(m,n))       !u on t-grid
      vit = (vice(m,n  )*dxh(m,n  )              &
           + vice(m,n-1)*dxh(m,n-1) )/(2.0*dx(m,n))       !v on t-grid  
   
      wind_ao(m,n) = sqrt( (uwnd(m,n)-ut )**2 + (vwnd(m,n)-vt )**2 )
      wind_ai(m,n) = sqrt( (uwnd(m,n)-uit)**2 + (vwnd(m,n)-vit)**2 )
      ioshear(m,n) = sqrt( (ut-uit)**2 + (vt-vit)**2 )

     endif

    enddo
  enddo
!$omp end parallel do

  call syncborder_real(windx, 1)
  call syncborder_real(windy, 1)
  call syncborder_real(ioshear, 1)

  if(periodicity_x/=0) then
    call cyclize_x(windx,nx,ny,1,mmm,mm)
    call cyclize_x(windy,nx,ny,1,mmm,mm)
    call cyclize_x(ioshear,nx,ny,1,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call cyclize_y(windx,nx,ny,1,nnn,nn)
    call cyclize_y(windy,nx,ny,1,nnn,nn)
    call cyclize_y(ioshear,nx,ny,1,nnn,nn)
  endif

!Computing fluxes over open water

!$omp parallel do private(m, n, wf, hf, sf, sh, lh, ev, lw ,sw, aopw) 
  do n=ny_start,ny_end
    do m=nx_start,nx_end 
      if(lu(m,n)>0.5) then
      
        call air_sea_turbulent_fluxes(wind_ao(m,n),      &   ! wind modulo, m/s
                                      slpr(m,n),         &   ! sea level pressure, Pa
                                      tatm(m,n),         &   ! air temperature, °C
                                      tt(m,n,1),         &   ! sea surface temp, °C
                            salref +  ss(m,n,1),         &   ! sea surface salinity, PSU
                                      qatm(m,n),         &   ! air specific humidity, kg/kg
                                       u_height,         &   ! height of wind datasets, m
                                       t_height,         &   ! height of tair datasets, m
                                       q_height,         &   ! height of qair datasets, m
                                             sh,         &   ! sensible heat flux, W/m^2
                                             ev,         &   ! evaporation rate, kg/m^2/s
                                             lh,         &   ! latent heat flux, W/m^2
                                      cd_ao(m,n)  )          ! drag coefficient [kg/m^2/s]
             
        aopw = 1.0-aistot(m,n)
        sh = sh*aopw
        ev = ev*aopw
        lh = lh*aopw
        lw=EmissWater*(lwr(m,n)-StephBoltz*(tt(m,n,1)+273.15)**4)*aopw
        sw=swr(m,n)*(1.0-AlbOpWater(m,n))*aopw

        sensheat(m,n)= sh
        latheat(m,n) = lh
        lw_bal(m,n)=   lw
        sw_bal_atm(m,n)=sw
        sw_bal_oc(m,n) =sw 

        wf= rain(m,n) + snow(m,n)*aopw + ev + runoff(m,n) + runoff_solid(m,n)
        hf= sh + lh + lw + sw -(snow(m,n)*aopw + runoff_solid(m,n))*lambda_f 
                
        wf_tot_atm(m,n)= wf 
        wf_tot_oc(m,n) = wf

        hf_tot_atm(m,n)= hf
        hf_tot_oc(m,n) = hf
        sf_tot_oc(m,n) = 0.0

      endif    
    enddo
  enddo
!$omp end parallel do
  
  call syncborder_real(cd_ao, 1)

  if(periodicity_x/=0) then
    call cyclize_x(cd_ao,nx,ny,1,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call cyclize_y(cd_ao,nx,ny,1,nnn,nn)
  endif

!Wind stress computation

 !$omp parallel do private(m,n,drag) 
  do n=ny_start,ny_end
    do m=nx_start,nx_end 
      
      if(lcu(m,n)>0.5) then
        drag=(cd_ao(m,n)+cd_ao(m+1,n))/2.0
        taux_atm(m,n)=drag*(windx(m,n)-uu(m,n,1))
        taux_oc(m,n)=taux_atm(m,n)
         tx_dif(m,n)= drag*windx(m,n)
        tx_coef(m,n)= drag
      endif    
      
      if(lcv(m,n)>0.5) then
        drag=(cd_ao(m,n)+cd_ao(m,n+1))/2.0
        tauy_atm(m,n)=drag*(windy(m,n)-vv(m,n,1))
        tauy_oc(m,n)=tauy_atm(m,n)
         ty_dif(m,n)= drag*windy(m,n)
        ty_coef(m,n)= drag
      endif 
    enddo
  enddo

!Computing total ice compactness, volume and pressure for dynamics block (total mass is computed before)
!$omp parallel do private(m,n,grad)
do n=ny_start-1, ny_end+1
 do m=nx_start-1, nx_end+1
   if(lu(m,n)>0.5) then
      cd_ai(m,n) = 0.0
      cd_io(m,n) = 0.0
      hf_sugar(m,n) = 0.0
   endif
 enddo
enddo
!$omp end parallel do 
     
 !sea ice thermodynamics
 if(ksw_ith>0) then
  call seaice_thermodynamics(tau) 
 endif

 !sea ice dynamics
 if(ksw_idn>0) then
  call seaice_dynamics(tau,nicestep) 
 endif

  !sea ice transport
 if(ksw_itr>0) then
  call trans_ice_mpdata(aice, uice,vice,tau,aux_array_itr1,aux_array_itr2,    &
                                            aux_array_itr3,aux_array_itr4)
  call trans_ice_mpdata(hice, uice,vice,tau,aux_array_itr1,aux_array_itr2,    &
                                            aux_array_itr3,aux_array_itr4)
  call trans_ice_mpdata(hsnow,uice,vice,tau,aux_array_itr1,aux_array_itr2,    &
                                            aux_array_itr3,aux_array_itr4)
  call trans_ice_mpdata(sice, uice,vice,tau,aux_array_itr1,aux_array_itr2,    &
                                            aux_array_itr3,aux_array_itr4)
 endif
 
  !$omp parallel do private(m,n,grad) 
  do n=ny_start,ny_end
    do m=nx_start,nx_end
      if(lu(m,n)>0.5) then
       do grad=1,mgrad
     
        if(hice(m,n,grad)<himin.or.aice(m,n,grad)<0.0) then
         wf_tot_oc(m,n)= wf_tot_oc(m,n)+hice(m,n,grad)*den_ice/tau+hsnow(m,n,grad)*den_snow/tau        
         sf_tot_oc(m,n)= sf_tot_oc(m,n)+sice(m,n,grad)*den_ice/tau        
!        hf_tot_oc(m,n)= hf_tot_oc(m,n)-(hice(m,n,grad)*den_ice+hsnow(m,n,grad)*den_snow)*lambda_f/tau

          hice(m,n,grad)=0.0
         hsnow(m,n,grad)=0.0
          aice(m,n,grad)=0.0
          sice(m,n,grad)=0.0
        endif

        if(sice(m,n,grad)<0.0) then
         sf_tot_oc(m,n)= sf_tot_oc(m,n)+sice(m,n,grad)*den_ice/tau        
         sice(m,n,grad)=0.0
        endif

        if(hsnow(m,n,grad)<hsmin) then
         wf_tot_oc(m,n)= wf_tot_oc(m,n)+hsnow(m,n,grad)*den_snow/tau        
!        hf_tot_oc(m,n)= hf_tot_oc(m,n)-hsnow(m,n,grad)*den_snow*lambda_f/tau
         hsnow(m,n,grad)=0.0
        endif
               
       enddo
      endif    
    enddo
  enddo

 !$omp parallel do private(m,n,grad,saice,frac) 
  do n=ny_start,ny_end
    do m=nx_start,nx_end 
      if(lu(m,n)>0.5) then
      
      saice=0.0
     do grad=1,mgrad
      saice=saice+aice(m,n,grad)
     enddo

     if(saice>aimin) then
      frac=min(saice, aimax)/saice
      aice(m,n,:)=aice(m,n,:)*frac
     endif

      endif    
    enddo
  enddo

  call syncborder_real(hice , mgrad)
  call syncborder_real(sice , mgrad)
  call syncborder_real(hsnow, mgrad)
  call syncborder_real(aice , mgrad)

  if(periodicity_x/=0) then
    call cyclize_x(hice ,nx,ny,mgrad,mmm,mm)
    call cyclize_x(sice ,nx,ny,mgrad,mmm,mm)
    call cyclize_x(hsnow,nx,ny,mgrad,mmm,mm)
    call cyclize_x(aice ,nx,ny,mgrad,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call cyclize_y(hice ,nx,ny,mgrad,nnn,nn)
    call cyclize_y(sice ,nx,ny,mgrad,nnn,nn)
    call cyclize_y(hsnow,nx,ny,mgrad,nnn,nn)
    call cyclize_y(aice ,nx,ny,mgrad,nnn,nn)
  endif

!$omp parallel do private(m,n,grad) 
    do n = ny_start-1,ny_end+1
      do m = nx_start-1,nx_end+1
        if (lu(m,n)>0.5) then
         aistot(m,n) = 0.0
          hitot(m,n) = 0.0
          hstot(m,n) = 0.0        
         do grad=1,mgrad
           aistot(m,n) = aistot(m,n) + aice(m,n,grad)
            hitot(m,n) =  hitot(m,n) + hice(m,n,grad)
            hstot(m,n) =  hstot(m,n) + hsnow(m,n,grad)
         enddo
         mistot(m,n) = den_ice*hitot(m,n)+den_snow*hstot(m,n)
        endif
      enddo
    enddo
!$omp end parallel do

endif

! 2nd and 3rd conditions (fluxes are used)
if(ksw_ssbc>=2) then

  sf_calc=0.0
  wf_calc=0.0

!fixed water budget
  if(variable_volume_budget==0) then

!$omp parallel do private(m,n,vflux,sf1,sf2) reduction(+:sf_calc)
   do n=ny_start,ny_end
     do m=nx_start,nx_end 
      if (lu(m,n)>0.5) then
           vflux= wf_tot_oc(m,n)/RefDen
           sf1 = -vflux*(salref+ss(m,n,1))
           sf2 = dkfs(m,n)*(sss_obs(m,n)-ss(m,n,1))
           sflux_adv(m,n) = sf1 
           sflux_dif(m,n) = sf_tot_oc(m,n)/RefDen + sf2
        vol_flux(m,n) = 0.0
          swflux(m,n) = rpart*sw_bal_oc(m,n)/(HeatCapWater*RefDen)       
           tflux_adv(m,n) = 0.0
           tflux_dif(m,n) = (hf_tot_oc(m,n)-rpart*sw_bal_oc(m,n)-hf_sugar(m,n))/(HeatCapWater*RefDen)     & 
                             + dkft(m,n)*(sst_obs(m,n)-tt(m,n,1))
           pt_forc_adv(m,n) = -vflux * pass_tracer(m,n,1)
           sf_calc = sf_calc +  (sf1 + sf2)*sqt(m,n)
      endif
     enddo
   enddo
!$omp end parallel do
   
   if(ksw_wflux>0) then
    buf_real4 = sf_calc
    call mpi_allreduce(buf_real4, sf_calc, 1, mpi_real4, mpi_sum, cart_comm, ierr)

    sf_ave=sf_calc/m_calc

!$omp parallel do private(m,n)
    do n=ny_start,ny_end
      do m=nx_start,nx_end 
       if (lu(m,n)>0.5) then
        sflux_dif(m,n) = sflux_dif(m,n) - sf_ave
       endif
      enddo
    enddo
!$omp end parallel do
   endif
  
  endif

!variable water budget
  if(variable_volume_budget>=1) then
 
  aux_array2d_01 = 0.0

    if(ksw_wflux>0) then
 !$omp parallel do private(m,n) reduction(+:wf_calc)
     do n=ny_start,ny_end
       do m=nx_start,nx_end 
        if (lu(m,n)>0.5) then
         aux_array2d_01(m,n) = - dkfs(m,n)*(sss_obs(m,n)-ss(m,n,1))/max(salref+ss(m,n,1),0.01)
         wf_calc = wf_calc + (wf_tot_atm(m,n)/RefDen + aux_array2d_01(m,n))*sqt(m,n)
        endif
       enddo
     enddo
 !$omp end parallel do
     buf_real4 = wf_calc
     call mpi_allreduce(buf_real4, wf_calc, 1, mpi_real4, mpi_sum, cart_comm, ierr)         
     wf_ave=wf_calc/m_calc
    endif
 
 !$omp parallel do private(m,n) 
    do n=ny_start,ny_end
      do m=nx_start,nx_end 
       if (lu(m,n)>0.5) then
         vol_flux(m,n) = wf_tot_oc(m,n)/RefDen + aux_array2d_01(m,n) - wf_ave
            sflux_adv(m,n) = - vol_flux(m,n)*salref
            sflux_dif(m,n) =   sf_tot_oc(m,n)/RefDen
           swflux(m,n) = rpart*sw_bal_oc(m,n)/(HeatCapWater*RefDen)       
            tflux_adv(m,n) =  vol_flux(m,n)*tt(m,n,1) 
            tflux_dif(m,n) = (hf_tot_oc(m,n)-rpart*sw_bal_oc(m,n)-hf_sugar(m,n))/(HeatCapWater*RefDen)       &
                       + dkft(m,n)*(sst_obs(m,n)-tt(m,n,1))
       endif
      enddo
    enddo
 !$omp end parallel do
 
  endif

  call syncborder_real(vol_flux, 1)

  if(periodicity_x/=0) then
    call cyclize_x(vol_flux,nx,ny,1,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call cyclize_y(vol_flux,nx,ny,1,nnn,nn)
  endif

endif

!$omp parallel do private(m,n,qu,qv)
  do n=ny_start,ny_end
    do m=nx_start,nx_end 

      if (lcu(m,n)>0.5) then
       qu=( vol_flux(m  ,n)*sqt(m  ,n)       &
          + vol_flux(m+1,n)*sqt(m+1,n) )/(2.0*squ(m,n))
       
       tx_adv(m,n) = qu*uu(m,n,1)
       tx_dif(m,n) = tx_dif(m,n)/RefDen
       tx_coef(m,n)= tx_coef(m,n)/RefDen
      endif

      if (lcv(m,n)>0.5) then
       qv=( vol_flux(m,n  )*sqt(m,n  )       &
          + vol_flux(m,n+1)*sqt(m,n+1) )/(2.0*sqv(m,n))
       ty_adv(m,n) = qv*vv(m,n,1)
       ty_dif(m,n) = ty_dif(m,n)/RefDen
       ty_coef(m,n)= ty_coef(m,n)/RefDen      
      endif
       
    enddo
  enddo
!$omp end parallel do

call syncborder_real(taux_oc, 1)
call syncborder_real(tauy_oc, 1)

if(periodicity_x/=0) then
  call cyclize_x(taux_oc,nx,ny,1,mmm,mm)
endif

if(periodicity_y/=0) then
  call cyclize_y(tauy_oc,nx,ny,1,nnn,nn)
endif

!Boundary condition for turbulent kinetic energy
if(iabs(ksw_vert)>=2) then
  !$omp parallel do private(m,n,tau2)
   do n=ny_start,ny_end
     do m=nx_start,nx_end 
       if(lu(m,n)>0.5) then
         tau2= (taux_oc(m,n)**2 + taux_oc(m-1,n)**2        &
               +tauy_oc(m,n)**2 + tauy_oc(m,n-1)**2)/2.0
         q2(m,n,1)= sqrt(tau2)/RefDen *B1_t**(2.0/3.0)
       endif
     enddo
   enddo
  !$omp end parallel do 
endif

endsubroutine sea_surface_fluxes
!==============================================================================================
subroutine sea_bottom_fluxes

integer m, n
real(4) u2, Cb, dzpzf

!$omp parallel do private(m,n,u2,Cb,dzpzf)          
        do n=ny_start,ny_end
         do m=nx_start,nx_end

          if(lu(m,n)>0.5) then
           u2 = (uu(m,n,nz)**2 + uu(m-1,n,nz)**2      &
                +vv(m,n,nz)**2 + vv(m,n-1,nz)**2 )/2.0
           dzpzf = max(hzt(nz+1)*hhq(m,n)/zbot, 3.0)
           Cb = max( min( (vonkarman/log(dzpzf))**2, Cb_max), Cb_min)
           bot_fric(m,n) = Cb * sqrt(u2 + Eb)
           q2(m,n,nz+1) = bot_fric(m,n) * sqrt(u2) * B1_t**(2.0/3.0)
          endif
          
         enddo
        enddo
 !$omp end parallel do        

call syncborder_real(bot_fric, 1)

if(periodicity_x/=0) then
  call cyclize_x(bot_fric,nx,ny,1,mmm,mm)
endif

if(periodicity_y/=0) then
  call cyclize_y(bot_fric,nx,ny,1,nnn,nn)
endif

endsubroutine sea_bottom_fluxes

endmodule flux_routes
