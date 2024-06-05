module seaice_routes
use constants
use basin_grid
use ocalg_routes
implicit none

contains 
!==========================================================
subroutine seaice_thermodynamics(tau)
use ocean_variables
use key_switches

real(4),intent(in):: tau

integer m,n,k,iter,grad,melt,mark
real(4) hi, hs, tf, tsurf, t_is, sal_oc, depth, href, hsug
real(4) kiceph, ksnowph, hbal, sink, spart, ai, as
real(4) dmsnow, dmice, dmoc, dqoc, dmsice, daice, dmi, dms, dqb
real(4) hsnew, hinew, hiold, qiw, heat_top, heat_bot, swi, swb, swpen_top, swpen_bot, vol
real(4) cd, cd_rt, ustar, u, log_arg, zstar, hi1
real(4) sh, lh, lw, ev, cdai
real(4) t1,t2,tmid
 

!$omp parallel do private(m, n, k, grad, sal_oc, href, mark,          &
!$omp     hbal,iter, tsurf,tf, u, heat_top, heat_bot, dmi, dms, dqb,  &
!$omp     dmsnow, dmice, sink, dmsice, hi, kiceph, hsug,              &
!$omp     swpen_top, swpen_bot, hs, spart, melt, ksnowph,             &
!$omp     cdai, t_is, ai, as, swb, swi,sh, lh, lw, ev, zstar, hi1,    &
!$omp     hsnew, hinew, dmoc, dqoc, ustar, qiw, cd, cd_rt, log_arg,   &
!$omp     hiold, daice, vol, depth, t1, t2, tmid)
do n=ny_start, ny_end              !loop on latitude
  do m=nx_start, nx_end            !loop on longitude
   if(lu(m,n)>0.5) then            !check water area

   sal_oc=max(ss(m,n,1)+salref,0.0)
   tf = freeze_temp(sal_oc)

   u = max(ioshear(m,n), eb_rt)
   
   href = ( href_min*(tatm(m,n)-tref_min)   &
           +href_max*(tref_max-tatm(m,n))  )/(tref_max-tref_min)
   href = max(min(href,href_max),href_min)
            
     do grad=1, mgrad               !loop over gradations

      heat_top=0.0
      heat_bot=0.0

!if ice is available      
      if(aice(m,n,grad)>aimin) then       !check of ice is available

!   snowfall over ice
       dmsnow=snow(m,n)*aice(m,n,grad)*tau
       hsnow(m,n,grad)=hsnow(m,n,grad)+dmsnow/den_snow
       wf_tot_atm(m,n)=wf_tot_atm(m,n)+dmsnow/tau

!snow to ice gravity transformation
       if(hsnow(m,n,grad)>hsmax*aice(m,n,grad)) then
        dmsnow=-gamma_stoice*tau*hsnow(m,n,grad)*den_snow
        dmice=-dmsnow
        hsnow(m,n,grad)=hsnow(m,n,grad)+dmsnow/den_snow
        hice(m,n,grad)=hice(m,n,grad)+dmice/den_ice
       endif

!snow to ice sinking transformation
       sink=(hsnow(m,n,grad)*den_snow+hice(m,n,grad)*den_ice)/RefDen
       
       if(sink>hice(m,n,grad)) then
      
        dmice= (sink-hice(m,n,grad))*den_ice
        dmsnow=-dmice
        dmsice= dmice*sal_oc*freeze_frac_sal
       
        hsnow(m,n,grad)=hsnow(m,n,grad)+dmsnow/den_snow
        hice(m,n,grad)=hice(m,n,grad)+dmice/den_ice
        sice(m,n,grad)=sice(m,n,grad)+dmsice/den_ice
        sf_tot_oc(m,n) = sf_tot_oc(m,n) - dmsice/tau

       endif
             
       hi=hice(m,n,grad)/aice(m,n,grad)
       kiceph=kice/hi
       swpen_top=i0*exp(-hi/h_swi)
       swpen_bot=swpen_top
       
       hs=hsnow(m,n,grad)/aice(m,n,grad)
       spart=hs/(hs+hsnow_steak)
    
 !----------------------------------------------------------------      
       if(hsnow(m,n,grad)>hsmin) then         !if ice is covered with snow
         
         melt=0
         
         ksnowph=ksnow/hs
         
         t1=-70.0
         t2= 1.0

         hbal= 1000.0
         iter=0

        do while(abs(t2-t1)>accur_tem.and.abs(hbal)>accur_flux.and.iter<=30)

         iter=iter+1
         tmid=(t1+t2)/2.0
         tsurf=tmid
         t_is=(tsurf*ksnowph+tf*kiceph)/(ksnowph+kiceph)         

         call  hbal_snow(hbal,   &
                 wind_ai(m,n),   &
                    slpr(m,n),   &
                    tatm(m,n),   &
                    qatm(m,n),   &
                     lwr(m,n),   &
                     swr(m,n),   &
                           hi,   &
                        spart,   &
                    swpen_top,   &
                    swpen_bot,   &
                      ksnowph,   &
                        tsurf,   &
                         t_is,   &
                           sh,   &
                           ev,   &
                           lh,   &
                         cdai,   &
                           lw,   &
                          swb,   &
                          swi,   &
                  AlbOpWater(m,n))         
                    
        if(hbal>0.0) then
          t1=tmid
        else
          t2=tmid
        endif
            
        enddo
if (iter>30) write(*,'(i2,2(a,i5),4(a,f12.5), a,i3)') 1, ' m=', m, ' n=', n, ' t1=', t1, ' t2=', t2, ' tsurf=', tsurf, ' hbal=', hbal, ' iter=', iter
        
        if(tsurf>tmelt) then   !recomputing fluxes and snow melting if temperature iv over tmelt
            melt=1
          
            tsurf=tmelt
            t_is=(tsurf*ksnowph+tf*kiceph)/(ksnowph+kiceph)  

         call  hbal_snow(hbal,   &
                 wind_ai(m,n),   &
                    slpr(m,n),   &
                    tatm(m,n),   &
                    qatm(m,n),   &
                     lwr(m,n),   &
                     swr(m,n),   &
                           hi,   &
                        spart,   &
                    swpen_top,   &
                    swpen_bot,   &
                      ksnowph,   &
                        tsurf,   &
                         t_is,   &
                           sh,   &
                           ev,   &
                           lh,   &
                         cdai,   &
                           lw,   &
                          swb,   &
                          swi,   &
                  AlbOpWater(m,n))   
                           
        endif
!end of recomputing fluxes
     
         sensheat(m,n)=sensheat(m,n)+sh*aice(m,n,grad)
          latheat(m,n)= latheat(m,n)+lh*aice(m,n,grad)
           lw_bal(m,n)=  lw_bal(m,n)+lw*aice(m,n,grad)
       sw_bal_atm(m,n)=  sw_bal_atm(m,n)+swb*aice(m,n,grad)
       hf_tot_atm(m,n)=  hf_tot_atm(m,n)+(sh+lh+lw+swb)*aice(m,n,grad)
        sw_bal_oc(m,n)=   sw_bal_oc(m,n)+swi*aice(m,n,grad)
        hf_tot_oc(m,n)=   hf_tot_oc(m,n)+swi*aice(m,n,grad)
            cd_ai(m,n)= cd_ai(m,n) + cdai*aice(m,n,grad)

!  snow sublimation/desublimation
         dmsnow=aice(m,n,grad)*ev*tau
         hsnew=hsnow(m,n,grad)+dmsnow/den_snow
         dmsnow=max(dmsnow,-hsnow(m,n,grad)*den_snow)
         hsnow(m,n,grad)=max(hsnew,0.0)         
         wf_tot_atm(m,n)=wf_tot_atm(m,n)+dmsnow/tau
         
! if all the snow is sublimated, remaining heat deficit sublimates ice
          dmice= min(hsnew,0.0)*den_snow
          hinew =hice(m,n,grad)+dmice/den_ice
          dmice=max(dmice,-hice(m,n,grad)*den_ice)
          hice(m,n,grad)=max(hinew,0.0) 
          wf_tot_atm(m,n)=wf_tot_atm(m,n)+dmice/tau         
          dqoc = min(hinew,0.0)*den_ice*lat_heat_subl(tmelt)
          dmoc = dqoc/lat_heat_vapor(tf)
          wf_tot_atm(m,n)= wf_tot_atm(m,n) + dmoc/tau
          wf_tot_oc(m,n) = wf_tot_oc(m,n) + dmoc/tau         
          hf_tot_oc(m,n) = hf_tot_oc(m,n) + dqoc/tau

          if(hinew<0.0) then
           sf_tot_oc(m,n)= sf_tot_oc(m,n) + sice(m,n,grad)*den_ice/tau
           sice(m,n,grad)=0.0
          endif
 
         if(melt==1) then

 ! snow melting
          dmsnow = - max(hbal,0.0)*aice(m,n,grad)/lambda_f*tau
          hsnew = hsnow(m,n,grad) + dmsnow/den_snow

!          if(hsnew<0) then
!           t_is=tmelt
!          endif

! if all the snow melted, remaining heat fuses underlying ice         
         dmsnow=max(dmsnow,-hsnow(m,n,grad)*den_snow)         
         hsnow(m,n,grad)=max(hsnew,0.0)
         heat_top=-min(hsnew,0.0)*den_snow*lambda_f
         wf_tot_oc(m,n)=wf_tot_oc(m,n)-dmsnow/tau
! end of snow melting
        endif

! ------------------------------------------------------------------------------       
       else        !if ice is not covered with snow

        melt=0
         
         t1=-70.0
         t2= 1.0

         hbal= 1000.0
         iter=0

        do while(abs(t2-t1)>accur_tem.and.abs(hbal)>accur_flux.and.iter<=30)
         
         iter=iter+1
         tmid=(t1+t2)/2.0
         t_is=tmid

         call  hbal_ice(hbal,   &
                wind_ai(m,n),   &
                   slpr(m,n),   &
                   tatm(m,n),   &
                   qatm(m,n),   &
                    lwr(m,n),   &
                    swr(m,n),   &
                          hi,   &
                   swpen_top,   &
                   swpen_bot,   &
                      kiceph,   &
                        t_is,   &
                          tf,   &
                          sh,   &
                          ev,   &
                          lh,   &
                        cdai,   &
                          lw,   &
                         swb,   &
                         swi,   &
                  AlbOpWater(m,n))         
        
        if(hbal>0.0) then
          t1=tmid
        else
          t2=tmid
        endif
           
        enddo
if(iter>30) write(*,'(i2,2(a,i5),4(a,f12.5), a,i3)') 2, ' m=', m, ' n=', n, ' t1=', t1, ' t2=', t2, ' t_is=', t_is, ' hbal=', hbal, ' iter=', iter
              
        if(t_is>tmelt) then
 !recomputing fluxes and ice melting         
            melt=1
            t_is=tmelt

         call hbal_ice(hbal,   &
               wind_ai(m,n),   &
                  slpr(m,n),   &
                  tatm(m,n),   &
                  qatm(m,n),   &
                   lwr(m,n),   &
                   swr(m,n),   &
                         hi,   &
                  swpen_top,   &
                  swpen_bot,   &
                     kiceph,   &
                       t_is,   &
                         tf,   &
                         sh,   &
                         ev,   &
                         lh,   &
                       cdai,   &
                         lw,   &
                        swb,   &
                        swi,   &
                  AlbOpWater(m,n))  
 !end of recomputing fluxes for ice melting  

        endif
 
         sensheat(m,n)=sensheat(m,n)+sh*aice(m,n,grad)
          latheat(m,n)= latheat(m,n)+lh*aice(m,n,grad)
           lw_bal(m,n)=  lw_bal(m,n)+lw*aice(m,n,grad)
       sw_bal_atm(m,n)=  sw_bal_atm(m,n)+swb*aice(m,n,grad)
       hf_tot_atm(m,n)=  hf_tot_atm(m,n)+(sh+lh+lw+swb)*aice(m,n,grad)
        sw_bal_oc(m,n)=   sw_bal_oc(m,n)+swi*aice(m,n,grad)
        hf_tot_oc(m,n)=   hf_tot_oc(m,n)+swi*aice(m,n,grad)
            cd_ai(m,n)=   cd_ai(m,n) + cdai*aice(m,n,grad)

!  ice sublimation/desublimation
         dmice=aice(m,n,grad)*ev*tau
         hinew=hice(m,n,grad)+dmice/den_ice
         dmice=max(dmice,-hice(m,n,grad)*den_ice)
         hice(m,n,grad)=max(hinew,0.0)
         wf_tot_atm(m,n)=wf_tot_atm(m,n)+dmice/tau
         dqoc = min(hinew,0.0)*den_ice*lat_heat_subl(tmelt)
         dmoc= dqoc/lat_heat_vapor(tf)        
         wf_tot_atm(m,n)= wf_tot_atm(m,n) + dmoc/tau         
          wf_tot_oc(m,n)= wf_tot_oc(m,n)  + dmoc/tau        
          hf_tot_oc(m,n)= hf_tot_oc(m,n)  + dqoc/tau

         if(hinew<0.0) then
          sf_tot_oc(m,n)= sf_tot_oc(m,n) + sice(m,n,grad)*den_ice/tau
          sice(m,n,grad)=0.0
         endif

 ! upper surface heat flux
          heat_top = hbal*aice(m,n,grad)*tau*float(melt)
! end of ice melting prepare (upper surface heat flux)
       
       endif           !end of if ice is covered with snow  

! flux on the bottom boundary of ice

       hi1= max(min(hi,hice_max_rough),hice_min_rough)
       zstar = z0star*hi1/hice_max_rough
      
      if(ksw_idn>0) then       
       depth = z(1)*hhq(m,n)
       log_arg = max(depth/zstar,3.0)
       cd = max( min( (vonkarman/log(log_arg))**2,cbi_max ), cbi_min )
       cd_rt = sqrt(cd)
       ustar = min( max( cd_rt*u,ustar_min),ustar_max)
       cd_io(m,n)=cd_io(m,n)+RefDen*cd_rt*ustar*aice(m,n,grad)      
      else
       ustar=ustar_ref 
      endif

       qiw=RefDen*HeatCapWater* sqrt(ustar/zstar)/b_io_hflux * max(tt(m,n,1)-tf,0.0)*aice(m,n,grad)
       heat_bot=(qiw +kiceph*(t_is-tf)*aice(m,n,grad))*tau
       hf_tot_oc(m,n)= hf_tot_oc(m,n) - qiw

       dmice= - (heat_top+heat_bot)/lambda_f
       hinew=hice(m,n,grad)+dmice/den_ice
       hiold=max(hice(m,n,grad),himin)
       dmice=max(dmice,-hice(m,n,grad)*den_ice)
       
       if(dmice>0) then  !ice freezing
         dmsice=dmice*sal_oc*freeze_frac_sal
!         daice =(1.0-aistot(m,n))*dmice/href/den_ice
!         daice =(1.0-aistot(m,n))*(1.0-exp(-dmice/den_ice/href))
          daice=0.0
         dmsnow=0.0
       else         !ice melting
         
         if(hinew<0.0) then !total melt
          dmoc=min(hinew,0.0)*den_ice          
          hf_tot_oc(m,n)= hf_tot_oc(m,n)-dmoc*lambda_f/tau 
         endif

         dmsice=dmice/hiold*sice(m,n,grad)
!lateral melting              
!        daice=cmelt*aice(m,n,grad)*dmice/hiold/den_ice        
         daice=aice(m,n,grad)*( (max(hinew,0.0)/hiold)**cmelt-1.0 )
         dmsnow= hsnow(m,n,grad)*daice/aice(m,n,grad)*den_snow

       endif
!End of melting

        wf_tot_oc(m,n)= wf_tot_oc(m,n)-dmice/tau-dmsnow/tau        
        sf_tot_oc(m,n)= sf_tot_oc(m,n)-dmsice/tau

        hsnow(m,n,grad)=hsnow(m,n,grad)+dmsnow/den_snow
        sice(m,n,grad)=sice(m,n,grad)+dmsice/den_ice
        hice(m,n,grad)=hice(m,n,grad)+dmice/den_ice       
        aice(m,n,grad)=aice(m,n,grad)+daice
 
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
      
      endif       !end of if ice is available 

     enddo        !End of the cycle over thickness gradations

!new ice formation
     
     hsug = min( hmax_sugar, 400.0/max(hitot(m,n)**2,himin**2) )
     dmi = 0.0
     dqb = 0.0
     dms = 0.0
     mark = 0

     do k=1,nz

!      depth=hhq(m,n)*zw(k+1) 
       depth=hhq_rest(m,n)*zw(k+1) 
        
       if(depth<hsug) then
        
        sal_oc=max(ss(m,n,k)+salref,0.0)
        tf = freeze_temp(sal_oc)
        vol=dz(k)*hhq(m,n)        

        if(tt(m,n,k)<tf) then

         mark = mark+1
         hbal=RefDen*HeatCapWater*(tf-tt(m,n,k))*vol
         dmice = hbal/lambda_f
         dmsice= sal_oc*freeze_frac_sal*dmice
         tt(m,n,k)=tf
         dmi = dmi + dmice
         dqb = dqb + hbal
         dms = dms +dmsice

        endif

       endif
     enddo
     
     if(mark>0) then 

       daice =(1.0-aistot(m,n))*(1.0-exp(-dmi/den_ice/href))
       
       hf_sugar(m,n)  = hf_sugar(m,n) + dqb/tau
       hf_tot_oc(m,n) = hf_tot_oc(m,n) + dqb/tau
       wf_tot_oc(m,n) = wf_tot_oc(m,n) - dmi/tau
       sf_tot_oc(m,n) = sf_tot_oc(m,n) -dms/tau

       hice(m,n,1)=hice(m,n,1)+dmi/den_ice
       aice(m,n,1)=min(aice(m,n,1)+daice,aimax)
       sice(m,n,1)=sice(m,n,1)+dms/den_ice
     endif

   endif          !end check water area
  enddo           !end loop on longitude
enddo             !end loop on latitude
!$omp end parallel do

 call syncborder_real(aice, mgrad)
 call syncborder_real(hice, mgrad)
 call syncborder_real(sice, mgrad)
 call syncborder_real(hsnow, mgrad)
 call syncborder_real(cd_ai, 1)
 call syncborder_real(cd_io, 1)
 call syncborder_real(tt, nz)

  if(periodicity_x/=0) then
   call cyclize_x( aice,nx,ny,mgrad,mmm,mm)
   call cyclize_x( hice,nx,ny,mgrad,mmm,mm)
   call cyclize_x( sice,nx,ny,mgrad,mmm,mm)
   call cyclize_x(hsnow,nx,ny,mgrad,mmm,mm)
   call cyclize_x(cd_ai,nx,ny,1,mmm,mm)
   call cyclize_x(cd_io,nx,ny,1,mmm,mm)
   call cyclize_x(tt,nx,ny,nz,mmm,mm)
  endif
  
  if(periodicity_y/=0) then
   call cyclize_y( aice,nx,ny,mgrad,nnn,nn)
   call cyclize_y( hice,nx,ny,mgrad,nnn,nn)
   call cyclize_y( sice,nx,ny,mgrad,nnn,nn)
   call cyclize_y(hsnow,nx,ny,mgrad,nnn,nn)
   call cyclize_y(cd_ai,nx,ny,1,nnn,nn)
   call cyclize_y(cd_io,nx,ny,1,nnn,nn)
   call cyclize_y(tt,nx,ny,nz,nnn,nn)
  endif
  
endsubroutine seaice_thermodynamics

!==================================================================
subroutine trans_ice_mpdata(ff,u,v,tau,uad,vad,fm,rhs)
include 'transport.fi'

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2), intent(in):: u,v
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad):: uad,vad,fm,rhs
real(4), intent(in):: tau

real(4) fx_p, fx_m, fy_p, fy_m
integer m,n,grad,iter

!Initiation of temporary arrays
!$omp parallel do private(m,n,grad) 
 do n=ny_start-1, ny_end+1
  do m=nx_start-1, nx_end+1
   do grad=1,mgrad
      uad(m,n,grad)=u(m,n)
      vad(m,n,grad)=v(m,n)
      rhs(m,n,grad)=0.0
       fm(m,n,grad)=ff(m,n,grad)
   enddo
  enddo
 enddo
!$omp end parallel do

!The basic transport scheme - upwind

!$omp parallel do private(m,n,grad, fx_p, fx_m, fy_p, fy_m) 
 do n=ny_start, ny_end
  do m=nx_start, nx_end
   if(lu(m,n)>0.5) then
     do grad=1, mgrad
      fx_p= ( u(m  ,n)+abs(u(m  ,n)) )*fm(m  ,n,grad)        &
          + ( u(m  ,n)-abs(u(m  ,n)) )*fm(m+1,n,grad)
      fx_m= ( u(m-1,n)+abs(u(m-1,n)) )*fm(m-1,n,grad)        &
          + ( u(m-1,n)-abs(u(m-1,n)) )*fm(m  ,n,grad)
      fy_p= ( v(m,n  )+abs(v(m,n  )) )*fm(m,n  ,grad)        &
          + ( v(m,n  )-abs(v(m,n  )) )*fm(m,n+1,grad)
      fy_m= ( v(m,n-1)+abs(v(m,n-1)) )*fm(m,n-1,grad)        &
          + ( v(m,n-1)-abs(v(m,n-1)) )*fm(m,n  ,grad)
      rhs(m,n,grad)= - ( fx_p*dyh(m,n)-fx_m*dyh(m-1,n)          &
                        +fy_p*dxh(m,n)-fy_m*dxh(m,n-1) )/2.0*tau
      ff(m,n,grad)=fm(m,n,grad)+rhs(m,n,grad)/sqt(m,n)
     enddo
   endif
  enddo
 enddo
!$omp end parallel do

 call syncborder_real(ff, mgrad)

 if(periodicity_x/=0) then
  call cyclize_x( ff ,nx,ny,mgrad,mmm,mm)
 endif

 if(periodicity_y/=0) then
  call cyclize_y( ff ,nx,ny,mgrad,nnn,nn)
 endif

!Correction of smoothed upwind solution by anti-diffusion fluxes
if(niter_mpdata>0) then

!iterations for reducing numerical dissipation 
 do iter=1, niter_mpdata
 
  call syncborder_real(rhs, mgrad)

 if(periodicity_x/=0) then
  call cyclize_x( rhs ,nx,ny,mgrad,mmm,mm)
 endif

 if(periodicity_y/=0) then
  call cyclize_y( rhs ,nx,ny,mgrad,nnn,nn)
 endif

 !Computing anti-diffusion velocities

   !$omp parallel do private(m,n,grad) 
    do n=ny_start, ny_end
     do m=nx_start, nx_end
   
      if(lcu(m,n)>0.5) then
        do grad=1, mgrad
         uad(m,n,grad)= ( abs(uad(m,n,grad))*(  ff(m+1,n,grad)- ff(m,n,grad) )      &
                            + uad(m,n,grad) *( rhs(m+1,n,grad)+rhs(m,n,grad) )/2.0/squ(m,n) )/max(ff(m+1,n,grad)+ff(m,n,grad),2e-10) 
        enddo
      endif
   
      if(lcv(m,n)>0.5) then
        do grad=1, mgrad
         vad(m,n,grad)= ( abs(vad(m,n,grad))*(  ff(m,n+1,grad)- ff(m,n,grad) )      &
                            + vad(m,n,grad) *( rhs(m,n+1,grad)+rhs(m,n,grad) )/2.0/sqv(m,n) )/max(ff(m,n+1,grad)+ff(m,n,grad),2e-10) 
        enddo
      endif 
   
     enddo
    enddo
   !$omp end parallel do 

   call syncborder_real(uad,mgrad)
   call syncborder_real(vad,mgrad)

   if(periodicity_x/=0) then
    call cyclize_x(uad,nx,ny,mgrad,mmm,mm)
   endif

   if(periodicity_y/=0) then
    call cyclize_y(vad,nx,ny,mgrad,nnn,nn)
   endif
 

!Initial values of the tracer
!$omp parallel do private(m,n,grad) 
 do n=ny_start-1, ny_end+1
  do m=nx_start-1, nx_end+1
   do grad=1,mgrad
     fm(m,n,grad)=ff(m,n,grad)
   enddo
  enddo
 enddo
!$omp end parallel do
 !The upwind correcting scheme with antidiffusion velocities

!$omp parallel do private(m,n,grad, fx_p, fx_m, fy_p, fy_m) 
 do n=ny_start, ny_end
  do m=nx_start, nx_end
   if(lu(m,n)>0.5) then
     do grad=1, mgrad
      fx_p= ( uad(m  ,n,grad)+abs(uad(m  ,n,grad)) )*fm(m  ,n,grad)        &
          + ( uad(m  ,n,grad)-abs(uad(m  ,n,grad)) )*fm(m+1,n,grad)
      fx_m= ( uad(m-1,n,grad)+abs(uad(m-1,n,grad)) )*fm(m-1,n,grad)        &
          + ( uad(m-1,n,grad)-abs(uad(m-1,n,grad)) )*fm(m  ,n,grad)
      fy_p= ( vad(m,n  ,grad)+abs(vad(m,n  ,grad)) )*fm(m,n  ,grad)        &
          + ( vad(m,n  ,grad)-abs(vad(m,n  ,grad)) )*fm(m,n+1,grad)
      fy_m= ( vad(m,n-1,grad)+abs(vad(m,n-1,grad)) )*fm(m,n-1,grad)        &
          + ( vad(m,n-1,grad)-abs(vad(m,n-1,grad)) )*fm(m,n  ,grad)
      rhs(m,n,grad)= - (fx_p*dyh(m,n)-fx_m*dyh(m-1,n)                      &
                       +fy_p*dxh(m,n)-fy_m*dxh(m,n-1) )/2.0*tau
      ff(m,n,grad)=fm(m,n,grad)+rhs(m,n,grad)/sqt(m,n)
     enddo
   endif
  enddo
 enddo
!$omp end parallel do

 call syncborder_real(ff, mgrad)

  if(periodicity_x/=0) then
   call cyclize_x( ff,nx,ny,mgrad,mmm,mm)
  endif
  
  if(periodicity_y/=0) then
   call cyclize_y( ff,nx,ny,mgrad,nnn,nn)
  endif
  
 enddo

endif

endsubroutine trans_ice_mpdata

!-------------------------------------------------------------------------------------------
subroutine seaice_dynamics(tau, nicestep)
use ocean_variables

real(4), intent(in):: tau
integer, intent(in):: nicestep

real(4) tau_inner, tdamp, coef1, coef2, coef0, denom
real(4) delta, epsd2, epst2, epss2, pres
real(4) drag_ai, drag_io, cdai, cdio, wnd, atot
real(4) rhs
integer m,n,k,l

tau_inner=tau/float(nicestep)
tdamp=tau/3.0
coef0=float(nicestep)
coef1 = tau_inner/(2.0*tdamp)
coef2 = coef1/extr2
denom = 1.0 + coef1

!computing needed values
        
!$omp parallel do private(m,n)
do n=ny_start-1, ny_end+1
 do m=nx_start-1, nx_end+1
   if(lu(m,n)>0.5) then
             pice(m,n) = Pice_cr*hitot(m,n)/max(aistot(m,n),aimin)*exp(-Cstar*(1.0-aistot(m,n)))   
       idyn_coefx(m,n) = 0.0
       idyn_coefy(m,n) = 0.0
   idyn_coefx_tot(m,n) = 0.0
   idyn_coefy_tot(m,n) = 0.0
        idyn_rhsx(m,n) = 0.0
        idyn_rhsy(m,n) = 0.0
    idyn_rhsx_tot(m,n) = 0.0
    idyn_rhsy_tot(m,n) = 0.0
          mistotu(m,n) = 0.0
          mistotv(m,n) = 0.0           
            piceh(m,n) = 0.0
            sig_d(m,n) = 0.0
            sig_t(m,n) = 0.0
            sig_s(m,n) = 0.0
   endif 
 enddo
enddo
!$omp end parallel do


!$omp parallel do private(m,n,drag_ai,drag_io,cdai, cdio, wnd)
do n=ny_start, ny_end
 do m=nx_start, nx_end
 
   if(lcu(m,n)>0.5) then
    
    mistotu(m,n)=max((mistot(m,n)+mistot(m+1,n))/2.0, ice_mass_min)

    cdai = (cd_ai(m,n) + cd_ai(m+1,n))/2.0 
    wnd = max( (wind_ai(m,n)+wind_ai(m+1,n))/2.0, 1.0 )
    cdai= max( cdai, aimin_dyn*airden_ref*cda_ref*wnd )
    drag_ai = cdai * windx(m,n)

    cdio = (cd_io(m,n) + cd_io(m+1,n))/2.0
    cdio = max( cdio, aimin_dyn*RefDen*sqrt(cdio_ref)*ustar_ref)
    drag_io = cdio * uu(m,n,1)  
    
    idyn_coefx(m,n) = 1.0/tau * (1.0+coef0) + (cdai + cdio)/mistotu(m,n)
    idyn_coefx_tot(m,n) = idyn_coefx(m,n)**2 + (rlh_s(m,n)**2+rlh_s(m,n-1)**2)/2.0                                        

    idyn_rhsx(m,n)=uice(m,n)/tau                                                      &
              - ( FreeFallAcc * (   ssh(m+1,n)-ssh(m,n)                               &
              +  (mistot(m+1,n)-mistot(m,n))/RefDen*float(variable_volume_budget) )   &
              +    (slpr(m+1,n)-slpr(m,n)  )/RefDen )/dxt(m,n)                        &
              + (drag_ai+drag_io)/mistotu(m,n)
   endif

   if(lcv(m,n)>0.5) then

    mistotv(m,n)=max((mistot(m,n)+mistot(m,n+1))/2.0, ice_mass_min)

    cdai = (cd_ai(m,n) + cd_ai(m,n+1))/2.0 
    wnd = max( (wind_ai(m,n)+wind_ai(m,n+1))/2.0, 1.0 )
    cdai= max( cdai, aimin_dyn*airden_ref*cda_ref*wnd )
    drag_ai = cdai * windy(m,n)

    cdio = (cd_io(m,n) + cd_io(m,n+1))/2.0
    cdio = max( cdio, aimin_dyn*RefDen*sqrt(cdio_ref)*ustar_ref)
    drag_io = cdio * vv(m,n,1)
    
    idyn_coefy(m,n) = 1.0/tau * (1.0+coef0) + (cdai + cdio)/mistotv(m,n)
    idyn_coefy_tot(m,n)=  idyn_coefy(m,n)**2 + (rlh_s(m,n)**2+rlh_s(m-1,n)**2)/2.0                                        

    idyn_rhsy(m,n)= vice(m,n)/tau                                                           &
               - ( FreeFallAcc * (   ssh(m,n+1)-ssh(m,n)                               &
               +  (mistot(m,n+1)-mistot(m,n))/RefDen*float(variable_volume_budget) )   &
               +    (slpr(m,n+1)-slpr(m,n)  )/RefDen )/dyt(m,n)                        &
               + (drag_ai+drag_io)/mistotv(m,n)
   endif
   
   if(luu(m,n)>0.5) then
    piceh(m,n) = ( pice(m,n) + pice(m,n+1) + pice(m+1,n) + pice(m+1,n+1))/4.0
   endif
    
 enddo
enddo
!$omp end parallel do 


do l=1,nicestep
!computing strain rate tensor components
!$omp parallel do private(m,n)
do n=ny_start,ny_end
 do m=nx_start,nx_end
    
    if(lu(m,n)>0.5) then
     eps_d(m,n)=( ( uice(m,n)*dyh(m,n) - uice(m-1,n)*dyh(m-1,n) )                   &
                 +( vice(m,n)*dxh(m,n) - vice(m,n-1)*dxh(m,n-1) ) )/sqt(m,n) 
 
     eps_t(m,n)=( ( uice(m,n)/dyh(m,n) - uice(m-1,n)/dyh(m-1,n) )*dy(m,n)**2        &
                 -( vice(m,n)/dxh(m,n) - vice(m,n-1)/dxh(m,n-1) )*dx(m,n)**2 )/sqt(m,n)
    endif
  
    if(luu(m,n)>0.5) then
     eps_s(m,n)=( ( uice(m,n+1)/dxt(m,n+1) - uice(m,n)/dxt(m,n) )*dxb(m,n)**2      &
                 +( vice(m+1,n)/dyt(m+1,n) - vice(m,n)/dyt(m,n) )*dyb(m,n)**2 )/sqh(m,n)
    endif

 enddo
enddo
!$omp end parallel do

call syncborder_real(eps_d, 1)
call syncborder_real(eps_t, 1)
call syncborder_real(eps_s, 1)

 if(periodicity_x/=0) then
  call cyclize_x(eps_d,nx,ny,1,mmm,mm)
  call cyclize_x(eps_t,nx,ny,1,mmm,mm)
  call cyclize_x(eps_s,nx,ny,1,mmm,mm)
 endif

 if(periodicity_y/=0) then
  call cyclize_y(eps_d,nx,ny,1,nnn,nn)
  call cyclize_y(eps_t,nx,ny,1,nnn,nn)
  call cyclize_y(eps_s,nx,ny,1,nnn,nn)
 endif
 
!computing stress tensor components
!$omp parallel do private(m,n,delta, epsd2, epst2, epss2, pres)
do n=ny_start,ny_end
 do m=nx_start,nx_end

    if(lu(m,n)>0.5) then
      epss2 = ( eps_s(m,n)**2 + eps_s(m,n-1)**2 + eps_s(m-1,n)**2 + eps_s(m-1,n-1)**2 )/4.0
      delta = sqrt( eps_d(m,n)**2  +(eps_t(m,n)**2+epss2)/extr2 )
!     pres = min(pice(m,n), ice_mass_damp*sqt(m,n)*tdamp/tau_inner**2*delta)
      pres = pice(m,n)            
      sig_d(m,n)=( sig_d(m,n) + pres*( eps_d(m,n)/(delta+def_rate_min) - 1.0 )*coef1  )/denom
      sig_t(m,n)=( sig_t(m,n) + pres * eps_t(m,n)/(delta+def_rate_min)        *coef2  )/denom
   
    endif
  
    if(luu(m,n)>0.5) then

      epsd2 = ( eps_d(m,n)**2 + eps_d(m,n+1)**2 + eps_d(m+1,n)**2 + eps_d(m+1,n+1)**2 )/4.0
      epst2 = ( eps_t(m,n)**2 + eps_t(m,n+1)**2 + eps_t(m+1,n)**2 + eps_t(m+1,n+1)**2 )/4.0
      
      delta = sqrt( epsd2 + (epst2 + eps_s(m,n)**2)/extr2 )
!     pres = min(piceh(m,n), ice_mass_damp*sqh(m,n)*tdamp/tau_inner**2*delta)
      pres = pice(m,n)
         
      sig_s(m,n)=( sig_s(m,n) + pres * eps_s(m,n)/(delta+def_rate_min)        *coef2 )/denom
      
    endif
 
 enddo
enddo
!$omp end parallel do

call syncborder_real(sig_d, 1)
call syncborder_real(sig_t, 1)
call syncborder_real(sig_s, 1)

 if(periodicity_x/=0) then
  call cyclize_x(sig_d,nx,ny,1,mmm,mm)
  call cyclize_x(sig_t,nx,ny,1,mmm,mm)
  call cyclize_x(sig_s,nx,ny,1,mmm,mm)
 endif

 if(periodicity_y/=0) then
  call cyclize_y(sig_d,nx,ny,1,nnn,nn)
  call cyclize_y(sig_t,nx,ny,1,nnn,nn)
  call cyclize_y(sig_s,nx,ny,1,nnn,nn)
 endif

!RHS terms
!$omp parallel do private(m,n,rhs)
do n=ny_start, ny_end
 do m=nx_start, nx_end
 
   if(lcu(m,n)>0.5) then
     
     rhs =  ( dy(m+1,n)**2*sig_t(m+1,n) - dy(m,n)**2*sig_t(m,n) )/(dyh(m,n)*squ(m,n))      &
         + (dxb(m,n)**2*sig_s(m,n) - dxb(m,n-1)**2*sig_s(m,n-1) )/(dxt(m,n)*squ(m,n))      &
         + (sig_d(m+1,n)-sig_d(m,n))/dxt(m,n) 

     idyn_rhsx_tot(m,n)=idyn_rhsx(m,n) + uice(m,n)/tau*coef0 + rhs/2.0/mistotu(m,n)

   endif

   if(lcv(m,n)>0.5) then
     
     rhs = -( dx(m,n+1)**2*sig_t(m,n+1) - dx(m,n)**2*sig_t(m,n) )/(dxh(m,n)*sqv(m,n))      &
         + (dyb(m,n)**2*sig_s(m,n) - dyb(m-1,n)**2*sig_s(m-1,n) )/(dyt(m,n)*sqv(m,n))      &
         + (sig_d(m,n+1)-sig_d(m,n))/dyt(m,n)
     
     idyn_rhsy_tot(m,n)=idyn_rhsy(m,n) + vice(m,n)/tau*coef0 + rhs/2.0/mistotv(m,n)

   endif
    
 enddo
enddo
!$omp end parallel do 

     call syncborder_real(idyn_rhsx_tot, 1)
     call syncborder_real(idyn_rhsy_tot, 1)

      if(periodicity_x/=0) then
       call cyclize_x(idyn_rhsx_tot,nx,ny,1,mmm,mm)
       call cyclize_x(idyn_rhsy_tot,nx,ny,1,mmm,mm)
      endif
 
      if(periodicity_y/=0) then
       call cyclize_y(idyn_rhsx_tot,nx,ny,1,nnn,nn)
       call cyclize_y(idyn_rhsy_tot,nx,ny,1,nnn,nn)
      endif

!Velocity computation 
!$omp parallel do private(m,n)
do n=ny_start, ny_end
 do m=nx_start, nx_end
 
   if(lcu(m,n)>0.5) then
    uice(m,n)=(    idyn_coefx(m,n)*idyn_rhsx_tot(m,n)                               &
              +  ( rlh_s(m,n  )*(idyn_rhsy_tot(m  ,n  )+idyn_rhsy_tot(m+1,n  ))     &
                 + rlh_s(m,n-1)*(idyn_rhsy_tot(m  ,n-1)+idyn_rhsy_tot(m+1,n-1)) )/4.0   )/idyn_coefx_tot(m,n)
   endif

   if(lcv(m,n)>0.5) then
    vice(m,n)=(    idyn_coefy(m,n)*idyn_rhsy_tot(m,n)                               &
              -  ( rlh_s(m  ,n)*(idyn_rhsx_tot(m  ,n  )+idyn_rhsx_tot(m  ,n+1))     &
                 + rlh_s(m-1,n)*(idyn_rhsx_tot(m-1,n  )+idyn_rhsx_tot(m-1,n+1)) )/4.0  )/idyn_coefy_tot(m,n)
   endif
    
 enddo
enddo
!$omp end parallel do 
     
     call syncborder_real(uice, 1)
     call syncborder_real(vice, 1)

      if(periodicity_x/=0) then
       call cyclize_x(uice,nx,ny,1,mmm,mm)
       call cyclize_x(vice,nx,ny,1,mmm,mm)
      endif
 
      if(periodicity_y/=0) then
       call cyclize_y(uice,nx,ny,1,nnn,nn)
       call cyclize_y(vice,nx,ny,1,nnn,nn)
      endif

enddo !end of internal time cycle

!$omp parallel do private(m,n,atot,cdai,cdio)
do n=ny_start, ny_end
 do m=nx_start, nx_end
   
   if(lcu(m,n)>0.5) then
    atot=(aistot(m,n)+aistot(m+1,n))/2.0
    cdai = (cd_ai(m,n) + cd_ai(m+1,n))/2.0 
    cdio = (cd_io(m,n) + cd_io(m+1,n))/2.0

    taux_atm(m,n)= taux_atm(m,n)*(1.0-atot)                      &
                  + cdai*(windx(m,n)-uice(m,n))

    taux_oc(m,n)= taux_oc(m,n)*(1.0-atot)                        &
                  + cdio*(uice(m,n)-uu(m,n,1))

     tx_dif(m,n)= tx_dif(m,n)*(1.0-atot)                         &
                  + cdio*uice(m,n)

    tx_coef(m,n)=tx_coef(m,n)*(1.0-atot)                         &
                  + cdio   
   endif
   
   if(lcv(m,n)>0.5) then
    atot=(aistot(m,n)+aistot(m,n+1))/2.0
    cdai = (cd_ai(m,n) + cd_ai(m,n+1))/2.0 
    cdio = (cd_io(m,n) + cd_io(m,n+1))/2.0

    tauy_atm(m,n)= tauy_atm(m,n)*(1.0-atot)                      &
                 + cdai*(windy(m,n)-vice(m,n))
    
    tauy_oc(m,n)= tauy_oc(m,n)*(1.0-atot)                        &
                 + cdio*(vice(m,n)-vv(m,n,1))
    
     ty_dif(m,n)= ty_dif(m,n)*(1.0-atot)                         &
                 + cdio*vice(m,n)

    ty_coef(m,n)=ty_coef(m,n)*(1.0-atot)                         &
                 + cdio  
   endif   
 
 enddo
enddo
!$omp end parallel do 

endsubroutine seaice_dynamics

endmodule seaice_routes
