module mixing_routes
use constants
use basin_grid
use ocalg_routes
implicit none

contains

!========================================================
subroutine vertical_frequencies(u,v,t,s,den_s,k_s,vbf2,shf2)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in) :: den_s,k_s,u,v,t,s
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: vbf2,shf2
 
integer m,n,k
 
real(4) u1,u2,v1,v2,r1,r2,pgrz,dns

!$omp parallel do private(m,n,k,u1,v1,u2,v2,r1,r2, pgrz,dns)
      do n=ny_start,ny_end
       do m=nx_start,nx_end

        if (lu(m,n)>0.5) then

           pgrz=RefDen*FreeFallAcc*hhq(m,n)
                       
           u2 = (u(m  ,n,1)*dyh(m  ,n)*hhu(m  ,n)    &
               + u(m-1,n,1)*dyh(m-1,n)*hhu(m-1,n))/2.0/hhq(m,n)/dy(m,n)       !u on t-grid
           v2 = (v(m,n  ,1)*hhv(m,n  )*dxh(m,n  )    &
               + v(m,n-1,1)*hhv(m,n-1)*dxh(m,n-1))/2.0/hhq(m,n)/dx(m,n)       !v on t-grid
                       
          do k=2,nz
                      
           r1 = density_unesco_tpot(t(m,n,k-1),s(m,n,k-1)+salref,       &
                                pgrz*zw(k),den_s(m,n,k-1),k_s(m,n,k-1))  
           r2 = density_unesco_tpot(t(m,n,k  ),s(m,n,k  )+salref,       &
                                pgrz*zw(k),den_s(m,n,k  ),k_s(m,n,k  )) 
           dns= 1000.0 + (r1+r2)/2.0

           u1 = u2
           v1 = v2
           
           u2 = (u(m  ,n,k)*dyh(m  ,n)*hhu(m  ,n)    &
               + u(m-1,n,k)*dyh(m-1,n)*hhu(m-1,n))/2.0/hhq(m,n)/dy(m,n)       !u on t-grid
           v2 = (v(m,n  ,k)*hhv(m,n  )*dxh(m,n  )    &
               + v(m,n-1,k)*hhv(m,n-1)*dxh(m,n-1))/2.0/hhq(m,n)/dx(m,n)       !v on t-grid

            vbf2(m,n,k) = FreeFallAcc/dns*(r2-r1)/(hzt(k)*hhq(m,n))
            shf2(m,n,k) = max((u2-u1)**2+(v2-v1)**2,0.0025e-4)/(hzt(k)*hhq(m,n))**2
    
          enddo
        endif

       enddo
      enddo
!$omp end parallel do

endsubroutine vertical_frequencies

!========================================================
subroutine ppmix(vbf2,shf2,anzu,anumaxu,anubgru,anzt,anumaxt,anubgrt)

real(4), intent(in):: anumaxu,anubgru,anumaxt,anubgrt
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in) :: vbf2,shf2
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: anzu,anzt
   
integer m,n,k
real(4) undimdepth, rich, nu_min_t, nu_min_u

!$omp parallel do private(m,n,k,undimdepth,rich,nu_min_t,nu_min_u)
      do n=ny_start,ny_end
       do m=nx_start,nx_end

        if (lu(m,n)>0.5) then
          undimdepth = dimdepth/hhq(m,n)
          do k=2,nz

           if(zw(k)<undimdepth) then
            nu_min_t = anubgrt + anumaxt/10.0
            nu_min_u = anubgru + anumaxu/10.0
           else
            nu_min_t = anubgrt
            nu_min_u = anubgru
           endif

 ! direct discrete form of Ri           
           rich = max( min( vbf2(m,n,k)/shf2(m,n,k),1e+05 ), 0.0)
           anzu(m,n,k) = min( anumaxu/(1.0+5.0*rich)**2 + nu_min_u, anumaxu)
           anzt(m,n,k) = min( anzu(m,n,k)/(1.0+5.0*rich)+ nu_min_t, anumaxt)
          
          enddo
        endif

       enddo
      enddo
!$omp end parallel do

      call syncborder_real(anzu, nz+1)

      if(periodicity_x/=0) then 
        call cyclize_x(anzu,nx,ny,nz+1,mmm,mm)
      endif

      if(periodicity_y/=0) then 
        call cyclize_y(anzu,nx,ny,nz+1,nnn,nn)
      endif

endsubroutine ppmix

!====================================================================================
subroutine stress_tension(u,v,str_t)
 
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in) :: u,v
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: str_t 
 
integer m,n,k

 !$omp parallel do private(m,n,k) 
 do n=ny_start, ny_end
   do m=nx_start, nx_end

    if(lu(m,n)>0.5) then
     do k=1,nz 
      str_t(m,n,k)=dy(m,n)/dx(m,n)*(u(m,n,k)/dyh(m,n)-u(m-1,n,k)/dyh(m-1,n))     &
                  -dx(m,n)/dy(m,n)*(v(m,n,k)/dxh(m,n)-v(m,n-1,k)/dxh(m,n-1)) 
     enddo
    endif

   enddo
 enddo
!$omp end parallel do

    call syncborder_real(str_t, nz)

    if(periodicity_x/=0) then
       call cyclize_x(str_t,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
       call cyclize_y(str_t,nx,ny,nz,nnn,nn)
    end if

endsubroutine stress_tension
!====================================================================================
subroutine stress_shear(u,v,str_s)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in) :: u,v
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: str_s
 
integer m,n,k

 !$omp parallel do private(m,n,k) 
 do n=ny_start, ny_end
   do m=nx_start, nx_end

    if(luu(m,n)>0.5) then
     do k=1,nz
      str_s(m,n,k)=dxb(m,n)/dyb(m,n)*(u(m,n+1,k)/dxt(m,n+1)-u(m,n,k)/dxt(m,n))    &
                  +dyb(m,n)/dxb(m,n)*(v(m+1,n,k)/dyt(m+1,n)-v(m,n,k)/dyt(m,n))    
     enddo
    endif

   enddo
 enddo
!$omp end parallel do

    call syncborder_real(str_s, nz)

    if(periodicity_x/=0) then
       call cyclize_x(str_s,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
       call cyclize_y(str_s,nx,ny,nz,nnn,nn)
    end if

endsubroutine stress_shear

!======================================================================
subroutine lateral_mix_coeff(ldiff_ts_upw, lvisc_upw, lvisc_4_upw,      &
                             ldiff_ts_smg, lvisc_smg, lvisc_4_smg,      &
                           str_t,str_s,amts_ref,amuv_ref,amuv4_ref,     &
                           amts,amuv,amuv4,u,v)

real(4), intent(in):: ldiff_ts_upw, lvisc_upw, lvisc_4_upw
real(4), intent(in):: ldiff_ts_smg, lvisc_smg, lvisc_4_smg
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in) :: amts_ref,amuv_ref,amuv4_ref,u,v
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(inout):: amts,amuv,amuv4,str_t,str_s

real(4) str_s2, mu_smag, mu4_smag, mu_upw, mu4_upw,u2,dr2
integer m,n,k

! Computing stress components
 call stress_tension(u,v,str_t)
 call stress_shear  (u,v,str_s)
 
 !$omp parallel do private(m,n,k,str_s2,mu_smag,mu4_smag,mu_upw,mu4_upw,u2,dr2) 
 do n=ny_start, ny_end
   do m=nx_start, nx_end

    if(lu(m,n)>0.5) then
     dr2=dx(m,n)**2+dy(m,n)**2 
     do k=1,nz
      u2 = (u(m,n,k)**2 + u(m-1,n,k)**2 + v(m,n,k)**2 + v(m,n-1,k)**2)/2.0
      mu_upw  = sqt(m,n)*sqrt(u2/dr2)/8.0
      mu4_upw = mu_upw*sqt(m,n)/8.0
      str_s2= ( str_s(m  ,n  ,k)**2 +str_s(m-1,n  ,k)**2    &
               +str_s(m  ,n-1,k)**2 +str_s(m-1,n-1,k)**2 )/4.0    
      mu_smag=sqt(m,n)*sqrt(str_t(m,n,k)**2 + str_s2)/8.0
      mu4_smag = mu_smag*sqt(m,n)/8.0
      amuv(m,n,k) = amuv_ref(m,n,k) + mu_upw*lvisc_upw    + mu_smag*lvisc_smg
      amts(m,n,k) = amts_ref(m,n,k) + mu_upw*ldiff_ts_upw + mu_smag*ldiff_ts_smg
      amuv4(m,n,k)=amuv4_ref(m,n,k) + sqrt(mu4_upw*lvisc_4_upw) + sqrt(mu4_smag*lvisc_4_smg)
     enddo
    endif

   enddo
 enddo
!$omp end parallel do

    call syncborder_real(amuv, nz)
    call syncborder_real(amuv4,nz)
    call syncborder_real(amts, nz)

    if(periodicity_x/=0) then
       call cyclize_x(amuv ,nx,ny,nz,mmm,mm)
       call cyclize_x(amuv4,nx,ny,nz,mmm,mm)
       call cyclize_x(amts ,nx,ny,nz,mmm,mm)
    end if

    if(periodicity_y/=0) then
       call cyclize_y(amuv ,nx,ny,nz,nnn,nn)
       call cyclize_y(amuv4,nx,ny,nz,nnn,nn)
       call cyclize_y(amts ,nx,ny,nz,nnn,nn)
    end if

endsubroutine lateral_mix_coeff

!========================================================================
subroutine Mellor_Yamada_gendis(q2,q2l,rhs_q2,rhs_q2l,coef_q2,coef_q2l,vbf2,shf2,    &
                               anzt,anzu,anubgru,anubgrt,anumaxu,anumaxt,botfr)

real(4), intent(in):: anumaxu,anubgru,anumaxt,anubgrt
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in) :: vbf2,shf2
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: anzu,anzt,rhs_q2,rhs_q2l,coef_q2,coef_q2l
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: q2, q2l
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2     ), intent(in):: botfr

real(4) q, lscl, rich, gdf
real(4) stab_t, stab_u, nu_min_t, nu_min_u
real(4) wall, lm1, dpth, l_max
integer m, n, k

!$omp parallel do private(m,n,k,q,l_max,lscl,rich,     &
!$omp         stab_t,stab_u,wall,lm1,dpth,gdf,nu_min_t,nu_min_u) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    
     
     if(lu(m,n)>0.5) then
                           
      do k=2,nz
        
        dpth=zw(k)*hhq(m,n)
        lm1=1.0/dpth + 1.0/(hhq(m,n)-dpth)        
              
        q2(M,N,K) = max( q2(M,N,K),tur_var_min)
       q2l(M,N,K) = max(q2l(M,N,K),tur_var_min)  
                
        q=sqrt(q2(m,n,k))
        lscl=q2l(m,n,k)/q2(m,n,k)
        
        if(vbf2(m,n,k)>0.0) then
         l_max= 0.53 * q/sqrt(vbf2(m,n,k))
         lscl=min(lscl,l_max)
        endif
        
        rich=min(-vbf2(m,n,k)*lscl**2/q2(m,n,k),0.028)
        stab_t=A2_t*(1.0-6.0*A1_t/B1_t)     &
               / (1.0 - (3.0*A2_t*B2_t+18.0*A1_t*A2_t)*rich)
        stab_u=( A1_t*(1.0-3.0*C1_t-6.0*A1_t/B1_t)          &
                +stab_t*((18.0*A1_t**2+9.0*A1_t*A2_t)*rich) ) &
                /(1.0-9.0*A1_t*A2_t*rich)

        if(dpth<=dimdepth) then
         nu_min_t = anubgrt + anumaxt/10.0
         nu_min_u = anubgru + anumaxu/10.0
        else
         nu_min_t = anubgrt
         nu_min_u = anubgru
        endif

        anzt(m,n,k)=min(nu_min_t+q*lscl*stab_t,anumaxt)
        anzu(m,n,k)=min(nu_min_u+q*lscl*stab_u,anumaxu)
        
        gdf=anzu(m,n,k)*shf2(m,n,k) - anzt(m,n,k)*vbf2(m,n,k)
        wall=1.0+E2_t*(lscl*lm1/vonkarman)**2

          rhs_q2(m,n,k) = 2.0*gdf
         coef_q2(m,n,k) = 2.0*q/B1_t/lscl        
         rhs_q2l(m,n,k) = gdf*E1_t*lscl
        coef_q2l(m,n,k) = wall*q/B1_t/lscl   

      enddo
       ! anzu(m,n,nz+1) = botfr(m,n)*hzt(nz+1)*hhq(m,n)
     endif

    enddo
   enddo
!$omp end parallel do        

    call syncborder_real(anzu, nz+1)
    
    if(periodicity_x/=0) then 
      call cyclize_x(anzu,nx,ny,nz+1,mmm,mm)
    endif

    if(periodicity_y/=0) then 
      call cyclize_y(anzu,nx,ny,nz+1,nnn,nn)
    endif

endsubroutine Mellor_Yamada_gendis

!=====================================================================
subroutine mld_density(den,mld)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in) :: den
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2   ), intent(inout):: mld

integer m,n,k
real(4) dro

!$omp parallel do private(m,n,k,dro)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
          k=2
          dro= den(m,n,k) - den(m,n,1)
          do while(dro<dro_mld.and.k<nz)
           k=k+1
           dro= den(m,n,k) - den(m,n,1)
          enddo
          mld(m,n) = max(hhq(m,n)*zw(k),dimdepth)
        endif
       enddo
      enddo
!$omp end parallel do

endsubroutine mld_density

!====================================================================
subroutine turb_vdiff(ff,tau,nu,factor_nu,rhs,coef)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: rhs,coef,nu
real(4), intent(in):: tau, factor_nu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: ff

integer m,n,k
real(8) a(nz+1),b(nz+1),c(nz+1),eta(nz+1),rksi(nz+1) 
real(8) bp, dp, dm

!$omp parallel do private(m,n,k,bp,dp,dm,a,b,c,eta,rksi)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
         bp = hhq_n(m,n)/tau

!Upper layer          
          k=2
          
          dm = (nu(m,n,k)+nu(m,n,k-1))/2.0*factor_nu/hhq_n(m,n)/dz(k-1)
          dp = (nu(m,n,k)+nu(m,n,k+1))/2.0*factor_nu/hhq_n(m,n)/dz(k  )

          c(k) = -dp
          a(k) = 0.0
          b(k) =  bp*hzt(k) + dp + dm + coef(m,n,k)*hhq_n(m,n)*hzt(k)
          eta(k) = bp*ff(m,n,k)*hzt(k) + dm*ff(m,n,k-1) + rhs(m,n,k)*hhq_n(m,n)*hzt(k)

! internal points.
         do k=3,nz-1

          dm = dp
          dp = (nu(m,n,k)+nu(m,n,k+1))/2.0*factor_nu/hhq_n(m,n)/dz(k  )
          
          c(k) = -dp
          a(k) = -dm
          b(k) =  bp*hzt(k) + dp + dm + coef(m,n,k)*hhq_n(m,n)*hzt(k)
          eta(k) = bp*ff(m,n,k)*hzt(k) + rhs(m,n,k)*hhq_n(m,n)*hzt(k)
         enddo

!Bottom layer
          k=nz

          dm = dp
          dp = (nu(m,n,k)+nu(m,n,k+1))/2.0*factor_nu/hhq_n(m,n)/dz(k  )
          
          c(k) = 0.0
          a(k) = -dm 
          b(k) =  bp*hzt(k) + dp + dm + coef(m,n,k)*hhq_n(m,n)*hzt(k)
          eta(k) = bp*ff(m,n,k)*hzt(k) + dp*ff(m,n,k+1) + rhs(m,n,k)*hhq_n(m,n)*hzt(k)

         call factor8(nz+1,a,b,c,eta,rksi,2,nz)
         
         do k=2,nz
          ff(m,n,k)=rksi(k)
         enddo

        endif
       enddo
      enddo
!$omp end parallel do

    call syncborder_real(ff, nz+1)

    if(periodicity_x/=0) then
      call cyclize_x(ff,nx,ny,nz+1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(ff,nx,ny,nz+1,nnn,nn)            
    end if

endsubroutine turb_vdiff

!====================================================================
subroutine turb_transport_cab_mod(ff,tau,fx,fx0,fy,fy0,fz,fz0,uh,vh,ww)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ), intent(in):: uh,vh
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(in):: ww
real(4), intent(in):: tau
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1):: fx,fx0,fy,fy0
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ):: fz,fz0

integer m,n,k, ind1,ind2
real(8) hmid, rhs,  up, um, vp, vm, flxz(nz)
real(4) u, v, w
real(8) cvert(0:nz+1), flx1(0:nz+1)

cvert(0) = 0.0
cvert(1:nz) = 1.0
cvert(nz+1) = 0.0

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      fx(m,n,:) = 0.0
      fy(m,n,:) = 0.0
      fz(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do  

!Predictor step

!$omp parallel do private(m,n,k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end
      
      if (lcu(m,n)>0.5) then
       do k=2,nz
        fx(m,n,k) = (ff(m,n,k)+ ff(m+1,n,k))/2.0
       enddo
      endif

      if (lcv(m,n)>0.5) then
       do k=2,nz
        fy(m,n,k) = (ff(m,n,k) + ff(m,n+1,k))/2.0
       enddo
      endif

      if (lu(m,n)>0.5) then
       do k=1,nz
        fz(m,n,k) = (ff(m,n,k) + ff(m,n,k+1))/2.0
       enddo
      endif
                
    enddo
   enddo
!$omp end parallel do  

    call syncborder_real(fx, nz+1)
    call syncborder_real(fy, nz+1)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz+1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz+1,nnn,nn)            
    end if

!$omp parallel do private(m,n,k,hmid,rhs,flxz,up,um,vp,vm) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lu(m,n)>0.5) then
       hmid=(hhq(m,n)+hhq_n(m,n))/2.0
       
       do k=1,nz
       flxz(k) = (ww(m,n,k)+ww(m,n,k+1))/2.0*fz(m,n,k)
       enddo

       do k=2,nz
        up=(uh(m  ,n,k-1)*dz(k-1)+uh(m  ,n,k)*dz(k))/2.0
        um=(uh(m-1,n,k-1)*dz(k-1)+uh(m-1,n,k)*dz(k))/2.0
        vp=(vh(m,n  ,k-1)*dz(k-1)+vh(m,n  ,k)*dz(k))/2.0
        vm=(vh(m,n-1,k-1)*dz(k-1)+vh(m,n-1,k)*dz(k))/2.0

        rhs=(up*dyh(m,n)*fx(m,n,k) - um*dyh(m-1,n)*fx(m-1,n,k)            &
            +vp*dxh(m,n)*fy(m,n,k) - vm*dxh(m,n-1)*fy(m,n-1,k))/sqt(m,n)  &
            + flxz(k) - flxz(k-1)
        ff(m,n,k) = ff(m,n,k)*(hhq(m,n)/hmid) - rhs*tau/2.0/hzt(k)/hmid
       enddo

      endif
                
    enddo
   enddo
!$omp end parallel do  

     call syncborder_real(ff, nz+1)
     
     if(periodicity_x/=0) then
       call cyclize_x(ff,nx,ny,nz+1,mmm,mm)
     end if
     
     if(periodicity_y/=0) then
       call cyclize_y(ff,nx,ny,nz+1,nnn,nn)
     end if

!Crossbeam step

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
      fx0(m,n,:) = fx(m,n,:)
      fy0(m,n,:) = fy(m,n,:)
      fz0(m,n,:) = fz(m,n,:)
    enddo
   enddo
!$omp end parallel do  


!$omp parallel do private(m,n,k,ind1,ind2,u,v,w,flx1) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end
      
      if (lcu(m,n)>0.5) then
       do k=2,nz
        u = (uh(m  ,n,k-1)*dz(k-1)+uh(m  ,n,k)*dz(k))/2.0
        ind1= nint( 0.5 - sign(0.5,u) )
        ind2= nint(  - sign(1.0,u) )
        fx(m,n,k) = (1.0 + lcu(m+ind2,n))*ff(m+ind1,n,k) - lcu(m+ind2,n)*fx0(m+ind2,n,k)
       enddo
      endif

      if (lcv(m,n)>0.5) then
       do k=2,nz
        v = (vh(m,n  ,k-1)*dz(k-1)+vh(m,n  ,k)*dz(k))/2.0
        ind1= nint( 0.5 - sign(0.5,v) )
        ind2= nint(  - sign(1.0,v) )
        fy(m,n,k) = (1.0 + lcv(m,n+ind2))*ff(m,n+ind1,k) - lcv(m,n+ind2)*fy0(m,n+ind2,k)
       enddo
      endif

      if (lu(m,n)>0.5) then
       flx1(0) = 0.0
       do k=1,nz
        flx1(k) = fz0(m,n,k)
       enddo
       flx1(nz+1) = 0.0
       
       do k=1,nz
        w = (ww(m,n,k)+ww(m,n,k+1))/2.0
        ind1= nint( 0.5 - sign(0.5,w) )
        ind2= nint(  - sign(1.0,w) )
        fz(m,n,k) = (1.0 + cvert(k+ind2))*ff(m,n,k+ind1) - cvert(k+ind2)*flx1(k+ind2)
       enddo
      endif
                
    enddo
   enddo
!$omp end parallel do  

    call syncborder_real(fx, nz+1)
    call syncborder_real(fy, nz+1)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz+1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz+1,nnn,nn)            
    end if

!Corrector step

!$omp parallel do private(m,n,k,hmid,rhs,flxz) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end

      if (lu(m,n)>0.5) then
       hmid=(hhq(m,n)+hhq_n(m,n))/2.0
       
       do k=1,nz
       flxz(k) = (ww(m,n,k)+ww(m,n,k+1))/2.0*fz(m,n,k)
       enddo
       
       do k=2,nz
        up=(uh(m  ,n,k-1)*dz(k-1)+uh(m  ,n,k)*dz(k))/2.0
        um=(uh(m-1,n,k-1)*dz(k-1)+uh(m-1,n,k)*dz(k))/2.0
        vp=(vh(m,n  ,k-1)*dz(k-1)+vh(m,n  ,k)*dz(k))/2.0
        vm=(vh(m,n-1,k-1)*dz(k-1)+vh(m,n-1,k)*dz(k))/2.0

        rhs=(up*dyh(m,n)*fx(m,n,k) - um*dyh(m-1,n)*fx(m-1,n,k)            &
            +vp*dxh(m,n)*fy(m,n,k) - vm*dxh(m,n-1)*fy(m,n-1,k))/sqt(m,n)  &
            + flxz(k) - flxz(k-1)
        ff(m,n,k) = ff(m,n,k)*(hmid/hhq_n(m,n))  - rhs*tau/2.0/hzt(k)/hhq_n(m,n) 
       enddo

      endif
                
    enddo
   enddo
!$omp end parallel do  

     call syncborder_real(ff, nz+1)
     
     if(periodicity_x/=0) then
       call cyclize_x(ff,nx,ny,nz+1,mmm,mm)
     end if
     
     if(periodicity_y/=0) then
       call cyclize_y(ff,nx,ny,nz+1,nnn,nn)
     end if

endsubroutine turb_transport_cab_mod

!==========================================================================================================
subroutine turb_sdiff(ff,fx,fy,tau,mu,factor_mu)

real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), intent(in):: mu
real(4), intent(in):: tau, factor_mu
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), intent(inout):: ff
real(4), dimension(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1):: fx,fy

integer m,n,k
real(8) rhs

!$omp parallel do private(m,n) 
   do n=bnd_y1,bnd_y2
    do m=bnd_x1,bnd_x2   
     fx(m,n,:) = 0.0
     fy(m,n,:) = 0.0
    enddo
   enddo
!$omp end parallel do  

!$omp parallel do private(m,n,k) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    

     if(lcu(m,n)>0.5) then
      do k=2,nz      
       fx(m,n,k) = - (mu(m,n,k)+mu(m+1,n,k)+mu(m,n,k-1)+mu(m+1,n,k-1))/4.0*factor_mu     &
                *hhu_n(m  ,n)*dyh(m  ,n)/dxt(m  ,n)*(ff(m+1,n,k)-ff(m  ,n,k))
      enddo
     endif
     
     if(lcv(m,n)>0.5) then
      do k=2,nz 
       fy(m,n,k) = - (mu(m,n,k)+mu(m,n+1,k)+mu(m,n,k-1)+mu(m,n+1,k-1))/4.0*factor_mu     &
                *hhv_n(m,n  )*dxh(m,n  )/dyt(m,n  )*(ff(m,n+1,k)-ff(m,n  ,k))
      enddo     
     endif

    enddo
   enddo
!$omp end parallel do        

    call syncborder_real(fx,nz+1)
    call syncborder_real(fy,nz+1)

    if(periodicity_x/=0) then
      call cyclize_x(fx,nx,ny,nz+1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(fy,nx,ny,nz+1,nnn,nn)            
    end if

!$omp parallel do private(m,n,k,rhs) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    

     if(lu(m,n)>0.5) then
      do k=2,nz
       rhs= fx(m,n,k)-fx(m-1,n,k)+fy(m,n,k)-fy(m,n-1,k)  
       ff(m,n,k) = ff(m,n,k) - tau*rhs/hhq_n(m,n)/sqt(m,n)  
      enddo
     endif

    enddo
   enddo
!$omp end parallel do        

    call syncborder_real(ff, nz+1)

    if(periodicity_x/=0) then
      call cyclize_x(ff,nx,ny,nz+1,mmm,mm)
    end if

    if(periodicity_y/=0) then
      call cyclize_y(ff,nx,ny,nz+1,nnn,nn)            
    end if

endsubroutine turb_sdiff

endmodule mixing_routes