module depth_routes
use constants
use basin_grid
use ocalg_routes
implicit none

contains

!============================================================
subroutine hh_init(hq, hu, hv, hh, sh, hq_r, hu_r, hv_r)
 
 real(8), dimension(bnd_x1:bnd_x2, bnd_y1:bnd_y2), intent(in):: sh
 real(4), dimension(bnd_x1:bnd_x2, bnd_y1:bnd_y2), intent(in):: hq_r
 real(4), dimension(bnd_x1:bnd_x2, bnd_y1:bnd_y2), intent(inout):: hq, hu, hv, hh, hu_r, hv_r
 
 real(4) slu
 integer m,n

!$omp parallel do private(m,n) 
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1
        if(lu(m,n)>0.5) then
          hq(m,n) =max(hq_r(m,n) + sh(m,n)*float(nonlinear_free_surface), hhqmin)
        endif
       end do
    end do
!$omp end parallel do

!$omp parallel do private(m,n,slu) 
      do n=ny_start-1,ny_end
       do m=nx_start-1,nx_end
        
        if(llu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
          slu=lu(m,n)+lu(m+1,n)
          hu_r(m,n)=(  hq_r(m  ,n)*lu(m  ,n)*sqt(m  ,n)   &
                    +  hq_r(m+1,n)*lu(m+1,n)*sqt(m+1,n) )/slu/squ(m,n)
            hu(m,n)=(  hq(m  ,n)*lu(m  ,n)*sqt(m  ,n)      &
                    +  hq(m+1,n)*lu(m+1,n)*sqt(m+1,n) )/slu/squ(m,n)
        endif

        if(llv(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
          slu=lu(m,n)+lu(m,n+1)
          hv_r(m,n)=(  hq_r(m,n  )*lu(m,n  )*sqt(m,n  )    &
                    +  hq_r(m,n+1)*lu(m,n+1)*sqt(m,n+1) )/slu/sqv(m,n)
          hv(m,n)=(  hq(m,n  )*lu(m,n  )*sqt(m,n  )        &
                  +  hq(m,n+1)*lu(m,n+1)*sqt(m,n+1) )/slu/sqv(m,n)
        endif

        if(luh(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
          slu=lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)
          hh(m,n)=(  hq(m  ,n  )*lu(m  ,n  )*sqt(m  ,n  )       &
                  +  hq(m+1,n  )*lu(m+1,n  )*sqt(m+1,n  )       &
                  +  hq(m  ,n+1)*lu(m  ,n+1)*sqt(m  ,n+1)       &
                  +  hq(m+1,n+1)*lu(m+1,n+1)*sqt(m+1,n+1) )/slu/sqh(m,n)
        endif

       end do
    end do
!$omp end parallel do

      call syncborder_real(hu_r, 1)
      call syncborder_real(hu,   1)
      call syncborder_real(hv_r, 1)
      call syncborder_real(hv,   1)
      call syncborder_real(hh,   1)
      
      if(periodicity_x/=0) then
        call cyclize_x(hu_r,nx,ny,1,mmm,mm)
        call cyclize_x(hu, nx,ny,1,mmm,mm)
        call cyclize_x(hv_r,nx,ny,1,mmm,mm)
        call cyclize_x(hv, nx,ny,1,mmm,mm)                          
        call cyclize_x(hh, nx,ny,1,mmm,mm)
      end if

      if(periodicity_y/=0) then
        call cyclize_y(hu_r,nx,ny,1,nnn,nn)
        call cyclize_y(hu, nx,ny,1,nnn,nn)
        call cyclize_y(hv_r,nx,ny,1,nnn,nn)
        call cyclize_y(hv, nx,ny,1,nnn,nn)                      
        call cyclize_y(hh, nx,ny,1,nnn,nn)
      end if

endsubroutine hh_init

!============================================================
subroutine hh_update(hq, hu, hv, hh, sh, hq_r)

 real(8), dimension(bnd_x1:bnd_x2, bnd_y1:bnd_y2), intent(in):: sh
 real(4), dimension(bnd_x1:bnd_x2, bnd_y1:bnd_y2), intent(in):: hq_r
 real(4), dimension(bnd_x1:bnd_x2, bnd_y1:bnd_y2), intent(inout):: hq, hu, hv, hh

 real(4) slu
 integer m,n

!$omp parallel do private(m,n) 
      do n=ny_start-1,ny_end+1
        do m=nx_start-1,nx_end+1
         if(lu(m,n)>0.5) then      
          hq(m,n) =max(hq_r(m,n) + sh(m,n)*float(nonlinear_free_surface), hhqmin)        
         endif
        end do
      end do
!$omp end parallel do

!$omp parallel do private(m,n,slu) 
      do n=ny_start-1,ny_end
       do m=nx_start-1,nx_end
        
        if(llu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
          slu=lu(m,n)+lu(m+1,n)
          hu(m,n)=( hq(m  ,n)*lu(m  ,n)*sqt(m  ,n)   &
                  + hq(m+1,n)*lu(m+1,n)*sqt(m+1,n) )/slu/squ(m,n)
        endif

        if(llv(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
          slu=lu(m,n)+lu(m,n+1)
          hv(m,n)=( hq(m,n  )*lu(m,n  )*sqt(m,n  )   &
                  + hq(m,n+1)*lu(m,n+1)*sqt(m,n+1) )/slu/sqv(m,n)
        endif

        if(luh(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
          slu=lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)
           hh(m,n)=( hq(m  ,n  )*lu(m  ,n  )*sqt(m  ,n)         &
                   + hq(m+1,n  )*lu(m+1,n  )*sqt(m+1,n)         &
                    +hq(m  ,n+1)*lu(m  ,n+1)*sqt(m  ,n+1)       &
                   + hq(m+1,n+1)*lu(m+1,n+1)*sqt(m+1,n+1) )/slu/sqh(m,n)
        endif

       end do
      end do
!$omp end parallel do

      call syncborder_real(hu, 1)
      call syncborder_real(hv, 1)
      call syncborder_real(hh, 1)

      if(periodicity_x/=0) then
        call cyclize_x(hu, nx,ny,1,mmm,mm)
        call cyclize_x(hv, nx,ny,1,mmm,mm)     
        call cyclize_x(hh, nx,ny,1,mmm,mm)
      end if

      if(periodicity_y/=0) then
        call cyclize_y(hu, nx,ny,1,nnn,nn)
        call cyclize_y(hv, nx,ny,1,nnn,nn)     
        call cyclize_y(hh, nx,ny,1,nnn,nn)
      end if

endsubroutine hh_update

!======================================================================
subroutine depth_ave(u,utr,mask,kbcl)
!-----------------------------------------------------------------------
! makes baroclinic and barotropic velocity component
! utr= <u>, u' = u - <u>.
 real(4) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
 real(4) mask(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(8) utr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),ff0
 integer m, n, k
 integer kbcl      ! key to remove barotropic component from the field: 0 - not, 1 - yes

!$omp parallel do private(m,n,k,ff0)
      !do n=bnd_y1+1,bnd_y2-1
      !do m=bnd_x1+1,bnd_x2-1
      do n = ny_start-1, ny_end+1
        do m = nx_start-1, nx_end+1
          if (mask(m,n)>0.5) then
            ff0 = 0.0d0
            do k=1,nz
              ff0 = ff0 + dble( u(m,n,k)*dz(k) )
            enddo
            utr(m,n) = ff0
            u(m,n,1:nz) = u(m,n,1:nz)-float(kbcl)*sngl(ff0)
          endif
        enddo
      enddo
!$omp end parallel do

endsubroutine depth_ave

endmodule depth_routes