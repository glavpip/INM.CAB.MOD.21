!---------------------module for definition of basin grid arrays-----------------
module basin_grid
use mpi_parallel_tools
implicit none

include 'basinpar.fi'

!arrays of basin grids
real(4), allocatable:: lu(:,:),        &  !mask of t-grid
                       lu1(:,:),       &  !mask of t-grid (1 everywhere)
                       luu(:,:),       &  !mask of h-grid (0 on boundary)
                       luh(:,:),       &  !mask of h-grid (1 on boundary)
                       lcu(:,:),       &  !mask of u-grid (0 on boundary)
                       lcv(:,:),       &  !mask of v-grid (0 on boundary)
                       llu(:,:),       &  !mask of u-grid (1 on boundary)
                       llv(:,:)           !mask of v-grid (1 on boundary)

integer, allocatable:: lbasins(:,:)       !integer mask

real(4), allocatable:: hhh(:,:),      &  !ocean depth on luh (h-points)
                       hhh_n(:,:),    &  !ocean depth on luh (h-points) at previous step
                   hhq_rest(:,:),     &  !ocean depth on lu  (t-points) at rest state
                   hhu_rest(:,:),     &  !ocean depth on lcu (u-points) at rest state
                   hhv_rest(:,:),     &  !ocean depth on lcv (v-points) at rest state
                       hhq(:,:),      &  !ocean depth on lu  (t-points) 
                       hhq_n(:,:),    &  !ocean depth on lu  (t-points) at previous step
                       hhu(:,:),      &  !ocean depth on lcu (u-points)
                       hhu_n(:,:),    &  !ocean depth on lcu (u-points) at previous step
                       hhv(:,:),      &  !ocean depth on lcv (v-points)
                       hhv_n(:,:),    &  !ocean depth on lcv (v-points) at pre-previous step
                       rlh_s(:,:),    &  !main (sin) coriolis parameter on h-points
                       rlh_c(:,:),    &  !2-nd (cos) coriolis parameter on h-points
                    z(:), zw(:),      &  !vertical sigma-levels (t-points and w-points)
                  hzt(:), dz(:),      &  !steps between t-levels and w-levels
            dxt(:,:),dyt(:,:),        &  !horizontal grid steps between   t-points (in meters)
            dx (:,:),dy (:,:),        &  !horizontal grid steps between u,v-points (in meters)
            dxh(:,:),dyh(:,:),        &  !horizontal grid steps between   h-points (in meters)
            dxb(:,:),dyb(:,:)            !horizontal grid steps between v,u-points (in meters)

real(4), allocatable:: sqt(:,:),      &  !grid area in T-points
                       squ(:,:),      &  !grid area in U-points
                       sqv(:,:),      &  !grid area in V-points
                       sqh(:,:)          !grid area in H-points

 real(8), allocatable::     xt(:),yt(:),        &  !horizontal t-grid            x- and y-coordinates (in degrees)
                            xu(:),yv(:)            !horizontal u-grid and v-grid x- and y-coordinates (in degrees)


 real(8), allocatable::   geo_lon_t(:,:),   &    !geographical longitudes of T-points
                          geo_lat_t(:,:),   &    !geographical latitudes  of T-points
                          geo_lon_u(:,:),   &    !geographical longitudes of U-points
                          geo_lat_u(:,:),   &    !geographical latitudes  of U-points
                          geo_lon_v(:,:),   &    !geographical longitudes of V-points
                          geo_lat_v(:,:),   &    !geographical latitudes  of V-points
                          geo_lon_h(:,:),   &    !geographical longitudes of H-points
                          geo_lat_h(:,:),   &    !geographical latitudes  of H-points
                       rotvec_coeff(:,:,:)       !cos and sin of angles between coordinate lines
 endmodule basin_grid
!---------------------end module for definition of basin grid arrays-----------------
