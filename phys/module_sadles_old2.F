!WRF:MiDEL_LAYER:PHYSICS
MODULE module_sadles
!
! SADLES: Simple Actuator Disc for LES
! 
! History:
!    SADLES 2.0 (Dec.2025): - Fixed grid structure access (grid%ht)
!                - Implemented Analytical Line Clipping for horizontal distribution
!                - Implemented Exact CDF integration for vertical distribution
!                - Corrected AGL logic (subtracting terrain height)
!                - Renamed yaw_ang to nac_head (Nacelle Heading),nac_head is defined in windturbines.txt or windturbines-ij.txt
!                - Added: config_flags%sadles_yawing, if =0, then nac_head is not updated, which enables wake defection
!    SADLES 1.2: - Remove direct option for induction factor
!    SADLES 1.0: - Initial version by Hai Bui, hai.bui@uib.no
!
! Reference: 
!    Bui, H., Bakhoday-Paskyabi, M., & Mohammadpour-Penchah, M. (2024). Implementation of a Simple Actuator Disk for Large-Eddy Simulation in the Weather Research and Forecasting Model (WRF-SADLES v1. 2) for wind turbine wake simulation. Geoscientific Model Development, 17(10), 4447-4465.
!

  USE module_driver_constants, ONLY : max_domains
  USE module_model_constants, ONLY :  piconst
  USE module_domain
  USE module_llxy
  USE module_configure, ONLY : grid_config_rec_type

  IMPLICIT NONE

  INTEGER           :: max_turbines, sadles_domains
  INTEGER, PARAMETER :: sad_segment=20

  ! Turbine Location (Real indices)
  REAL, DIMENSION(:,:), ALLOCATABLE :: tb_i, tb_j  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: tb_types
  
  ! Ambient grid points 
  INTEGER, DIMENSION(:,:), ALLOCATABLE ::  tb_k, amb_spd  
  
  ! Dims: maxdoms 
  INTEGER, DIMENSION(max_domains) :: sadles_opt, nturbines

  ! Wind turbines specification
  INTEGER, PARAMETER :: max_tb_types = 50, max_tb_steps = 50    
  REAL, DIMENSION(max_tb_types) :: tb_hubheight, tb_diameter, tb_sthcoef, &
                                   tb_npower, tb_nradius, tb_area 
  INTEGER, DIMENSION(max_tb_types) :: tb_steps
  REAL, DIMENSION(max_tb_types,max_tb_steps) :: tb_ws, tb_thcoef, tb_power, &
                                                tb_Pcoef, tb_hs, tb_aif   

  ! Actuator disc discretization 
  INTEGER :: sad_max_points, row_max_points 

  ! Dimension: (max_turbines, max_domains, sad_max_points)
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: sad_i, sad_k, sad_j 
  REAL, DIMENSION(:,:,:), ALLOCATABLE    :: sad_Pow, sad_AreaFact, sad_Area           
  INTEGER, DIMENSION(:,:), ALLOCATABLE   :: sad_npoints           
  REAL,  DIMENSION(:,:), ALLOCATABLE     :: th_coef, Pw_coef, ai_fact 

  ! Temporary Row Buffers
  integer, dimension(:), ALLOCATABLE  ::  row_i, row_j
  real, dimension(:), ALLOCATABLE     ::  row_Area, row_AreaFact
  
  ! Orientation: nac_head = Absolute direction turbine is facing (Math angle: 0=E, 90=N)
  real, dimension(:), ALLOCATABLE     ::  nac_head 

  ! Output Buffers
  INTEGER, PARAMETER :: out_buffer = 100, foutput=20
  INTEGER :: i_buffer
  REAL, DIMENSION(:,:,:), ALLOCATABLE  :: out_Pow, out_MisAng, out_HubSpd
  REAL, DIMENSION(out_buffer) :: out_hours

#if defined(DM_PARALLEL) && !defined(STUBMPI)
      INCLUDE "mpif.h"
#endif

CONTAINS

! ====================================================================
! Main Subroutine: SADLES Tendency
! ====================================================================
SUBROUTINE  sadles_tendency(id, grid, dx, z3d, hgt, &  ! <--- Added hgt argument
               config_flags &
              ,u,v,mut                  &
              ,ru_tendf,rv_tendf        &
              ,ims,ime,jms,jme,kms,kme  &
              ,ips,ipe,jps,jpe          &
              )  
   IMPLICIT NONE

  INTEGER, INTENT(IN) :: id
  TYPE(domain) , TARGET          :: grid
  TYPE (grid_config_rec_type) :: config_flags
  REAL,INTENT(IN)  :: dx
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: u,v, z3d
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: mut
  ! ADDED: Terrain Height Input (Pass grid%ht from caller)
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: hgt 
  
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf, rv_tendf

  LOGICAL, EXTERNAL :: wrf_dm_on_monitor

  ! Local variables
  integer :: itb, ntb, ttype, ir,jr,kr, ipt, ii,kk,jj
  real ::  nrad, speedn, speedt, angle, cosab, Amb_speedn
  real ::  hours, ai_fact_safe, u_hub, v_hub, u_tmp, v_tmp, local_sum
  integer :: ierror, sadcom
  
  ! String declarations (Restored)
  CHARACTER(LEN=256) :: tb_num, input, output, fformat

  hours = (grid%itimestep-1)*grid%dt/3600.  

  ! Update output buffer index
#if defined(DM_PARALLEL)  && !defined(STUBMPI)
  CALL WRF_GET_DM_COMMUNICATOR(sadcom)
  IF ( wrf_dm_on_monitor() ) THEN
     if (i_buffer.le.out_buffer) then
        i_buffer=i_buffer+1
     else
        i_buffer=1
     endif
     out_hours(i_buffer) = hours
  END IF
  CALL wrf_dm_bcast_integer(i_buffer,1)
#else
  if (i_buffer.le.out_buffer) then
        i_buffer=i_buffer+1
  else
        i_buffer=1
  endif
  out_hours(i_buffer) = hours
#endif

  ! 1. Adjust Yaw
  if (config_flags%sadles_yawing .eq. 1) then
     call  Yaw_turbines(id,grid,dx,z3d, hgt, u,v    &
                    ,ims,ime,jms,jme,kms,kme  &
                    ,ips,ipe,jps,jpe )
  end if
 
  ! 2. Recalculate Geometric Intersection 
  call  calc_AreaFact(id,grid, dx,z3d, hgt        &
                    ,ims,ime,jms,jme,kms,kme  &
                    ,ips,ipe,jps,jpe )
 
  ! 3. Calculate Coefficients based on Hub Speed
  ntb = nturbines(id)
  th_coef(1:ntb,id) = 0.
  Pw_coef(1:ntb,id) = 0.
  ai_fact(1:ntb,id) = 0.
  out_HubSpd(1:ntb,id,i_buffer) = 0.
  out_MisAng(1:ntb,id,i_buffer) = 0.
  out_Pow(1:ntb,id,i_buffer) = 0.

   do itb=1,ntb
      ! Use NINT to access specific grid columns to match parallel ownership
      ir = nint(tb_i(itb,id))
      jr = nint(tb_j(itb,id))
      kr = tb_k(itb,id) 
      ttype = tb_types(itb,id)
      
      if ( ( ir .le. ipe-1 ) .and.  ( ir .ge. ips ) .and. (jr .le. jpe-1 ) .and. (jr.ge. jps ) ) then
         u_hub = 0.5*(u(ir,kr,jr) + u(ir+1,kr,jr)) 
         v_hub = 0.5*(v(ir,kr,jr) + v(ir,kr,jr+1))
         out_HubSpd(itb,id,i_buffer) = sqrt( u_hub**2 + v_hub**2) 
         
         angle = atan2(-v_hub, -u_hub) 
         cosab = cos(angle - nac_head(itb))
         speedn = out_HubSpd(itb,id,i_buffer) * cosab
         
         ai_fact(itb,id) = Axial_induction(id,itb,speedn)
         th_coef(itb,id) = Thrust_coef(id,itb,speedn)             
         Pw_coef(itb,id) = Power_coef(id,itb,speedn)  
         out_MisAng(itb,id,i_buffer) = (angle - nac_head(itb)) * 180./piconst
      end if 
   end do

#if defined(DM_PARALLEL)  && !defined(STUBMPI)
  CALL MPI_ALLREDUCE(out_HubSpd(1:ntb,id,i_buffer), out_HubSpd(1:ntb,id,i_buffer), ntb, MPI_REAL, MPI_SUM, sadcom, ierror)
  CALL MPI_ALLREDUCE(out_MisAng(1:ntb,id,i_buffer), out_MisAng(1:ntb,id,i_buffer), ntb, MPI_REAL, MPI_SUM, sadcom, ierror)
  CALL MPI_ALLREDUCE(th_coef(1:ntb,id), th_coef(1:ntb,id), ntb, MPI_REAL, MPI_SUM, sadcom, ierror)
  CALL MPI_ALLREDUCE(ai_fact(1:ntb,id), ai_fact(1:ntb,id), ntb, MPI_REAL, MPI_SUM, sadcom, ierror)
#endif

  ! 4. Apply Tendencies
  do itb=1,ntb
      ir = nint(tb_i(itb,id))
      jr = nint(tb_j(itb,id))
      ttype = tb_types(itb,id)
      nrad = tb_nradius(ttype)
      
      if ( ( ir-nrad .le. ipe ) .and.  ( ir+nrad .ge. ips ) .and. (jr-nrad .le. jpe ) .and. (jr+nrad.ge. jps ) ) then
         
         ! Safety clamp for induction factor
         ai_fact_safe = min(ai_fact(itb,id), 0.5)

         do ipt=1, sad_npoints(itb,id)
            ii = sad_i(itb,id,ipt)
            kk = sad_k(itb,id,ipt)
            jj = sad_j(itb,id,ipt)
            
            if ( ( ii .le. ipe ) .and.  ( ii .ge. ips ) .and. (jj .le. jpe ) .and. (jj.ge.jps ) ) then
               u_tmp = 0.5*(u(ii,kk,jj) + u(ii+1,kk,jj))
               v_tmp = 0.5*(v(ii,kk,jj) + v(ii,kk,jj+1))
               angle = atan2(-v_tmp, -u_tmp)     
               
               speedn = sqrt(u_tmp**2 + v_tmp**2) * cos(angle - nac_head(itb))  
               Amb_speedn = speedn / (1.0 - ai_fact_safe)

               ! Drag Force Application (Factor 0.25 = 0.5 Physics * 0.5 Grid Splitting)
               ru_tendf(ii,kk,jj)   = ru_tendf(ii,kk,jj)   + 0.25 * mut(ii,jj) * th_coef(itb,id) * (Amb_speedn**2) * cos(nac_head(itb)) * sad_AreaFact(itb,id,ipt) 
               ru_tendf(ii+1,kk,jj) = ru_tendf(ii+1,kk,jj) + 0.25 * mut(ii,jj) * th_coef(itb,id) * (Amb_speedn**2) * cos(nac_head(itb)) * sad_AreaFact(itb,id,ipt) 
               rv_tendf(ii,kk,jj)   = rv_tendf(ii,kk,jj)   + 0.25 * mut(ii,jj) * th_coef(itb,id) * (Amb_speedn**2) * sin(nac_head(itb)) * sad_AreaFact(itb,id,ipt)
               rv_tendf(ii,kk,jj+1) = rv_tendf(ii,kk,jj+1) + 0.25 * mut(ii,jj) * th_coef(itb,id) * (Amb_speedn**2) * sin(nac_head(itb)) * sad_AreaFact(itb,id,ipt) 

               sad_Pow(itb,id,ipt) = 0.5 * Pw_coef(itb,id) * sad_Area(itb,id,ipt) * Amb_speedn**3 
            end if
         end do
      end if
       
#if defined(DM_PARALLEL)  && !defined(STUBMPI)
      local_sum = SUM(sad_Pow(itb,id,:))  
      CALL MPI_ALLREDUCE( local_sum, out_Pow(itb,id,i_buffer), 1, MPI_REAL, MPI_SUM, sadcom, ierror)
#else
      out_Pow(itb,id,i_buffer) = SUM(sad_Pow(itb,id,:))
#endif
  end do

  ! 5. Write Output
#if defined(DM_PARALLEL)  && !defined(STUBMPI)
  IF ( wrf_dm_on_monitor() ) THEN
#endif
     if (i_buffer.eq.out_buffer) then
         write(tb_num,'(I3)'),nturbines(id)
         write(*,*)"SADLES: writing turbine output"
         fformat="(f12.7,"//trim(tb_num)//"(f8.2))"

         write(output,"('sadles_MisAng.d',i0.2)") id
         OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND',ACTION="WRITE")
         do i_buffer=1,out_buffer
             write(foutput,fformat) out_hours(i_buffer), out_MisAng(1:ntb,id,i_buffer)
         enddo  
         CLOSE(foutput)

         write(output,"('sadles_HubSpd.d',i0.2)") id
         OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND',ACTION="WRITE")
         do i_buffer=1,out_buffer
             write(foutput,fformat) out_hours(i_buffer), out_HubSpd(1:ntb,id,i_buffer)
         enddo
         CLOSE(foutput)

         write(output,"('sadles_Power.d',i0.2)") id
         OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND',ACTION="WRITE")
         do i_buffer=1,out_buffer
             write(foutput,fformat) out_hours(i_buffer), out_Pow(1:ntb,id,i_buffer)/1000. 
         enddo
         CLOSE(foutput)

         i_buffer = 0
     endif
#if defined(DM_PARALLEL)   && !defined(STUBMPI)   
  END IF
  CALL wrf_dm_bcast_integer(i_buffer,1)
#endif

END SUBROUTINE sadles_tendency

! ====================================================================
! Yaw Control 
! ====================================================================
  subroutine Yaw_turbines(id,grid, dx, z3d, hgt, u,v    & 
              ,ims,ime,jms,jme,kms,kme  &
              ,ips,ipe,jps,jpe          &
              )
  Implicit None
  TYPE(domain) , TARGET          :: grid
  INTEGER, INTENT(IN) :: id
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
  REAL, INTENT(IN) :: dx
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: z3d,u,v
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: hgt

  REAL, DIMENSION(kms:kme)  :: z
  ! Local Real coordinates for NINT usage
  real :: ir, jr 
  integer :: k,nz,ntb, itb, ttype, ii, jj
  real :: tb_z, uhub, vhub, hgt_at_turb
  integer :: khub
  
  ntb = nturbines(id) 
  nz = kme-kms
  
  do itb=1,ntb
      ir = tb_i(itb,id)
      jr = tb_j(itb,id)
      ttype=tb_types(itb,id)
      tb_z = tb_hubheight(ttype)
      
      ! Parallel Ownership: Use NINT logic
      if ( ( nint(ir) .le. ipe ) .and.  ( nint(ir) .gt. ips ) .and. &
           ( nint(jr) .le. jpe ) .and.  ( nint(jr) .gt. jps ) ) then
          
          ii = nint(ir)
          jj = nint(jr)
          hgt_at_turb = hgt(ii, jj)
          z(:) = z3d(ii,:,jj) - hgt_at_turb
          
          khub=-1
          do k=1, nz
              if ( (z(k).ge.tb_z) .and. (khub.eq.-1) ) then
                 khub = k
                 exit
              end if
          enddo
          uhub = 0.5*( u(ii,khub,jj) + u(ii+1,khub,jj))
          vhub = 0.5*( v(ii,khub,jj) + v(ii,khub,jj+1))
          nac_head(itb) = atan2(-vhub,-uhub)
       end if 
  end do    
  end subroutine Yaw_turbines


! ====================================================================
! Discretization Driver
! ====================================================================
  SUBROUTINE  calc_AreaFact(id,grid, dx, z3d, hgt & 
              ,ims,ime,jms,jme,kms,kme  &
              ,ips,ipe,jps,jpe          &
              )
  Implicit None
  TYPE(domain) , TARGET          :: grid
  INTEGER, INTENT(IN) :: id
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
  REAL, INTENT(IN) :: dx
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: z3d
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: hgt

  REAL, DIMENSION(kms:kme)  :: tb_Areaz, tb_Faz, tb_Rz, z
  integer :: k,nz, itb, ii,jj, ttype, kbot,khub,ktop
  real :: tb_z, tb_r, ir, jr, angle, &
          ti0,ti1,tj0,tj1, hgt_at_turb
  integer :: row_n_points
  integer :: icell,ierror, ntb
  integer :: sadcom
  
  nz = kme-kms
  ntb = nturbines(id)
  
  sad_i(1:ntb,id,:)=0
  sad_k(1:ntb,id,:)=0
  sad_j(1:ntb,id,:)=0
  sad_AreaFact(1:ntb,id,:) = 0.
  sad_Area(1:ntb,id,:) = 0.
  sad_Pow(1:ntb,id,:)=0.
  sad_npoints(1:ntb,id) =0
  tb_k(1:ntb,id)  = 0

  do itb=1,ntb
      ttype=tb_types(itb,id)
      ir = tb_i(itb,id)
      jr = tb_j(itb,id)
      tb_z = tb_hubheight(ttype)
      tb_r = tb_diameter(ttype)/2.0
      
      ! Use NINT for ownership logic
      if ( ( nint(ir) .le. ipe ) .and.  ( nint(ir) .gt. ips ) .and. &
           ( nint(jr) .le. jpe ) .and.  ( nint(jr) .gt. jps ) ) then
          
          hgt_at_turb = hgt(nint(ir), nint(jr))
          z(:) = z3d(nint(ir),:,nint(jr)) - hgt_at_turb
          
          CALL vert_AreaFact(tb_z, tb_r, z, nz, dx, tb_Areaz, tb_Faz, tb_Rz, kbot, khub, ktop)
          
          angle = nac_head(itb)
          tb_k(itb,id) = khub

          do k=kbot,ktop
             ! Line coordinates based on Geometric Radius
             ti0 = ir - tb_Rz(k)*sin(-angle)
             ti1 = ir + tb_Rz(k)*sin(-angle)
             tj0 = jr - tb_Rz(k)*cos(angle)
             tj1 = jr + tb_Rz(k)*cos(angle)
             
             call horz_AreaSplit(tb_Faz(k),ti0,tj0,ti1,tj1,row_max_points,row_n_points,row_i,row_j,row_AreaFact)
             call horz_AreaSplit(tb_Areaz(k),ti0,tj0,ti1,tj1,row_max_points,row_n_points,row_i,row_j,row_Area)
             
             do icell=1,row_n_points
                 sad_npoints(itb,id) = sad_npoints(itb,id) +1
                 if (sad_npoints(itb,id).gt.sad_max_points) then
                    write(*,*) "FATAL: sad_max_points too small."
                    stop
                 end if
                 sad_i(itb,id,sad_npoints(itb,id)) = row_i(icell)
                 sad_j(itb,id,sad_npoints(itb,id)) = row_j(icell)
                 sad_k(itb,id,sad_npoints(itb,id)) = k
                 sad_Area(itb,id,sad_npoints(itb,id)) = row_Area(icell)
                 sad_AreaFact(itb,id,sad_npoints(itb,id)) = row_AreaFact(icell)
             end do 
          end do
      endif
  end do

#if defined(DM_PARALLEL)  && !defined(STUBMPI)
  CALL WRF_GET_DM_COMMUNICATOR(sadcom)
  CALL MPI_ALLREDUCE(sad_npoints(1:ntb,id), sad_npoints(1:ntb,id),  ntb, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_i(1:ntb,id,:), sad_i(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_k(1:ntb,id,:), sad_k(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_j(1:ntb,id,:), sad_j(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_Pow(1:ntb,id,:), sad_Pow(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_AreaFact(1:ntb,id,:), sad_AreaFact(1:ntb,id,:), ntb*sad_max_points, MPI_REAL, MPI_SUM,  sadcom,ierror)
#endif  
  
  END SUBROUTINE calc_AreaFact

! ====================================================================
! Helper Functions
! ====================================================================
 real function Thrust_coef(id,tbid,speed)
     implicit none
     integer, intent(in) :: id,tbid
     real, intent(in) :: speed
     REAL, DIMENSION(max_tb_steps) :: ref_speed, ref_coef
     integer :: steps,ttype
     ttype = tb_types(tbid,id) 
     ref_coef = tb_thcoef(ttype,:)
     ref_speed = tb_hs(ttype,:)   
     steps=tb_steps(ttype)
     Thrust_coef = lin_interp(ref_coef,ref_speed,steps,speed)     
  end function Thrust_coef
 
 real function Power_coef(id,tbid,speed)
     implicit none
     integer, intent(in) :: id,tbid
     real, intent(in) :: speed
     REAL, DIMENSION(max_tb_steps) :: ref_speed, ref_coef
     integer :: steps,ttype
     ttype = tb_types(tbid,id) 
     ref_coef = tb_Pcoef(ttype,:)
     ref_speed = tb_hs(ttype,:)   
     steps=tb_steps(ttype)
     Power_coef = lin_interp(ref_coef,ref_speed,steps,speed)     
  end function Power_coef

 real function Axial_induction(id,tbid,speed)
     implicit none
     integer, intent(in) :: id,tbid
     real, intent(in) :: speed
     REAL, DIMENSION(max_tb_steps) :: ref_speed, ref_coef
     integer :: steps,ttype
     ttype = tb_types(tbid,id)
     ref_coef = tb_aif(ttype,:)
     ref_speed = tb_hs(ttype,:)  
     steps=tb_steps(ttype)
     Axial_induction = lin_interp(ref_coef,ref_speed,steps,speed)     
  end function Axial_induction

  ! ------------------------------------------------------------------
  ! Calculate Cumulative Circle Area from Bottom (-R) to z
  ! ------------------------------------------------------------------
  real function circle_segment(r,z)
     implicit none
     real, intent(in) :: r,z
     real :: term1, term2, offset, z_clamped, my_pi
     my_pi = 3.14159265

     if (z .le. -r) then
        circle_segment = 0.0
        return
     end if
     if (z .ge. r) then
        circle_segment = my_pi * r**2
        return
     end if

     z_clamped = max(-r, min(r, z))
     term1 = z_clamped * sqrt(r**2 - z_clamped**2)
     term2 = r**2 * asin(z_clamped / r)
     offset = my_pi * r**2 / 2.0
     circle_segment = term1 + term2 + offset
  end function circle_segment


  ! ------------------------------------------------------------------
  ! Calculate Vertical Area Distribution
  ! ------------------------------------------------------------------
  subroutine vert_AreaFact(tb_z, tb_r, z, nz, dx, Areas, Fas, Rs, kbot, khub, ktop)
     implicit none
     real, intent(in) :: tb_z, tb_r, dx
     integer, intent(in) :: nz
     real, dimension(nz), intent(in) :: z
     real, dimension(nz-1), intent(out) :: Areas, Fas, Rs
     integer, intent(out) :: kbot, khub, ktop

     real :: t_zbot, t_ztop, dx2
     real :: z_layer_bot, z_layer_top, area_bot, area_top
     real :: z_ov_min, z_ov_max, z_mid_ov, dist_mid
     integer :: k

     Fas = 0.0; Rs = 0.0; Areas = 0.0
     dx2 = dx * dx
     t_zbot = tb_z - tb_r
     t_ztop = tb_z + tb_r

     kbot = -1; khub = -1; ktop = -1
     do k = 1, nz
        if (kbot .eq. -1 .and. z(k) .ge. t_zbot) kbot = max(1, k-1)
        if (khub .eq. -1 .and. z(k) .ge. tb_z) khub = max(1, k-1)
        if (ktop .eq. -1 .and. z(k) .ge. t_ztop) then
           ktop = k-1; exit
        end if
     end do
     if (ktop .eq. -1) ktop = nz - 1
     if (kbot .eq. -1) kbot = 1

     do k = kbot, ktop
         if (k .lt. 1 .or. k .ge. nz) cycle
         z_layer_bot = z(k); z_layer_top = z(k+1)

         area_bot = circle_segment(tb_r, z_layer_bot - tb_z)
         area_top = circle_segment(tb_r, z_layer_top - tb_z)
         Areas(k) = max(0.0, area_top - area_bot)

         if ((z_layer_top - z_layer_bot) .gt. 0.0) then
             Fas(k) = Areas(k) / (z_layer_top - z_layer_bot) / dx2
         end if

         z_ov_min = max(z_layer_bot, t_zbot)
         z_ov_max = min(z_layer_top, t_ztop)
         if (z_ov_max .gt. z_ov_min) then
             z_mid_ov = (z_ov_min + z_ov_max) * 0.5
             dist_mid = abs(z_mid_ov - tb_z)
             if (dist_mid .lt. tb_r) then
                 Rs(k) = sqrt(tb_r**2 - dist_mid**2)
             else; Rs(k) = 0.0; end if
         else; Rs(k) = 0.0; end if
     end do
     Rs = Rs / dx 
  end subroutine vert_AreaFact

  ! ------------------------------------------------------------------
  ! Horizontal Discretization
  ! ------------------------------------------------------------------
  subroutine horz_AreaSplit(Area,ir0,jr0,ir1,jr1,max_points,n_points,grid_i,grid_j,grid_Area)
     implicit none
     real,intent(in) :: Area,ir0,jr0,ir1,jr1
     integer,intent(in) :: max_points
     integer,intent(out) :: n_points
     integer,dimension (max_points), intent(out) :: grid_i,grid_j
     real,dimension (max_points), intent(out)  :: grid_Area

     real :: di, dj, t, t_mid, i_mid, j_mid, temp_t
     real, dimension(200) :: t_vals 
     integer :: nt, k, i, swapped, istart, iend, jstart, jend
     real, parameter :: eps = 1.0e-6

     n_points = 0; grid_Area = 0.
     if (Area .le. 1.0e-9) return

     di = ir1 - ir0; dj = jr1 - jr0
     
     nt = 2; t_vals(1) = 0.0; t_vals(2) = 1.0

     if (abs(di) .gt. eps) then
         istart = floor(min(ir0, ir1)) + 1
         iend   = floor(max(ir0, ir1))
         do k = istart, iend
             t = (real(k) - ir0) / di
             if (t .gt. eps .and. t .lt. (1.0 - eps)) then
                 nt = nt + 1; if (nt .le. 200) t_vals(nt) = t
             end if
         end do
     end if

     if (abs(dj) .gt. eps) then
         jstart = floor(min(jr0, jr1)) + 1
         jend   = floor(max(jr0, jr1))
         do k = jstart, jend
             t = (real(k) - jr0) / dj
             if (t .gt. eps .and. t .lt. (1.0 - eps)) then
                 nt = nt + 1; if (nt .le. 200) t_vals(nt) = t
             end if
         end do
     end if

     do i = 1, nt-1
         swapped = 0
         do k = 1, nt-i
             if (t_vals(k) .gt. t_vals(k+1)) then
                 temp_t = t_vals(k); t_vals(k) = t_vals(k+1); t_vals(k+1) = temp_t
                 swapped = 1
             end if
         end do
         if (swapped .eq. 0) exit
     end do

     do k = 1, nt-1
         t = t_vals(k+1) - t_vals(k)
         if (t .gt. 1.0e-7) then
             t_mid = (t_vals(k) + t_vals(k+1)) * 0.5
             i_mid = ir0 + t_mid * di
             j_mid = jr0 + t_mid * dj
             n_points = n_points + 1
             if (n_points .gt. max_points) then
                 n_points = max_points; exit
             end if
             grid_i(n_points) = floor(i_mid)
             grid_j(n_points) = floor(j_mid)
             grid_Area(n_points) = Area * t
         end if
     end do
  end subroutine horz_AreaSplit

  ! ... (Rest of Init and Utils) ...
  SUBROUTINE  init_module_sadles(id,config_flags,xlong,xlat,sadles_initialized, &
                                  ims,ime,jms,jme,kms,kme,its,ite,jts,jte,ids,ide,jds,jde)
   IMPLICIT NONE
   integer id,ims,ime,jms,jme,kms,kme,ids,ide,jds,jde
   integer its,ite,jts,jte
   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
   TYPE (grid_config_rec_type) :: config_flags
   TYPE (PROJ_INFO) :: ts_proj
   logical :: sadles_initialized
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   integer :: max_dom
   real :: max_radius,min_dx,min_dz
   integer :: tmp_n
   real :: tmp_lat,tmp_lon, tmp_i, tmp_j,dx
   INTEGER i,j,k,i_tb, itype, ntb
   CHARACTER(LEN=256) :: tb_num,input,output,fformat
   
   max_dom = config_flags%max_dom
   dx = config_flags%dx
   max_radius = config_flags%sadles_maxradius
   min_dx = config_flags%sadles_mindx
   min_dz = config_flags%sadles_mindz
   sadles_domains = max_dom - id + 1
   sad_max_points = ceiling( (3.14159*max_radius**2)/(min_dx*min_dz) ) * 2
   row_max_points = ceiling(max_radius/min_dx)*4
   
#if defined(DM_PARALLEL)  && !defined(STUBMPI)
   IF ( wrf_dm_on_monitor() ) THEN
#endif
       if ( config_flags%windfarm_ij .eq. 1 ) then
         open(70,file='windturbines-ij.txt',form='formatted',status='old')
       else
         Call setup_wrf_mapproj(config_flags,ts_proj,xlong,xlat, &
                                          ims,ime,jms,jme,its,ite,jts,jte,ids,ide,jds,jde)
         open(70,file='windturbines.txt',form='formatted',status='old')
       end if

       max_turbines=0
 10    read(70,*,end=100) 
         max_turbines=max_turbines+1
         goto 10
 100   continue
       rewind (70)
#if defined(DM_PARALLEL)   && !defined(STUBMPI)      
    END IF
    CALL wrf_dm_bcast_integer(max_turbines,1)
#endif
    
    if (.not. sadles_initialized) then
        allocate (tb_types(max_turbines,id:max_dom), &
                  tb_i(max_turbines,id:max_dom), &
                  tb_j(max_turbines,id:max_dom), &
                  tb_k(max_turbines,id:max_dom),  &
                  amb_spd(max_turbines,id:max_dom),  &
                  sad_npoints(max_turbines,id:max_dom), &
                  th_coef(max_turbines,id:max_dom), &
                  Pw_coef(max_turbines,id:max_dom), &
                  ai_fact(max_turbines,id:max_dom), &
                  nac_head(max_turbines)) 

        allocate ( sad_i(max_turbines,id:max_dom,sad_max_points), & 
                   sad_k(max_turbines,id:max_dom,sad_max_points), &
                   sad_j(max_turbines,id:max_dom,sad_max_points), &
                   sad_Pow(max_turbines,id:max_dom,sad_max_points), &
                   sad_Area(max_turbines,id:max_dom,sad_max_points), &
                   sad_AreaFact(max_turbines,id:max_dom,sad_max_points))

        allocate ( out_Pow(max_turbines,id:max_dom,out_buffer), &
                   out_MisAng(max_turbines,id:max_dom,out_buffer), &
                   out_HubSpd(max_turbines,id:max_dom,out_buffer) )
        allocate ( row_i(row_max_points), row_j(row_max_points), &
                   row_Area(row_max_points),row_AreaFact(row_max_points))
        
        tb_i=-9999; tb_j=-9999; sadles_initialized=.true.
    endif
    
#if defined(DM_PARALLEL)    && !defined(STUBMPI)  
   IF ( wrf_dm_on_monitor() ) THEN
#endif
       i_tb=1; nturbines(id) = 0
       do i=1,max_turbines
          if ( config_flags%windfarm_ij .eq. 1 ) then
              read(70,*) tb_i(i,id), tb_j(i,id), tb_types(i,id), nac_head(i)
              nac_head(i)= nac_head(i)*piconst/180.
              nturbines(id)=nturbines(id)+1
          else
              read(70,*)  tmp_lat,tmp_lon,tmp_n
              CALL latlon_to_ij(ts_proj, tmp_lat,tmp_lon,  tmp_i, tmp_j)
              if ((tmp_i .le. ide) .and. (tmp_i.ge.ids) .and. (tmp_j .le. jde) .and. (tmp_j.ge.jds)) then
                  tb_i(i_tb,id) = tmp_i; tb_j(i_tb,id) = tmp_j 
                  tb_types(i_tb,id) = tmp_n; i_tb = i_tb+1 
                  nturbines(id)=nturbines(id)+1
              end if 
          endif
       end do
       close(70)

       tb_steps=0
       do i_tb=1,nturbines(id)
          itype = tb_types(i_tb,id)
          if (tb_steps(itype) .eq. 0 ) then 
             write(tb_num,*) itype; tb_num=trim(adjustl(tb_num))
             input="wind-turbine-"//trim(tb_num)//".tbl"
             OPEN(file=trim(input),unit=19,FORM='FORMATTED',STATUS='OLD')
             READ (19,*,ERR=132)tb_steps(itype)
             read(19,*,ERR=132)tb_hubheight(itype),tb_diameter(itype),tb_sthcoef(itype),tb_npower(itype)
             tb_nradius(itype) = tb_diameter(itype)/2/dx
             tb_area(itype) = piconst*tb_diameter(itype)**2/4. 
             do k=1,tb_steps(itype)
                 read(19,*,ERR=132) tb_ws(itype,k), tb_thcoef(itype,k), tb_power(itype,k) 
                 tb_power(itype,k) = tb_power(itype,k)*1000. 
                 tb_hs(itype,k) = 0.5*tb_ws(itype,k)*(1 + SQRT( 1 -  tb_thcoef(itype,k)  ) ) 
                 tb_aif(itype,k) = 0.5*(1 - SQRT( 1 - tb_thcoef(itype,k) ) )             
                 tb_Pcoef(itype,k) = tb_power(itype,k)/ (0.5*tb_area(itype)*tb_ws(itype,k)**3) 
             enddo
             close (19)
           end if
       end do
 132   continue
       
       write(output,"('sadles_info.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE') 
       write(foutput,"('N_turbines:',i3)") nturbines(id)
       write(tb_num,'(I3)'),nturbines(id)
       fformat="('turbine_type:',"//trim(tb_num)//"(I6))"
       write(foutput,fformat) tb_types(1:nturbines(id),id)
       fformat="('turbine_i:   ',"//trim(tb_num)//"(f6.2))"
       write(foutput,fformat) tb_i(1:nturbines(id),id)
       fformat="('turbine_j:   ',"//trim(tb_num)//"(f6.2))"
       write(foutput,fformat) tb_j(1:nturbines(id),id)
       close(foutput)

       write(output,"('sadles_MisAng.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE'); CLOSE(foutput)
       write(output,"('sadles_HubSpd.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE'); CLOSE(foutput)
       write(output,"('sadles_Power.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE'); CLOSE(foutput)
#if defined(DM_PARALLEL)     && !defined(STUBMPI)    
   END IF
   CALL wrf_dm_bcast_integer(nturbines,max_domains)
   CALL wrf_dm_bcast_integer(tb_types,max_turbines*sadles_domains)
   CALL wrf_dm_bcast_real(tb_i,max_turbines*sadles_domains)
   CALL wrf_dm_bcast_real(tb_j,max_turbines*sadles_domains)
   CALL wrf_dm_bcast_real(nac_head,max_turbines)
   CALL wrf_dm_bcast_real(tb_hubheight,max_tb_types)
   CALL wrf_dm_bcast_real(tb_diameter,max_tb_types)
   CALL wrf_dm_bcast_real(tb_area,max_tb_types)
   CALL wrf_dm_bcast_real(tb_nradius,max_tb_types)
   CALL wrf_dm_bcast_real(tb_npower,max_tb_types)
   CALL wrf_dm_bcast_real(tb_sthcoef,max_tb_types)
   CALL wrf_dm_bcast_real(tb_npower,max_tb_types)
   CALL wrf_dm_bcast_integer(tb_steps,max_tb_types)
   CALL wrf_dm_bcast_real(tb_ws,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_hs,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_aif,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_thcoef,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_power,max_tb_types*max_tb_steps)
#endif   
   i_buffer = 0
   sadles_opt(id)=config_flags%sadles_opt
  END SUBROUTINE init_module_sadles

  real function lin_interp(Y,X,Nx,x0)
   implicit none
   integer, intent(In) :: Nx
   real, dimension(Nx), intent(In) :: X, Y
   real, intent(In) :: x0
   integer :: i1,i2
   if (x0 .le.X(1)) then
      lin_interp = 0 
   else if (x0 .ge.X(Nx)) then
      lin_interp = 0 
   else
      do i2=2,Nx
         if (X(i2).ge.x0) then
            i1 = i2-1; exit
         end if
      end do
      lin_interp = Y(i1) + (x0 - X(i1))*( Y(i2)-Y(i1) ) / ( X(i2)-X(i1) )
   end if
  end function

  SUBROUTINE setup_wrf_mapproj(config_flags,ts_proj,xlong,xlat, &
                                            ims,ime,jms,jme,its,ite,jts,jte,ids,ide,jds,jde)
       integer ims,ime,jms,jme,ids,ide,jds,jde,its,ite,jts,jte
       TYPE (grid_config_rec_type) :: config_flags
       TYPE (PROJ_INFO) :: ts_proj
       REAL, DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
       REAL :: known_lat, known_lon
       CALL map_init(ts_proj)
       known_lat = xlat(its,jts)
       known_lon = xlong(its,jts)
       IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj, truelat1=config_flags%truelat1, &
                      lat1=known_lat, lon1=known_lon, knowni=REAL(its), knownj=REAL(jts), dx=config_flags%dx)
       ELSE IF (config_flags%map_proj == PROJ_LC) THEN
         CALL map_set(PROJ_LC, ts_proj, truelat1=config_flags%truelat1, truelat2=config_flags%truelat2, &
                      stdlon=config_flags%stand_lon, lat1=known_lat, lon1=known_lon, knowni=REAL(its), knownj=REAL(jts), dx=config_flags%dx)
       ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj, truelat1=config_flags%truelat1, stdlon=config_flags%stand_lon, &
                      lat1=known_lat, lon1=known_lon, knowni=REAL(its), knownj=REAL(jts), dx=config_flags%dx)
       END IF
  END SUBROUTINE setup_wrf_mapproj

SUBROUTINE geos_nudge ( grid , config_flags ,                        & 
                          ids, ide, jds, jde, kds, kde,              & 
                          ims, ime, jms, jme, kms, kme,              & 
                          ips, ipe, jps, jpe,                        & 
                          dt ,                                       & 
                          u_field , v_field,                         & 
                          u_base, v_base                          )    
   IMPLICIT NONE
   TYPE(domain), INTENT(IN)                           :: grid
   TYPE(grid_config_rec_type), INTENT(IN)             :: config_flags
   INTEGER, INTENT(IN)                                :: ids, ide, jds, jde, kds, kde
   INTEGER, INTENT(IN)                                :: ims, ime, jms, jme, kms, kme
   INTEGER, INTENT(IN)                                :: ips, ipe, jps, jpe
   REAL, INTENT(IN)                                   :: dt
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN) :: u_field, v_field
   REAL, DIMENSION(kms:kme), INTENT(INOUT)            :: u_base, v_base
   INTEGER :: i, j, k_ref, local_point_count, global_point_count
   REAL    :: local_sum_u, local_sum_v, global_sum_u, global_sum_v
   REAL    :: average_u_raw, average_v_raw, f, alpha, delta_Ug, delta_Vg, tendency_u, tendency_v, error_u, error_v
   LOGICAL, SAVE :: first_call = .TRUE.
   REAL, SAVE    :: u_hub_smooth = 0.0, v_hub_smooth = 0.0
   REAL, SAVE    :: u_hub_smooth_past = 0.0, v_hub_smooth_past = 0.0
   REAL, SAVE    :: time_since_last_update = 0.0
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
#if defined(DM_PARALLEL)  && !defined(STUBMPI)
   INTEGER :: ierror, sadcom
#endif
   k_ref = config_flags%nudge_lev
   f = config_flags%ideal_f
#if defined(DM_PARALLEL)  && !defined(STUBMPI)
      CALL WRF_GET_DM_COMMUNICATOR(sadcom)
      local_sum_u = SUM( u_field(ips:ipe, k_ref, jps:jpe) )
      local_sum_v = SUM( v_field(ips:ipe, k_ref, jps:jpe) )
      local_point_count = (ipe - ips + 1) * (jpe - jps + 1)
      CALL MPI_Allreduce( local_sum_u, global_sum_u, 1, MPI_REAL, MPI_SUM, sadcom, ierror )
      CALL MPI_Allreduce( local_sum_v, global_sum_v, 1, MPI_REAL, MPI_SUM, sadcom, ierror )
      CALL MPI_Allreduce( local_point_count, global_point_count, 1, MPI_INTEGER, MPI_SUM, sadcom, ierror )
#else
      global_sum_u = SUM( u_field(ids:ide, k_ref, jds:jde) )
      global_sum_v = SUM( v_field(ids:ide, k_ref, jds:jde) )
      global_point_count = (ide - ids + 1) * (jde - jds + 1)
#endif
   IF (global_point_count > 0) THEN
      average_u_raw = global_sum_u / REAL(global_point_count)
      average_v_raw = global_sum_v / REAL(global_point_count)
   ELSE
      average_u_raw = 0.0; average_v_raw = 0.0
   END IF
   IF (first_call) THEN
      u_hub_smooth = average_u_raw; v_hub_smooth = average_v_raw
      u_hub_smooth_past = average_u_raw; v_hub_smooth_past = average_v_raw
      time_since_last_update = 0.0; first_call = .FALSE.
   ELSE
      IF (config_flags%nudge_tsmooth > 0.0) THEN
          alpha = dt / config_flags%nudge_tsmooth
          u_hub_smooth = alpha * average_u_raw + (1.0 - alpha) * u_hub_smooth
          v_hub_smooth = alpha * average_v_raw + (1.0 - alpha) * v_hub_smooth
      ELSE
          u_hub_smooth = average_u_raw; v_hub_smooth = average_v_raw
      END IF
   END IF
   time_since_last_update = time_since_last_update + dt
   IF (time_since_last_update >= config_flags%nudge_tfilter) THEN
      tendency_u = (u_hub_smooth - u_hub_smooth_past) / time_since_last_update
      tendency_v = (v_hub_smooth - v_hub_smooth_past) / time_since_last_update
      error_u = (u_hub_smooth - config_flags%nudge_u) / config_flags%nudge_tau
      error_v = (v_hub_smooth - config_flags%nudge_v) / config_flags%nudge_tau
      delta_Ug = - (1.0 / f) * (tendency_v + error_v)
      delta_Vg =   (1.0 / f) * (tendency_u + error_u)
      u_base = u_base + delta_Ug; v_base = v_base + delta_Vg
      time_since_last_update = 0.0
      u_hub_smooth_past = u_hub_smooth; v_hub_smooth_past = v_hub_smooth
   END IF
END SUBROUTINE geos_nudge

END MODULE module_sadles
