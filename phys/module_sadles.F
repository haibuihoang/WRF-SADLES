!WRF:MiDEL_LAYER:PHYSICS


MODULE module_sadles
!
! An Simple Actuator Disc for LES (SADLES) application where one turbine maybe
! represented by one or several gridpoints
! History:
!    First version: WRFV4.3: use a simple actuator disc model with thrust coeffient 
!                   By Hai Bui, hai.bui@uib.no
!
!
! Reference: to be submitted to GMD
!
!

  USE module_driver_constants, ONLY : max_domains
  USE module_model_constants, ONLY :  piconst
  USE module_domain
  USE module_llxy
  USE module_dm, ONLY : wrf_dm_min_real
  USE module_configure, ONLY : grid_config_rec_type

  IMPLICIT NONE


  INTEGER           :: max_turbines, sadles_domains
  INTEGER, PARAMETER :: sad_segment=20
  !dims:max_turbines,max_doms
  REAL, DIMENSION(:,:), ALLOCATABLE :: tb_i, tb_j  !Note: REAL indices, e.g. 40.43, 50.12 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: tb_types  !

  INTEGER, DIMENSION(:,:), ALLOCATABLE ::  amb_i, amb_j, amb_k, amb_spd  !the ambient grid points (1 radius in front of the turbine)
  !Dims: maxdoms 
  INTEGER, DIMENSION(max_domains) :: sadles_opt, nturbines

  ! Wind turbines specification
  INTEGER, PARAMETER :: max_tb_types = 50, max_tb_steps = 50    
  REAL, DIMENSION(max_tb_types) :: tb_hubheight,tb_diameter,tb_sthcoef,tb_npower, tb_nradius, tb_area  ! tb_nradius: radius/dx
  INTEGER, DIMENSION(max_tb_types) :: tb_steps
  REAL, DIMENSION(max_tb_types,max_tb_steps) :: tb_ws,  tb_thcoef, tb_power, tb_ws1, tb_aif   ! tb_ws1 is the estimated windspeed at the AD, tb_aif: axial induction factor
    
  ! Actuator disc discretization 

  ! Max point per turbines: suppose we have a 60m radius, and 20m grid in both
  ! horionzontal and vertical, sad_max_points is about pi*R^2/(dx**2) ~ 29
  ! points, we choose 50 points for safety
  INTEGER :: sad_max_points,row_max_points  ! Max point per turbine, max_point per row
  ! Dimension: max_turbines,id:max_domains,sad_max_points
  ! These will be calculated in: sadles_yaw subroutine
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: sad_i, sad_k, sad_j !Dimension(n_domains,n_turbines,sad_max_points) 
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: sad_AreaFact           !Dimension(n_domains,n_turbines,sad_max_points)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: sad_npoints !Dimension(n_domains,n_turbines)
  REAL,  DIMENSION(:,:), ALLOCATABLE   :: th_coef,ai_fact, tke_coef

  integer, dimension(:), ALLOCATABLE  ::  row_i,row_j
  real, dimension(:), ALLOCATABLE ::  row_Area

  !These the buffer for writing the turbine output, should not write out every timestep!
  INTEGER, PARAMETER :: out_buffer = 100, foutput=20
  INTEGER :: i_buffer
  REAL, DIMENSION(:,:,:), ALLOCATABLE  :: out_Pow, out_AmbSpd, out_HubSpd   !Dims: (n_domains,n_turbines,out_buffer)
  REAL, DIMENSION(out_buffer) :: out_hours
                             
#if defined(DM_PARALLEL) && !defined(STUBMPI)
      INCLUDE "mpif.h"
#endif

  ! Note:
  !
  ! the dry air collum mass are grid%muu and grid%muv at u and v grid points
  !

CONTAINS


!Must call sadles_yaw first!
  SUBROUTINE  sadles_tendency(id, grid, ctke_fact &
              ,u,v,mut                  &
              ,ru_tendf,rv_tendf        &
              ,tke_tend                 &   ! Is this coupled or not???
              ,ims,ime,jms,jme,kms,kme  &
              ,ips,ipe,jps,jpe          &
         
              )  
!
!
!
  INTEGER, INTENT(IN) :: id
  TYPE(domain) , TARGET          :: grid
  REAL, INTENT(IN) :: ctke_fact
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: u,v
  REAL, DIMENSION(ims:ime,jms:jme), INTENT(IN) :: mut
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf,rv_tendf,tke_tend

  LOGICAL, EXTERNAL :: wrf_dm_on_monitor

! local variables
  integer :: itb, ntb, ttype, ir,jr,kr, sadcom,ierror, ipt, ii,kk,jj
  real ::  nrad,speed
  real :: C_pow, C_tke
  real, parameter :: rho_std=1.225
  CHARACTER*256 tb_num,input,output,fformat
  real :: hours, ai_fact_temp
  integer :: opt

  opt=sadles_opt(id)  ! =1 or 2

  !write(*,*) "SADLESS tendency..."
  !write(*,*) "opt=",opt
  !write(*,*) "ips,ipe,jps,jpe",ips,ipe,jps,jpe
  



  hours = (grid%itimestep-1)*grid%dt/3600.  
  !write(*,*)'sadles_tendency: hours=', hours," <--",grid%itimestep,grid%dt
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

  !Get the ambient wind speed 

  !Loop for all turbine and calc. Thrust coeficient
  ntb = nturbines(id)
  th_coef(1:ntb,id) = 0.
  tke_coef(1:ntb,id) = 0.
  ai_fact(1:ntb,id) = 0.
  out_HubSpd(1:ntb,id,i_buffer) = 0.
  out_AmbSpd(1:ntb,id,i_buffer) = 0.
  out_Pow(1:ntb,id,i_buffer) = 0.


! First get the speed at hubheight
  do itb=1,ntb
      ir = tb_i(itb,id)
      jr = tb_j(itb,id)
      kr = amb_k(itb,id)
      ttype = tb_types(itb,id)
      if ( ( ir .le. ipe ) .and.  ( ir .ge. ips ) .and. (jr .le. jpe ) .and. (jr.ge. jps ) ) then
          out_HubSpd(itb,id,i_buffer) = sqrt( u(ir,kr,jr)**2 + v(ir,kr,jr)**2)
          ai_fact(itb,id)  = Axial_induction(id,itb,out_HubSpd(itb,id,i_buffer))
          if (opt.eq.2) then              
              th_coef(itb,id)  = Thrust_coef(id,itb,out_HubSpd(itb,id,i_buffer),2)             
              out_AmbSpd(itb,id,i_buffer) = out_HubSpd(itb,id,i_buffer)/(1-ai_fact(itb,id)) 
              out_Pow(itb,id,i_buffer) = turbine_power(id,itb,out_AmbSpd(itb,id,i_buffer))
              C_pow =  2.*out_Pow(itb,id,i_buffer) /(rho_std*tb_area(ttype)*out_AmbSpd(itb,id,i_buffer)**3 )
              tke_coef(itb,id) = ctke_fact*(th_coef(itb,id) - C_pow)
          end if
      end if 
  end do
  CALL MPI_ALLREDUCE(out_HubSpd(1:ntb,id,i_buffer), out_HubSpd(1:ntb,id,i_buffer),  ntb,MPI_REAL, MPI_SUM,  sadcom,ierror)

!for opt=1 only
  if (opt.eq.1) then
     do itb=1,ntb
        ir = amb_i(itb,id)
        jr = amb_j(itb,id)
        kr = amb_k(itb,id)
        ttype = tb_types(itb,id)
        if ( ( ir .le. ipe ) .and.  ( ir .ge. ips ) .and. (jr .le. jpe ) .and. (jr.ge. jps ) ) then
           out_AmbSpd(itb,id,i_buffer) = sqrt( u(ir,kr,jr)**2 + v(ir,kr,jr)**2)
           ai_fact_temp = ( out_AmbSpd(itb,id,i_buffer)  - out_HubSpd(itb,id,i_buffer) )/ out_AmbSpd(itb,id,i_buffer)
           if (ai_fact_temp .le. ai_fact(itb,id)) then
               ai_fact(itb,id) =  ai_fact_temp
           end if
   
           th_coef(itb,id)  = Thrust_coef(id,itb,out_AmbSpd(itb,id,i_buffer),1)
           out_Pow(itb,id,i_buffer) = turbine_power(id,itb,out_AmbSpd(itb,id,i_buffer))
           C_pow =  2.*out_Pow(itb,id,i_buffer) /(rho_std*tb_area(ttype)*out_AmbSpd(itb,id,i_buffer)**3 )
           tke_coef(itb,id) = ctke_fact*(th_coef(itb,id) - C_pow)
         end if
      end do 
  end if


  CALL MPI_ALLREDUCE(th_coef(1:ntb,id), th_coef(1:ntb,id),  ntb,MPI_REAL, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(tke_coef(1:ntb,id), tke_coef(1:ntb,id),  ntb,MPI_REAL, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(ai_fact(1:ntb,id), ai_fact(1:ntb,id),  ntb,MPI_REAL, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(out_AmbSpd(1:ntb,id,i_buffer), out_AmbSpd(1:ntb,id,i_buffer),  ntb,MPI_REAL, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(out_Pow(1:ntb,id,i_buffer), out_Pow(1:ntb,id,i_buffer),  ntb,MPI_REAL, MPI_SUM,  sadcom,ierror)
  !At this point, we have: Ct, Ctke, a

  !Check
  !write(*,*) "Turbine thrust, Ambient speed,  coeficient domain #",id
  !do itb=1,ntb
  !   write(*,*) "#",itb,"-->",th_coef(itb,id),out_AmbSpd(itb,id,i_buffer)
  !end do


  IF ( wrf_dm_on_monitor() ) THEN
     if (i_buffer.eq.out_buffer) then
         write(*,*) "Sadles: writing buffer to file"
         write(tb_num,'(I3)'),nturbines(id)
         fformat="(f12.7,"//trim(tb_num)//"(f6.2))"

         write(output,"('sadles_AmbSpd.d',i0.2)") id
         OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND',ACTION="WRITE")
         write(*,*)fformat
         do i_buffer=1,out_buffer
             write(foutput,fformat) out_hours(i_buffer), out_AmbSpd(1:ntb,id,i_buffer)
         enddo  
         CLOSE(foutput)

         write(output,"('sadles_HubSpd.d',i0.2)") id
         OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND',ACTION="WRITE")
         write(*,*)fformat
         do i_buffer=1,out_buffer
             write(foutput,fformat) out_hours(i_buffer), out_HubSpd(1:ntb,id,i_buffer)
         enddo
         CLOSE(foutput)

         fformat="(f12.7,"//trim(tb_num)//"(f8.2))"
         write(output,"('sadles_Power.d',i0.2)") id
         OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='OLD',POSITION='APPEND',ACTION="WRITE")
         write(*,*)fformat
         do i_buffer=1,out_buffer
             write(foutput,fformat) out_hours(i_buffer), out_Pow(1:ntb,id,i_buffer)
         enddo
         CLOSE(foutput)

         i_buffer = 0
   
     endif
  END IF
  CALL wrf_dm_bcast_integer(i_buffer,1)



  !Now update tendency for u,v
  do itb=1,ntb
      ir = tb_i(itb,id)
      jr = tb_j(itb,id)
      ttype = tb_types(itb,id)
      nrad = tb_nradius(ttype)
      ! If turbine is within the current memory tile:
      if ( ( ir-nrad .le. ipe ) .and.  ( ir+nrad .ge. ips ) .and. (jr-nrad .le. jpe ) .and. (jr+nrad.ge. jps ) ) then
         write(*,*)"Update tendency for turbine #",itb, th_coef(itb,id) 
         !th_coef(itb,id)  
         do ipt=1, sad_npoints(itb,id)
            ii = sad_i(itb,id,ipt)
            kk = sad_k(itb,id,ipt)
            jj = sad_j(itb,id,ipt)
            !Only update within the current tile
            if ( ( ii .le. ipe ) .and.  ( ii .ge. ips ) .and. (jj .le. jpe ) .and. (jj.ge.jps ) ) then
               speed = sqrt( u(ii,kk,jj)**2 + v(ii,kk,jj)**2)
               ru_tendf(ii,kk,jj) =  ru_tendf(ii,kk,jj) - 0.5*mut(ii,jj)*th_coef(itb,id)*speed* u(ii,kk,jj)*sad_AreaFact(itb,id,ipt) / (1 - ai_fact(itb,id))**2. 
               rv_tendf(ii,kk,jj) =  rv_tendf(ii,kk,jj) - 0.5*mut(ii,jj)*th_coef(itb,id)*speed* v(ii,kk,jj)*sad_AreaFact(itb,id,ipt) / (1 - ai_fact(itb,id))**2.
               tke_tend(ii,kk,jj) =  tke_tend(ii,kk,jj) + 0.5*mut(ii,jj)*tke_coef(itb,id)*(speed**3) *sad_AreaFact(itb,id,ipt) / (1 - ai_fact(itb,id))**3.
            end if
         end do
      end if
  end do



  !Stop for debug
  !stop
  END SUBROUTINE sadles_tendency



 !Calculate Thrust coef from table
 real function Thrust_coef(id,tbid,speed,opt)
     implicit none
     integer, intent(in) :: id,tbid,opt
     real, intent(in) :: speed
     REAL, DIMENSION(max_tb_steps) :: ref_speed, ref_coef
     integer :: steps,ttype
     ttype = tb_types(tbid,id) 
     ref_coef = tb_thcoef(ttype,:)
     if (opt.eq.1) then
         ref_speed = tb_ws(ttype,:)  
     else
         ref_speed = tb_ws1(ttype,:)  
     end if
     steps=tb_steps(ttype)
     Thrust_coef = lin_interp(ref_coef,ref_speed,steps,speed)     !Now we need to interpolate from these reference speed and coeficent      
  end function Thrust_coef
 
 !Calculate Thrust coef from table, opt=2 only
 real function Axial_induction(id,tbid,speed)
     implicit none
     integer, intent(in) :: id,tbid
     real, intent(in) :: speed
     REAL, DIMENSION(max_tb_steps) :: ref_speed, ref_coef
     integer :: steps,ttype
     ttype = tb_types(tbid,id)
     ref_coef = tb_aif(ttype,:)
     ref_speed = tb_ws1(ttype,:)  !Use hub speed 
     steps=tb_steps(ttype)
     Axial_induction = lin_interp(ref_coef,ref_speed,steps,speed)     !Now we need to interpolate from these reference speed and coeficent      
  end function Axial_induction
  
 !Calculate power
 real function turbine_power(id,tbid,speed)
     implicit none
     integer, intent(in) :: id,tbid
     real, intent(in) :: speed
     REAL, DIMENSION(max_tb_steps) :: ref_speed, ref_power
     integer :: steps,ttype
     ttype = tb_types(tbid,id)
     ref_power = tb_power(ttype,:)
     ref_speed = tb_ws(ttype,:)  !Use ambient speed 
     steps=tb_steps(ttype)
     turbine_power = lin_interp(ref_power,ref_speed,steps,speed)     
  end function turbine_power




!
! The simple method is just direct the turbine against the spontaneous wind
! direction
!
!
! Note: what if the turbine is near a tile border?
! Can we broadcast from a normal cpu?
!

  SUBROUTINE  sadles_yaw(id,grid, dx    &
              ,z3d,u,v                  &
              ,ims,ime,jms,jme,kms,kme  &
              ,ips,ipe,jps,jpe          &
              )
  Implicit None

  !TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
  TYPE(domain) , TARGET          :: grid
  INTEGER, INTENT(IN) :: id
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ips,ipe,jps,jpe
  REAL, INTENT(IN) :: dx
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(IN) :: z3d,u,v

  !This contain vertical discretized areas and radius
  !Note, tb_Rz is converted to index unit
  REAL, DIMENSION(kms:kme)  :: tb_Faz, tb_Rz,z

  integer :: k,nz, itb, ii,jj, ttype, kbot,khub,ktop
  real :: tb_z, tb_r, ir, jr, &
          uhub, vhub, uamb, vamb, &    ! uamb, vamb are abmbient velocity measusued 1 diameter in front of turbine
          angle,                  &    ! wind direction (not WMO convention)
          ti0,ti1,tj0,tj1              ! coordinates for turbline line


  integer :: row_n_points
  integer :: icell,ierror, ntb
  integer :: sadcom
  

  nz = kme-kms
  !write(*,*) "SADLESS Yawing..."
  !write(*,*) "Memory dimensions",ims,ime,jms,jme
  !write(*,*) "Patch dimensions",ips,ipe,jps,jpe
  !write(*,*) "nturbines=",nturbines(id)
  !write(*,*) "dx",dx


  

  !reinit value (needed for MPI_ALLREDUCE to work correctly)
  ntb = nturbines(id)
  sad_i(1:ntb,id,:)=0
  sad_k(1:ntb,id,:)=0
  sad_j(1:ntb,id,:)=0
  sad_AreaFact(1:ntb,id,:) = 0.
  sad_npoints(1:ntb,id) =0
  amb_i(1:ntb,id)  = 0
  amb_k(1:ntb,id)  = 0
  amb_j(1:ntb,id) = 0


  !Loop for all turbine 
  do itb=1,ntb
      ttype=tb_types(itb,id)
      ir = tb_i(itb,id)
      jr = tb_j(itb,id)
      tb_z = tb_hubheight(ttype)
      tb_r = tb_diameter(ttype)/2
      ! If the turbine is within the current memory tile:
      if ( ( ir .le. ipe+0.5 ) .and.  ( ir .gt. ips-0.5 ) .and. (jr .le. jpe+0.5 ) .and.  (jr.gt. jps-0.5 ) ) then
          !write(*,*) "Within Cell:Turbine #",itb,ir,jr
          z(:) = z3d(int(ir),:,int(jr))
          CALL vert_AreaFact(tb_z, tb_r, z, nz, dx, tb_Faz, tb_Rz,kbot,khub,ktop)
          uhub = u(int(ir),khub,int(jr))
          vhub = v(int(ir),khub,int(jr))
          angle = atan2(-vhub,-uhub)
               
          !write(*,*)"After calling vert_AreaFact: kbot,khub,ktop=",kbot,khub,ktop   
          !write(*,*) "Turbine vertical discretized area factors and normalized radius:"
          !write(*,*) "uhub,vhub,angle:",uhub,vhub,angle
          !do k=kbot,ktop
          !   write(*,*) k, tb_Faz(k), tb_Rz(k)
          !end do

         !Calculate ambient coorinate for calc. ambient wind speed and Thrust
         !coefiction later
          amb_k(itb,id) = khub
          amb_i(itb,id) = int(  ir + 4*tb_r*cos(angle)/dx )
          amb_j(itb,id) = int(  jr + 4*tb_r*sin(angle)/dx )
          !Notice: ambient i,j can be out of the tile!
          !write(*,*) "Ambient i,j,k:",amb_i(itb,id),amb_j(itb,id),amb_k(itb,id)
          !we need to broadcast the Cthrust once calcualted

          !Now calculate the 3D grid-based Area factor
          do k=kbot,ktop
             ti0 = ir - tb_Rz(k)*sin(-angle)
             ti1 = ir + tb_Rz(k)*sin(-angle)
             tj0 = jr - tb_Rz(k)*cos(angle)
             tj1 = jr + tb_Rz(k)*cos(angle)
             call horz_AreaSplit(tb_Faz(k),ti0,tj0,ti1,tj1,row_max_points,row_n_points,row_i,row_j,row_Area)
             do icell=1,row_n_points
                 sad_npoints(itb,id) = sad_npoints(itb,id) +1
                 if (sad_npoints(itb,id).gt.sad_max_points) then
                    write(*,*) "ERROR in sadles_yaw: sad_max_points is too small"
                    write(*,*) "sad_max_points=",sad_max_points
                    write(*,*) "Should increase sadless_maxradius or decrease sadles_mindx,sadles_mindz"
                    stop
                 end if
                 sad_i(itb,id,sad_npoints(itb,id)) = row_i(icell)
                 sad_j(itb,id,sad_npoints(itb,id)) = row_j(icell)
                 sad_k(itb,id,sad_npoints(itb,id)) = k
                 sad_AreaFact(itb,id,sad_npoints(itb,id)) = row_Area(icell)
             end do 
          end do
      endif
      
  end do

  !Check if the information is passed successfully
  !write(*,*)"Before MPI_ALLREDUCE...."
  !write(*,*)"Sadles_yaw: turbines in domain#",id
  !do itb=1,ntb
  !    !write(*,*) "Turbine#",itb,sad_npoints(itb,id),amb_i(itb,id),amb_k(itb,id), amb_j(itb,id)
  !    do icell=1, sad_npoints(itb,id)
  !        !write(*,*)icell,"-->",sad_i(itb,id,icell),sad_k(itb,id,icell),sad_j(itb,id,icell),sad_AreaFact(itb,id,icell)
  !    end do
  !end do



  !Pass the information of all turbines to all processes
  !this only work if there is no overlape of the processing
  CALL WRF_GET_DM_COMMUNICATOR(sadcom)
  CALL MPI_ALLREDUCE(sad_npoints(1:ntb,id), sad_npoints(1:ntb,id),  ntb, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(amb_i(1:ntb,id), amb_i(1:ntb,id), ntb, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(amb_j(1:ntb,id), amb_j(1:ntb,id), ntb, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(amb_k(1:ntb,id), amb_k(1:ntb,id), ntb, MPI_INTEGER, MPI_SUM,  sadcom,ierror)

  CALL MPI_ALLREDUCE(sad_i(1:ntb,id,:), sad_i(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_k(1:ntb,id,:), sad_k(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_j(1:ntb,id,:), sad_j(1:ntb,id,:), ntb*sad_max_points, MPI_INTEGER, MPI_SUM,  sadcom,ierror)
  CALL MPI_ALLREDUCE(sad_AreaFact(1:ntb,id,:), sad_AreaFact(1:ntb,id,:), ntb*sad_max_points, MPI_REAL, MPI_SUM,  sadcom,ierror)

  !Check if the information is passed successfully
  !write(*,*)"After MPI_ALLREDUCE...."
  !write(*,*)"Sadles_yaw: turbines in domain#",id
  !do itb=1,ntb
  !    write(*,*) "Turbine #",itb,sad_npoints(itb,id),amb_i(itb,id),amb_k(itb,id), amb_j(itb,id)
  !    do icell=1, sad_npoints(itb,id)
  !        write(*,*)icell,"-->",sad_i(itb,id,icell),sad_k(itb,id,icell),sad_j(itb,id,icell),sad_AreaFact(itb,id,icell)
  !    end do
  !end do
  !stop
  
  END SUBROUTINE sadles_yaw

 

  

!
! Note: Some code adapted from module_wind_fitch
!
!
  SUBROUTINE init_module_sadles(id,config_flags,xlong,xlat,sadles_initialized,&
                                            ims,ime,jms,jme,kms,kme,its,ite,jts,jte,ids,ide,jds,jde)
!
  IMPLICIT NONE
!
   integer id,ims,ime,jms,jme,kms,kme,ids,ide,jds,jde
   integer its,ite,jts,jte
   REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
   TYPE (grid_config_rec_type) :: config_flags
   TYPE (PROJ_INFO) :: ts_proj
   logical :: sadles_initialized
!  
   LOGICAL, EXTERNAL :: wrf_dm_on_monitor
   integer :: max_dom
   real :: max_radius,min_dx,min_dz
!local variables
   integer :: tmp_n
   real :: tmp_lat,tmp_lon, tmp_i, tmp_j,dx
   INTEGER i,j,k,i_tb, itype, ntb
   CHARACTER*256 tb_num,input,output,fformat
   real, dimension(max_tb_steps) :: aif  !axial induction factor 

!use the same option of fitch scheme
   write(*,*)"init_module_sadles"
  
   max_dom = config_flags%max_dom
   dx = config_flags%dx
   max_radius = config_flags%sadles_maxradius
   min_dx = config_flags%sadles_mindx
   min_dz = config_flags%sadles_mindz
   sadles_domains = max_dom - id + 1
   sad_max_points = ceiling( (3.14159*max_radius**2)/(min_dx*min_dz) )
   row_max_points = ceiling(max_radius/min_dx)*2
   write(*,*)"sad_max_points=",sad_max_points
   write(*,*)"row_max_points=",row_max_points
   IF ( wrf_dm_on_monitor() ) THEN
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
    END IF


    CALL wrf_dm_bcast_integer(max_turbines,1)

    write(*,*)"max_turbines=",max_turbines
!-- Allocate arrays for turbines list
!To save the parallel communication, 
! We don't need to allocde mememory for all domains (max_domains)
! But just from the first domain that have turbine id:max_domains)
! Also, it's neccessary to filter out turbines not needed to reduce max_turbines
    if (.not. sadles_initialized) then
        allocate (tb_types(max_turbines,id:max_dom), &
                  tb_i(max_turbines,id:max_dom), &
                  tb_j(max_turbines,id:max_dom), &
                  amb_i(max_turbines,id:max_dom), &
                  amb_j(max_turbines,id:max_dom), &
                  amb_k(max_turbines,id:max_dom),  &
                  amb_spd(max_turbines,id:max_dom),  &
                  sad_npoints(max_turbines,id:max_dom), &
                  th_coef(max_turbines,id:max_dom), &
                  tke_coef(max_turbines,id:max_dom), &
                  ai_fact(max_turbines,id:max_dom))

        allocate ( sad_i(max_turbines,id:max_dom,sad_max_points), & 
                   sad_k(max_turbines,id:max_dom,sad_max_points), &
                   sad_j(max_turbines,id:max_dom,sad_max_points), &
                   sad_AreaFact(max_turbines,id:max_dom,sad_max_points))

        allocate ( out_AmbSpd(max_turbines,id:max_dom,out_buffer), &
                   out_Pow(max_turbines,id:max_dom,out_buffer), &
                   out_HubSpd(max_turbines,id:max_dom,out_buffer) )
        allocate ( row_i(row_max_points), row_j(row_max_points), &
                   row_Area(row_max_points))


!        write(*,*)"sadles_initialized from domain id ",id,' -> sadles_domains=',sadles_domains
        tb_i=-9999
        tb_j=-9999
        sadles_initialized=.true.
    endif


!-- Read turbine list
   IF ( wrf_dm_on_monitor() ) THEN
       i_tb=1
       nturbines(id) = 0
       write(*,*) "Sadles_init: Turbines in domain #",id
       !write(*,*) " max_turbines=",max_turbines
       do i=1,max_turbines
          if ( config_flags%windfarm_ij .eq. 1 ) then
              read(70,*) tb_i(i,id), tb_j(i,id), tb_types(i,id)
              nturbines(id)=nturbines(id)+1
          else
              read(70,*)  tmp_lat,tmp_lon,tmp_n
              CALL latlon_to_ij(ts_proj, tmp_lat,tmp_lon,  tmp_i, tmp_j)
              if ((tmp_i .le. ide) .and. (tmp_i.ge.ids) .and. (tmp_j .le. jde) .and. (tmp_j.ge.jds)) then
                  !write(*,*)i_tb,tmp_lat,tmp_lon,tmp_n
                  tb_i(i_tb,id) = tmp_i
                  tb_j(i_tb,id) = tmp_j 
                  tb_types(i_tb,id) = tmp_n
                  i_tb = i_tb+1 
                  nturbines(id)=nturbines(id)+1
              end if 
          endif
       end do

       write(*,*) "Turbines in domain#",id,": ",nturbines(id),"/",max_turbines
       do i_tb=1,  nturbines(id)
          write(*,*)'SADLES Turbine #',i_tb,':, (i,j) = (',tb_i(i_tb,id),',', tb_j(i_tb,id),'); Type= ',tb_types(i_tb,id)
       end do
       close(70)

!-- Read turbine specifications for each type
       tb_steps=0
       do i_tb=1,nturbines(id)
          itype = tb_types(i_tb,id)
          if (itype.gt.max_tb_types) then
              write(*,*) "Error: turbine type must smaller than ", max_tb_types
             stop
          end if

          if (tb_steps(itype) .eq. 0 ) then    !Only read once!


             write(tb_num,*) itype
             tb_num=trim(adjustl(tb_num))
             input="wind-turbine-"//trim(tb_num)//".tbl"
 !            write(*,*)"Opening file:>",input,"<"
             OPEN(file=trim(input),unit=19,FORM='FORMATTED',STATUS='OLD')
             READ (19,*,ERR=132)tb_steps(itype)

             if (tb_steps(itype).gt.max_tb_steps) then
                write(*,*) "Error: turbine step must smaller than ",max_tb_steps
                stop
             end if

             read(19,*,ERR=132)tb_hubheight(itype),tb_diameter(itype),tb_sthcoef(itype),tb_npower(itype)
             tb_nradius(itype) = tb_diameter(itype)/dx
             tb_area(itype) = 3.14159265359*tb_nradius(itype)**2/4. 
             do k=1,tb_steps(itype)
                 read(19,*,ERR=132) tb_ws(itype,k), tb_thcoef(itype,k), tb_power(itype,k)
                 !Calculate windspeed at AD
                 tb_ws1(itype,k) = 0.5*tb_ws(itype,k)*(1 + SQRT( 1 -  tb_thcoef(itype,k)  ) )
                 tb_aif(itype,k) = 0.5*(1 - SQRT( 1 - tb_thcoef(itype,k) ) )
             enddo
             close (19)

!             write(*,*) "Turbine type#",itype,tb_hubheight(itype),tb_diameter(itype),tb_sthcoef(itype),tb_npower(itype),tb_steps(itype)
           end if
       end do
 132   continue
      
! Write turbine information to file

       write(output,"('sadles_info.d',i0.2)") id
       write(*,*) "Writing ",output

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

       write(output,"('sadles_AmbSpd.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE')
       CLOSE(foutput)             
       write(output,"('sadles_HubSpd.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE')
       CLOSE(foutput)             
       write(output,"('sadles_Power.d',i0.2)") id
       OPEN(file=output,unit=foutput,FORM='FORMATTED',STATUS='REPLACE')
       CLOSE(foutput)             
 
   END IF
!   write(*,*)"...Broadcasting data:",SIZE(tb_types)," ",max_turbines*sadles_domains
   

   CALL wrf_dm_bcast_integer(nturbines,max_domains)
   CALL wrf_dm_bcast_integer(tb_types,max_turbines*sadles_domains)
   CALL wrf_dm_bcast_real(tb_i,max_turbines*sadles_domains)
   CALL wrf_dm_bcast_real(tb_j,max_turbines*sadles_domains)

   CALL wrf_dm_bcast_real(tb_hubheight,max_tb_types)
   CALL wrf_dm_bcast_real(tb_diameter,max_tb_types)
   CALL wrf_dm_bcast_real(tb_area,max_tb_types)
   CALL wrf_dm_bcast_real(tb_nradius,max_tb_types)
   CALL wrf_dm_bcast_real(tb_npower,max_tb_types)
   CALL wrf_dm_bcast_real(tb_sthcoef,max_tb_types)
   CALL wrf_dm_bcast_real(tb_npower,max_tb_types)
   CALL wrf_dm_bcast_integer(tb_steps,max_tb_types)

   CALL wrf_dm_bcast_real(tb_ws,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_ws1,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_aif,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_thcoef,max_tb_types*max_tb_steps)
   CALL wrf_dm_bcast_real(tb_power,max_tb_types*max_tb_steps)

   i_buffer = 0
   sadles_opt(id)=config_flags%sadles_opt
  END SUBROUTINE init_module_sadles

!map projection
  SUBROUTINE setup_wrf_mapproj(config_flags,ts_proj,xlong,xlat, &
                                            ims,ime,jms,jme,its,ite,jts,jte,ids,ide,jds,jde)
       integer ims,ime,jms,jme,ids,ide,jds,jde
       integer its,ite,jts,jte
       TYPE (grid_config_rec_type) :: config_flags
       TYPE (PROJ_INFO) :: ts_proj
       REAL,     DIMENSION( ims:ime , jms:jme ) , INTENT(IN) :: xlong,xlat
       REAL :: known_lat, known_lon

         CALL map_init(ts_proj)

         known_lat = xlat(its,jts)
         known_lon = xlong(its,jts)

      ! Mercator
      IF (config_flags%map_proj == PROJ_MERC) THEN
         CALL map_set(PROJ_MERC, ts_proj,               &
                      truelat1 = config_flags%truelat1, &
                      lat1     = known_lat,             &
                      lon1     = known_lon,             &
                      knowni   = REAL(its),             &
                      knownj   = REAL(jts),             &
                      dx       = config_flags%dx)

      ! Lambert conformal
      ELSE IF (config_flags%map_proj == PROJ_LC) THEN
         CALL map_set(PROJ_LC, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      truelat2 = config_flags%truelat2,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = REAL(its),              &
                      knownj   = REAL(jts),              &
                      dx       = config_flags%dx)
!      ! Polar stereographic
      ELSE IF (config_flags%map_proj == PROJ_PS) THEN
         CALL map_set(PROJ_PS, ts_proj,                  &
                      truelat1 = config_flags%truelat1,  &
                      stdlon   = config_flags%stand_lon, &
                      lat1     = known_lat,              &
                      lon1     = known_lon,              &
                      knowni   = REAL(its),              &
                      knownj   = REAL(jts),              &
                      dx       = config_flags%dx)
      END IF
  END SUBROUTINE setup_wrf_mapproj

  
 !Calculate the circle segment area for z to the top of the cirle (z is
 !relative to the circle origo)
 real function circle_segment(r,z)
     implicit none
     real, intent(in) :: r,z
     real :: area
     real :: r1,ang
     if (abs(z).gt.r) then
        area = 0.
     else
        r1 = sqrt(r**2 - z**2)
        ang = acos(z/r)
        area = ang*r**2 - z*r1
     end if
     circle_segment=area
     return 
  end function circle_segment



  !Calculate between veritcal level
  !tb_z, tb_r are turbine hubheight and radius, z are height of full-level,
  !which is assume to change very little, we can recalculate when yawing as well
  !areas are between the vertical levels
  !Rs are the Radius of the approximated segment at each half level
  !  then normalized by dx


  subroutine  vert_AreaFact(tb_z,tb_r, z, nz, dx, Fas,Rs,kbot,khub,ktop)
     implicit none
     real, intent(in) :: tb_z, tb_r, dx
     integer, intent(in) :: nz
     real, dimension(nz),intent(in) :: z
     real, dimension(nz-1), intent(out) :: Fas, Rs
     integer, intent(out) :: kbot,khub,ktop
 
     real :: t_zbot,t_ztop, zr1,zr2,dx2, area1, area2
     real, dimension(nz-1) :: z1s, z2s, areas
     
     integer :: k     

     Fas=0.
     Rs=0.
     t_zbot = tb_z - tb_r
     t_ztop = tb_z + tb_r
     dx2 = dx*dx
     kbot=-1
     khub=-1
     ktop=-1
     do k=1,nz
        if ( (z(k).ge.t_zbot) .and. (kbot.eq.-1) ) kbot = k
        if ( (z(k).ge.tb_z) .and. (khub.eq.-1) ) khub = k
        if ( (z(k).ge.t_ztop) .and. (ktop.eq.-1) ) then
           ktop = k-1
           exit
        end if
     end do
     z1s=0.
     z2s=0.
     do k=kbot,ktop
         z1s(k) = z(k)
         z2s(k) = z(k+1)
         if (k.eq.kbot) then
            z1s(k) = t_zbot
         else if (k.eq.ktop) then
            z2s(k) = t_ztop
         end if
         area1 = circle_segment(tb_r,z1s(k)-tb_z)
         area2 = circle_segment(tb_r,z2s(k)-tb_z)
         areas(k) = area1 - area2
         Fas(k) = areas(k) / (z(k+1) - z(k)) / dx2
         Rs(k) = areas(k) / (z2s(k)-z1s(k)) / 2
     end do 
     Rs = Rs/dx  !Normalized by dx

  end subroutine vert_AreaFact

  !Discretiziing Area into horizontal grid
  !return to number of grid point (n_points), their associate i,j, and
  !splited
  !area
  subroutine horz_AreaSplit(Area,ir0,jr0,ir1,jr1,max_points,n_points,grid_i,grid_j,grid_Area)
     implicit none
     real,intent(in) :: Area,ir0,jr0,ir1,jr1
     integer,intent(in) :: max_points
     integer,intent(out) :: n_points
     integer,dimension (max_points), intent(out) :: grid_i,grid_j
     real,dimension (max_points), intent(out)  :: grid_Area

     integer, parameter :: n_seq = 50
     real, dimension(n_seq) :: ir,jr
     integer :: ip,i,j,curr_i,curr_j,curr_n

     call linspace(ir0,ir1,n_seq,ir)
     call linspace(jr0,jr1,n_seq,jr)
     curr_i=-1
     curr_j=-1
     n_points=0
     grid_Area = 0.
     do ip=1,n_seq
        i=int(ir(ip))
        j=int(jr(ip))
        if ((i.ne.curr_i).or.(j.ne.curr_j)) then  ! a new cell
           n_points=n_points+1
           grid_i(n_points) = i
           grid_j(n_points) = j
           curr_i = i
           curr_j = j
        end if
        grid_Area(n_points) =  grid_Area(n_points) + Area/n_seq
     end do
  end subroutine horz_AreaSplit

 !npoint>=2
  subroutine linspace(xstart,xend,npoint,arr)
     implicit none
     real, intent(in) :: xstart,xend
     integer, intent(in) :: npoint
     real, dimension(npoint),intent(out) :: arr
     integer :: i
     real :: dx
     dx = (xend-xstart)/(npoint-1)
     arr = (/(xstart+(i-1)*dx, i=1,npoint)/)
  end subroutine linspace

!simple linear interpolation
  real function lin_interp(Y,X,Nx,x0)
   implicit none
   integer, intent(In) :: Nx
   real, dimension(Nx), intent(In) :: X, Y
   real, intent(In) :: x0
   integer :: i1,i2
   if (x0 .le.X(1)) then
      !lin_interp=Y(1) 
      lin_interp = 0 ! cut in speed
   else if (x0 .ge.X(Nx)) then
      !lin_interp=Y(Nx)
      lin_interp = 0 ! cut out speed
   else
      do i2=2,Nx        !Simple search, would be fast if using, for example binary search
         if (X(i2).ge.x0) then
            i1 = i2-1
            exit
         end if
      end do
      lin_interp = Y(i1) + (x0 - X(i1))*( Y(i2)-Y(i1) ) / ( X(i2)-X(i1) )
   end if
   return
   end function


END MODULE module_sadles



