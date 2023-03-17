!WRF:MiDEL_LAYER:PHYSICS

MODULE module_cellpert
!
! A simple cell perturbation method
!
!                   By Hai Bui, hai.bui@uib.no
!
!
  USE module_driver_constants, ONLY : max_domains
  USE module_model_constants, ONLY :  piconst
  USE module_domain
!
  USE module_configure, ONLY : grid_config_rec_type

  IMPLICIT NONE
  
  INTEGER, parameter :: cpert_size = 8
  INTEGER, parameter :: cpert_border = 24
  INTEGER, parameter :: Max_Nx_cpert=150,Max_Nz_cpert=40,Max_Ny_cpert=150  ! Equivalant Max_Nx,Max_Ny=1600,
  REAL, DIMENSION(Max_Nx_cpert,Max_Nz_cpert,Max_Ny_cpert) :: cpert  ! dims: max_doms, Ny_cpert, Nx_cpert 
            
                   
#if defined(DM_PARALLEL) && !defined(STUBMPI)
      INCLUDE "mpif.h"
#endif

  ! Note:
  !
  ! the dry air collum mass are grid%muu and grid%muv at u and v grid points
  !

CONTAINS


subroutine gen_th_cpert(config_flags,cell_pert_th, t_2, &
               ids,  ide,  jds,  jde,  kds,  kde,   &
               ims,ime, jms,jme, kms,kme,           &
               ips,ipe,jps,jpe)

  TYPE (grid_config_rec_type) :: config_flags
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_2  ! This is WRF grid! 
  REAL, DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(OUT) :: cell_pert_th  ! This is WRF grid! 
  INTEGER, INTENT(IN) :: ids,ide,jds,jde,kds,kde
  INTEGER, INTENT(IN) :: ims,ime,jms,jme,kms,kme
  INTEGER, INTENT(IN) :: ips,ipe,jps,jpe

  INTEGER :: k1, k2, i, j, k, ip,jp
  LOGICAL, EXTERNAL :: wrf_dm_on_monitor

  !calculate 
  INTEGER :: Nx_cpert, Nz_cpert, Ny_cpert
  real :: kweight
  
  k1 = config_flags%cell_pert_k1
  k2 = config_flags%cell_pert_k2

!First generate the random perturbation
  IF ( wrf_dm_on_monitor() ) THEN
     Nx_cpert = CEILING((ide-ids+1.)/cpert_size)
     Ny_cpert = CEILING((jde-jds+1.)/cpert_size)
     call RANDOM_NUMBER(cpert(1:Nx_cpert,1:k2,1:Ny_cpert))      
     cpert(1:Nx_cpert,1:k2,1:Ny_cpert) = 2*( cpert(1:Nx_cpert,1:k2,1:Ny_cpert) - 0.5) * config_flags%cell_pert_magnitude
  END IF
  CALL wrf_dm_bcast_integer(Nx_cpert,1)
  CALL wrf_dm_bcast_integer(Ny_cpert,1)
  !write(*,*) "Nx_cpert, Ny_cpert:",Nx_cpert, Ny_cpert
  !write(*,*) "Broadcast cpert.... ",Nx_cpert*k2*Ny_cpert
  CALL wrf_dm_bcast_real(cpert(1:Nx_cpert,1:k2,1:Ny_cpert),Nx_cpert*k2*Ny_cpert)

!Now assign cpert to cell_pert_th
  DO i=ips,ipe
    DO j=jps,jpe
       IF (  ((i .le. cpert_border)     .and. (1.eq.config_flags%cell_pert_xs)) .or. &
             ((i .ge. ide-cpert_border) .and. (1.eq.config_flags%cell_pert_xe)) .or. &
             ((j .le. cpert_border)     .and. (1.eq.config_flags%cell_pert_ys)) .or. &
             ((j .ge. jde-cpert_border) .and. (1.eq.config_flags%cell_pert_ye))  ) THEN

          ip = ceiling(1.*i/cpert_size)
          jp = ceiling(1.*j/cpert_size)
          DO k=1,k1-1
              cell_pert_th(i,k,j) = cpert(ip,k,jp)
              t_2(i,k,j) = t_2(i,k,j) + cpert(ip,k,jp)
          END DO
          DO k=k1,k2
              kweight = cos(0.5*piconst*(k-k1)/(k2-k1))**2
              cell_pert_th(i,k,j) = cpert(ip,k,jp)*kweight
              t_2(i,k,j) = t_2(i,k,j) + cpert(ip,k,jp)*kweight
          END DO
           
       END IF
    END DO
  END DO

end subroutine gen_th_cpert

END MODULE module_cellpert

