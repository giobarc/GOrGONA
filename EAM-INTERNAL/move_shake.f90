SUBROUTINE MOVE_SHAKE(eps)
  
  !***************************
  !This move consists in displacing every atom of the cluster
  !within a spherical shell centered in its initial position
  !whose minimum and maximum radii are written in the input file
  !***************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  !local variables
  Integer :: ind
  Real(8) :: deltarho,deltatheta,deltaphi
  Real*8 :: eps
  
  do ind=1,nat3d
!	zin(ind) = 0.0
!     	do while(zin(ind) < 0.005) 
     		!choice of the displacement
     		call ranmar(zzz)
     		deltarho=rho_inf_sk+zzz*(rho_sup_sk-rho_inf_sk)
     		!ATTENZIONE: tolgo riga che non capisco
     		deltarho=(rho_inf_sk*rho_inf_sk*rho_inf_sk+ (deltarho-rho_inf_sk)*&
          		(rho_sup_sk*rho_sup_sk+rho_sup_sk*rho_inf_sk+rho_inf_sk*rho_inf_sk))**(1./3.)
	
     		call ranmar(zzz)
!     		deltatheta=pi*zzz
     		deltatheta=acos(1.-2.*zzz)
     		call ranmar(zzz)
     		deltaphi=2*pi*zzz


     		xin(ind) = xin(ind)+deltarho*sin(deltatheta)*cos(deltaphi)
     		yin(ind) = yin(ind)+deltarho*sin(deltatheta)*sin(deltaphi)
     		zin(ind) = zin(ind)+deltarho*cos(deltatheta)
		
!	enddo     
  enddo
  
  
END SUBROUTINE MOVE_SHAKE
