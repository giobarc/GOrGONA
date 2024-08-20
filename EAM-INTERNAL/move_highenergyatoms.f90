SUBROUTINE MOVE_HIGHENERGYATOMS(eps)

  !***************************
  !This move consists in:
  !1) evaluating the minimum number of pv in the cluster
  !2) if this number is .ge. 5, the shake move is applied
  !3) if this number is .lt. 5, there are two possibilities:
  !4a) npvicmin .le. 3 --> all atoms with npvic=3 are displaced
  !according to the cdm-centered move. 
  !4b) npvicmin .eq. 4 --> atoms with npvic=4 are displaced
  !with the probability pvic4 according to the cdm-centered move.
  !***************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  !local variables
  Integer :: ind,ind1,tentativi
  Real(8) :: deltarho,deltatheta,deltaphi
  Real(8) :: xcm,ycm,zcm,eps
  Real(8) :: distcm,maxdist,pvic4
  Real(8) :: x_old,y_old,z_old
  Real(8) :: ener_threshold
  
  pvic4=0.25
  ener_threshold=-2.2 !atoms with ener > ener_threshold are considered as high-energy atoms


  !------------------------------------------------------------------------------
  !move bonds if there are atoms with high energy (ener_atom>-2.2 eV)
  !move is performed on the weakly bonded atoms with 0.5 probability
  !probability pvic4
  do ind1=1,nat3d 
  call ranmar(zzz)
  if ((ener_atom(ind1) .gt. ener_threshold).and.(zzz.lt.0.5)) then
     x_old=xin(ind1)
     y_old=yin(ind1)
     z_old=zin(ind1)
     !calculate the center of mass
     xcm=0._8
     ycm=0._8
     zcm=0._8
     do ind=1,nat3d
        xcm=xcm+xin(ind)
        ycm=ycm+yin(ind)
        zcm=zcm+zin(ind)
     enddo
     xcm=xcm/float(nat3d)
     ycm=ycm/float(nat3d)
     zcm=zcm/float(nat3d)
     ! calculate the maximal distance from the center of mass
     maxdist=0.d0
     do ind=1,nat3d
        distcm=(xin(ind)-xcm)*(xin(ind)-xcm)+&
             (yin(ind)-ycm)*(yin(ind)-ycm)+&
             (zin(ind)-zcm)*(zin(ind)-zcm)
        if (distcm.gt.maxdist) maxdist=distcm
     enddo
     maxdist=sqrt(maxdist)
     ! choose a position randomly within a sphere of radius maxdist
        
        call ranmar(zzz)
        deltarho=(0.7+(0.3)*zzz)*maxdist
        call ranmar(zzz)
        deltatheta=pi*zzz
        call ranmar(zzz)
        deltaphi=2*pi*zzz
        

        xin(ind1) = xcm+deltarho*sin(deltatheta)*cos(deltaphi)
        yin(ind1) = ycm+deltarho*sin(deltatheta)*sin(deltaphi)
        zin(ind1) = zcm+deltarho*cos(deltatheta)
                
  !---------------------------------------------------------------------------------
  !move shake if ind1 is a high-energy atom:     
  else ! su atomi con bassa energia applico la mossa shake
        !choice of the displacement
        call ranmar(zzz)
        deltarho=rho_inf_he+zzz*(rho_sup_he-rho_inf_he)
        deltarho=(rho_inf_he*rho_inf_he*rho_inf_he+ (deltarho-rho_inf_he)*&
             (rho_sup_he*rho_sup_he+rho_sup_he*rho_inf_he+rho_inf_he*rho_inf_he))**(1./3.)
        call ranmar(zzz)
        deltatheta=pi*zzz
        call ranmar(zzz)
        deltaphi=2*pi*zzz
	
	if( zin(ind1)+deltarho*cos(deltatheta) <= 0 ) deltarho=(eps-zin(ind1))/cos(deltatheta)

        xin(ind1) = xin(ind1)+deltarho*sin(deltatheta)*cos(deltaphi)
        yin(ind1) = yin(ind1)+deltarho*sin(deltatheta)*sin(deltaphi)
        zin(ind1) = zin(ind1)+deltarho*cos(deltatheta)
 
  endif !su atomi con alta energia
  enddo
  
END SUBROUTINE MOVE_HIGHENERGYATOMS
