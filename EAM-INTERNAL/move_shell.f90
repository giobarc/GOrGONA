SUBROUTINE MOVE_SHELL(eps)

  !***************************
  !This move consists in:
  !1) randomly choosing an atom (iatom) of the cluster, among the less coordinated (number of nn < 9)
  !2) locating the maximum distance between atoms and cdm (maxdist)
  !3) moving iatom in a randomly chosen position within a shell of fixed thickness (shell_thickness=1.5 A)
  !centered in cdm and of minimum radius=maxdist
   !***************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  !local variables
  Integer :: ind1,ind,iatom
  Integer :: pochivicini,pochivicini_counter,pochivicini_array(nat3d)
  Integer :: aa
  Real(8) :: deltarho,deltatheta,deltaphi
  Real(8) :: xcm,ycm,zcm,eps
  Real(8) :: distcm,maxdist
  Real(8) :: x_old,y_old,z_old
  
  !evaluation of the number of nn of each atom:
  call comparison(pochivicini)
  !looking for atoms with less than 9 nn:
  pochivicini_counter=0
  
  do ind=1,nat3d
   if(quantivicini(ind) .lt. 9)then
   pochivicini_counter=pochivicini_counter+1
   pochivicini_array(pochivicini_counter)=ind
   endif
  enddo
      
  !random extraction of one atom between the less bonded atoms: 
  !this will be the atom to displace
  call ranmar(zzz)
  aa=int(pochivicini_counter*zzz)+1
  iatom=pochivicini_array(min(aa,pochivicini_counter))
  
  
  !saving old position of iatom
  x_old=xin(iatom)
  y_old=yin(iatom)
  z_old=zin(iatom)
  
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
  
  ! choose a position randomly within a shell shell_thickness thick, 
  !whose maximum radius is maxdist
      
     !a move is attempted to:
     call ranmar(zzz)
     deltarho=maxdist-zzz*shell_thickness
     call ranmar(zzz)
     deltatheta=pi*zzz
     call ranmar(zzz)
     deltaphi=2*pi*zzz
     
 
     xin(iatom) = xcm+deltarho*sin(deltatheta)*cos(deltaphi)
     yin(iatom) = ycm+deltarho*sin(deltatheta)*sin(deltaphi)
     zin(iatom) = zcm+deltarho*cos(deltatheta)
     
     
  
END SUBROUTINE MOVE_SHELL
