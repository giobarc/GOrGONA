SUBROUTINE MOVE_BALL(eps)

  !***************************
  !This move consists in:
  !1) randomly choosing an atom (iatom) of the cluster
  !2) locating the maximum distance between atoms and cdm (maxdist)
  !3) moving iatom in a randomly chosen position within a sphere
  !centered in cdm and of radius=maxdist
  !***************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  !local variables
  Integer :: iatom,ind
  Real(8) :: deltarho,deltatheta,deltaphi
  Real(8) :: xcm,ycm,zcm,eps
  Real(8) :: distcm,maxdist
  Real(8) :: x_old,y_old,z_old
      
  !random choice of the atom to displace
  call ranmar(zzz)
  iatom=int(nat3d*zzz)+1
  iatom=min(iatom,nat3d)
  
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
  
  ! choose a position randomly within a sphere of radius maxdist
      call ranmar(zzz)
     deltarho=zzz*maxdist
     call ranmar(zzz)
     deltatheta=pi*zzz
     call ranmar(zzz)
     deltaphi=2*pi*zzz
     
     xin(iatom) = xcm+deltarho*sin(deltatheta)*cos(deltaphi)
     yin(iatom) = ycm+deltarho*sin(deltatheta)*sin(deltaphi)
     zin(iatom) = zcm+deltarho*cos(deltatheta)
     
     
  
END SUBROUTINE MOVE_BALL
