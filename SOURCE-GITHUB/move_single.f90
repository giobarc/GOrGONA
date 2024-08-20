SUBROUTINE MOVE_SINGLE(eps)

  !***************************
  !This move consists in displacing one atom (randomly chosen) of the cluster
  !within a spherical shell centered in its initial position
  !whose minimum and maximum radii are written in the input file
  !***************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  !local variables
  Integer :: iatom
  Real(8) :: deltarho,deltatheta,deltaphi
  Real*8 :: eps
  
  !choice of the atom to displace
  call ranmar(zzz)
  iatom=nat3d*zzz+1
  iatom=min(iatom,nat3d)
  
  !choice of the displacement
  call ranmar(zzz)
  deltarho=rho_inf_sk+zzz*(rho_sup_sg-rho_inf_sg)
  call ranmar(zzz)
  deltatheta=pi*zzz
  call ranmar(zzz)
  deltaphi=2*pi*zzz

  if( zin(iatom)+deltarho*cos(deltatheta) <= 0 ) deltarho=(eps-zin(iatom))/cos(deltatheta)
  
  xin(iatom) = xin(iatom)+deltarho*sin(deltatheta)*cos(deltaphi)
  yin(iatom) = yin(iatom)+deltarho*sin(deltatheta)*sin(deltaphi)
  zin(iatom) = zin(iatom)+deltarho*cos(deltatheta)
  
END SUBROUTINE MOVE_SINGLE
