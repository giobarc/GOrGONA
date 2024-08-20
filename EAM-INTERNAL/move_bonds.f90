SUBROUTINE MOVE_BONDS(eps)

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
  Integer :: ind,ind1,ind2,iatom,i_chosen,npvicmin,natomi_bonds,natomvmin
  Integer :: iposition(nat3d),npvic(nat3d),ivatom(nat3d)
  Real(8) :: deltarho,deltatheta,deltaphi
  Real(8) :: xcm,ycm,zcm
  Real(8) :: distcm,maxdist,pvic4,nnmax2
  Real(8) :: x_old,y_old,z_old,displa
  Real(8) :: dij,eps


  pvic4=1.d0
 
  nnmax2=dmax1(nn(1),nn(2))*dmax1(nn(1),nn(2))
!--------------------------------------------------------------------
!  !calculation of npvicmin
  !
  npvic(:)=0
  do ind1=1,nat3d
     do ind2=ind1+1,nat3d
        dij=(xin(ind1)-xin(ind2))*(xin(ind1)-xin(ind2))+&
             (yin(ind1)-yin(ind2))*(yin(ind1)-yin(ind2))+&
             (zin(ind1)-zin(ind2))*(zin(ind1)-zin(ind2))
        !         if ((ind1.eq.1)) write (*,*) 'distanza',ind2,sqrt(dij)
        if (dij.le.nnmax2) then
           npvic(ind1)=npvic(ind1)+1
           npvic(ind2)=npvic(ind2)+1
        endif
     enddo
  enddo
  npvicmin=100
  do ind1=1,nat3d
     if (npvic(ind1).lt.npvicmin) npvicmin=npvic(ind1)
  enddo
! list of atoms with npvicmin
  natomvmin=0
  do ind1=1,nat3d
     if (npvic(ind1).eq.npvicmin) then
        natomvmin=natomvmin+1
        ivatom(natomvmin)=ind1
     endif 
  enddo
  call ranmar(zzz)
  i_chosen=int(zzz*natomvmin)+1 
  i_chosen=min(natomvmin,i_chosen)
  iatom=ivatom(i_chosen)
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
  call ranmar(zzz)
  deltarho=(0.8d0+0.2d0*zzz)*maxdist
  call ranmar(zzz)
  deltatheta=acos(1.-2.*zzz)
  call ranmar(zzz)
  deltaphi=2*pi*zzz
  xin(iatom) = xcm+deltarho*sin(deltatheta)*cos(deltaphi)
  yin(iatom) = ycm+deltarho*sin(deltatheta)*sin(deltaphi)
  zin(iatom) = zcm+deltarho*cos(deltatheta)

  
END SUBROUTINE MOVE_BONDS
