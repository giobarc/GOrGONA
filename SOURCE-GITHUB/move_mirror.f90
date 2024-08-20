SUBROUTINE MOVE_MIRROR
  
  !***************************
  !This move consists rotating randomly the cluster and cutting upper and lower slices which
  !are specularly reflected
  !***************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  !local variables
  Integer :: ind,isotto,isopra
  Real(8) :: alphaang,betaang,gammaang,xff,yff,zff,zmin,zmax
  Real(8) :: xcenter,ycenter,zcenter,zinf,zsup

! center the cluster in zero
xcenter=0.d0
ycenter=0.d0
zcenter=0.d0
do ind=1,nat3d
   xcenter=xcenter+xin(ind)
   ycenter=ycenter+yin(ind)
   zcenter=zcenter+zin(ind)
enddo
xcenter=xcenter/dfloat(nat3d)
ycenter=ycenter/dfloat(nat3d)
zcenter=zcenter/dfloat(nat3d)
do ind=1,nat3d
   xin(ind)=xin(ind)-xcenter
   yin(ind)=yin(ind)-ycenter
   zin(ind)=zin(ind)-zcenter
enddo

!rotate the cluster randomly
 !choose the angles
 call ranmar(zzz)
 alphaang=2*pi*zzz
 call ranmar(zzz)
 betaang=2*pi*zzz
 call ranmar(zzz)
 gammaang=2*pi*zzz

!rotate all atoms

do ind=1,nat3d
   !rotate about z
   xff=xin(ind)
   yff=yin(ind)
   xin(ind)=dcos(alphaang)*xff-dsin(alphaang)*yff
   yin(ind)=dsin(alphaang)*xff+dcos(alphaang)*yff
   !rotate about y
   xff=xin(ind)
   zff=zin(ind)
   xin(ind)=dcos(betaang)*xff+dsin(betaang)*zff
   zin(ind)=-dsin(betaang)*xff+dcos(betaang)*zff
   !rotate about x
   yff=yin(ind)
   zff=zin(ind)
   yin(ind)=dcos(gammaang)*yff-dsin(gammaang)*zff
   zin(ind)=dsin(gammaang)*yff+dcos(gammaang)*zff
enddo

!choose zcut randomly

if (zcutflag.eq.1) then    
   call ranmar(zzz) 
   zcut=zzz
endif


!calculate zmin and zmax
!zmin=0.d0
zmax=0.d0
!deltaz=0.d0
do ind=1,nat3d
   if (dabs(zin(ind)).gt.zmax) zmax=dabs(zin(ind))   
enddo
zsup=(1.d0-zcut)*zmax
!zinf=(1.d0-zcut)*zmin
isotto=0
isopra=0
!specular reflection of two slices
do ind=1,nat3d 
   if (zin(ind).lt.-zsup) then 
      zin(ind)=-zin(ind)+deltaz
      isotto=isotto+1
   else     
      if (zin(ind).gt.zsup) then
         zin(ind)=-zin(ind)-deltaz
         isopra=isopra+1
      endif
   endif        
enddo
write (44,'(i8,3f12.6,2i6)') ind_tot,zmax,zcut,zsup,isotto,isopra
!write (55,'(3e14.6)') alphaang,betaang,gammaang
 flush(44)
END SUBROUTINE MOVE_MIRROR
