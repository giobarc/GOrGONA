SUBROUTINE MOVE_BROWNIAN
  
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
  Integer :: ind,ipas_bro,i,j
  Integer :: ndc2(nat3d)
  logical :: first
  Real(8) :: vxx(nat3d),vyy(nat3d),vzz(nat3d),xv(nat3d),yv(nat3d),zv(nat3d)
  Real(8) :: dij2,dxij,dyij,dzij
  Real(8) :: cbltz,mass_bro,fric,tstep_bro,fact_bro_1,fact_bro_2,factv

  cbltz=1.38d-23

  fric=1.d12
  tstep_bro=5d-15
  mass_bro=1.67d-25
! tutte le quantità sono in unità del sistema internazionale, la massa è unica e di riferimento (un pò meno della
! massa dell'argento)
  first = .true.

  
  factv=dsqrt(cbltz*temp_brow/mass_bro)
  fact_bro_1=dsqrt(2.d0*cbltz*temp_brow*fric*tstep_bro/(mass_bro))*1d10 !Angs/sec
  fact_bro_2=(16.d0/arete(1))*tstep_bro/mass_bro
!  write (*,*) fact_bro_1,fact_bro_2
!  stop
  do ind=1,nat3d
     call gauss
     vxx(ind)=g1*factv
     vyy(ind)=g2*factv
     vzz(ind)=g3*factv
  enddo
  call bigvoi(first)
  first = .false.
  do ipas_bro=1,npas_bro
     do ind=1,nat3d
        xat(ind)=xin(ind)/arete(1)
        yat(ind)=yin(ind)/arete(1)
        zat(ind)=zin(ind)/arete(1)
     enddo
        call bigvoi(first)
     	call force_rgl
!	if(potential .eq. '_lj') call force_lj
     do ind=1,nat3d
        call gauss
!        write (*,*)g1,g2,g3
        vxx(ind)=vxx(ind)+fact_bro_1*g1-fric*tstep_bro*vxx(ind)+fact_bro_2*fx(ind)
	vyy(ind)=vyy(ind)+fact_bro_1*g2-fric*tstep_bro*vyy(ind)+fact_bro_2*fy(ind)
        vzz(ind)=vzz(ind)+fact_bro_1*g3-fric*tstep_bro*vzz(ind)+fact_bro_2*fz(ind)

        xin(ind)=xin(ind)+vxx(ind)*tstep_bro
        yin(ind)=yin(ind)+vyy(ind)*tstep_bro
        zin(ind)=zin(ind)+vzz(ind)*tstep_bro

     enddo    
  enddo

  
END SUBROUTINE MOVE_BROWNIAN
