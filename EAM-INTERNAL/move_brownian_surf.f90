SUBROUTINE MOVE_BROWNIAN_SURF
  
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
  Integer :: ind,ind1,ind2,ipas_bro,i,j,nsurf,npviclim
  Integer :: ndc2(nat3d)
  Integer :: npvic(nat3d)
  logical :: first
  Real(8) :: vxx(nat3d),vyy(nat3d),vzz(nat3d),xv(nat3d),yv(nat3d),zv(nat3d)
  Real(8) :: dij2,dij,dxij,dyij,dzij,nnmax2
  Real(8) :: cbltz,mass_bro,fric,tstep_bro,fact_bro_1,fact_bro_2,factv

  cbltz=1.38d-23
  fric=1.d12
  tstep_bro=5d-15
  mass_bro=1.67d-25
! tutte le quantità sono in unità del sistema internazionale, la massa è unica e di riferimento (un pò meno della
! massa dell'argento)
  first = .true.

  
  factv=dsqrt(cbltz*temp_browsurf/mass_bro)
  fact_bro_1=dsqrt(2.d0*cbltz*temp_brow*fric*tstep_bro/(mass_bro))*1d10 !Angs/sec
  fact_bro_2=(16.d0/arete(1))*tstep_bro/mass_bro
!  write (*,*) fact_bro_1,fact_bro_2
!  stop
  nnmax2=dmax1(nn(1),nn(2))*dmax1(nn(1),nn(2))
!conta dei primivicini
  npviclim=11
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
  nsurf=0
  do ind=1,nat3d
     if (npvic(ind).lt.npviclim) then
     call gauss
     vxx(ind)=g1*factv
     vyy(ind)=g2*factv
     vzz(ind)=g3*factv
     nsurf=nsurf+1
     endif
  enddo
!  write (*,*) "nsurf=",nsurf
  call bigvoi(first)
  first = .false.
! calcolo numero vicini
 do ind=1,nat3d
    do ind2=ind+1,nat3d
       ddd=(xat(ind)-xat(ind1))
    enddo
 enddo
  do ipas_bro=1,npas_brosurf
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
        if (npvic(ind).lt.npviclim) then
        xin(ind)=xin(ind)+vxx(ind)*tstep_bro
        yin(ind)=yin(ind)+vyy(ind)*tstep_bro
        zin(ind)=zin(ind)+vzz(ind)*tstep_bro
        endif
     enddo    
  enddo
!        do ind=1,nat3d
!        write (*,*)spec(ind),xin(ind),yin(ind),zin(ind)
!        enddo
!  stop 
!     do ind=1,nat3d
!        ndc2(ind)=0
!     enddo
!      do i=1,nat3d-1
!        do j=i+1,nat3d
!           dxij=xv(j)-xv(i)
!           dyij=yv(j)-yv(i)
!           dzij=zv(j)-zv(i)
!
!           dij2=dxij*dxij+dyij*dyij+dzij*dzij
!
!           if (dij2.le.dc2*dc2) then
!              ndc2(i)=ndc2(i)+1
!              ndc2(j)=ndc2(j)+1
!           endif
!           
!
!        enddo
!      enddo



!     do ind=1,nat3d
!        if (ndc2(ind).eq.0) write (599,*) 'senza vicini', ind_tot, ind, dc2
!        write (599,*) 'senza vicini', ind_tot, ind, dc2
!        call flush(599)
!     enddo
  
END SUBROUTINE MOVE_BROWNIAN_SURF
