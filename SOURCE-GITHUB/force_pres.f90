SUBROUTINE force_pres
USE FUNCTION_PAR
USE BH_MOD
USE HISTO
implicit none
! Local variables

Real(8) :: den(nmax0),frx(nmax0),fry(nmax0),frz(nmax0)
Real(8) :: x(nmax0),y(nmax0),z(nmax0)
Real(8) :: ebi,eri,eneri,for,forsudik,denik,presri,presbi,prsr,prsb
Real(8) :: fbx,fby,fbz,f
Real(8) :: xik,yik,zik,dik
Real(8) :: s00,s01,s02,s11,s12,s22
Real(8) :: s00b,s01b,s02b,s11b,s12b,s22b
Real(8) :: s00r,s01r,s02r,s11r,s12r,s22r
Real(8) :: auto1,auto2,auto3
Real(8) :: dik0,psudik0,dueapsudik0,dueq,qsi2,qsi2qsudik0,pexp,qexp,espo
Real(8) :: dikm,dikm2,dikm3,dikm4,dikm5
Integer :: i,j,k,itypik

!f(dabs(zat(1))>3275786477704)write(*,*)'0000000000 zat(1)=',zat(1)
!write(*,*)'force_rgl'
!write(*,*)'cutoff_start ',cutoff_start
!write(*,*)'cutoff_end ',cutoff_end


!write(*,*)a3(1),a3(2),a3(3)
!write(*,*)a4(1),a4(2),a4(3)
!write(*,*)a5(1),a5(2),a5(3)
!write(*,*)x3(1),x3(2),x3(3)
!write(*,*)x4(1),x4(2),x4(3)
!write(*,*)x5(1),x5(2),x5(3)

!     ener_debug = ener
      ener=0.d0

      do i=1,nat3d
       dfx(i)=fx(i)
       dfy(i)=fy(i)
       dfz(i)=fz(i)
       x(i)=xat(i)
       y(i)=yat(i)
       z(i)=zat(i)
      enddo

       do 10 i=1,nat3d
        pres_atom(i)=0.d0
        den(i)=0.d0
        frx(i)=0.d0
        fry(i)=0.d0
        frz(i)=0.d0
        ebi=0.d0
        eri=0.d0
        eneri=0.d0
        presri=0.d0
        presbi=0.d0
        s00=0.d0
        s01=0.d0
        s02=0.d0
        s11=0.d0
        s12=0.d0
        s22=0.d0
        s00b=0.d0
        s01b=0.d0
        s02b=0.d0
        s11b=0.d0
        s12b=0.d0
        s22b=0.d0
        s00r=0.d0
        s01r=0.d0
        s02r=0.d0
        s11r=0.d0
        s12r=0.d0
        s22r=0.d0

        IF(nvois(i).eq.0) then
        	frx(i)=0.d0
       		fry(i)=0.d0
        	frz(i)=0.d0
	ELSE !there exist neighbours
	  	do 20 j=1,nvois(i)
           		k=ivois(j,i)

		        if((itype(i).eq.1).and.(itype(k).eq.1)) then
            			 itypik=1   ! stesso metallo A
		        else if((itype(i).eq.2).and.(itype(k).eq.2)) then
             			 itypik=2   ! stesso metallo B
          		else
             			itypik=3  ! interazione A-B
          		endif

           !!  itypik=1

			dik0=dist(itypik)
 		        psudik0=p(itypik)/dik0
          		dueapsudik0=2.d0*a(itypik)*psudik0
          		dueq=2.d0*q(itypik)
          		qsi2=qsi(itypik)*qsi(itypik)
          		qsi2qsudik0=qsi2*q(itypik)/dik0

          		xik=x(k)-x(i)
          		yik=y(k)-y(i)
          		zik=z(k)-z(i)
          		dispx(j,i)=xik
          		dispy(j,i)=yik
          		dispz(j,i)=zik
          		dik=dsqrt(xik*xik+yik*yik+zik*zik)
          		if (dik.lt.cutoff_start(itypik)) then
             			espo=dik/dik0-1.d0
             			pexp=dexp(-p(itypik)*espo)
             			qexp=dexp(-dueq*espo)
             			den(i)=qsi2*qexp+den(i)
             			for=pexp*dueapsudik0
             			eri=eri+a(itypik)*pexp
             			prsr=p(itypik)*a(itypik)*pexp/dik0
             			prsb=-q(itypik)*qsi2*qexp/dik0
          		else
             			dikm=dik-cutoff_end(itypik)
             			dikm2=dikm*dikm
             			dikm3=dikm2*dikm
             			dikm4=dikm3*dikm
             			dikm5=dikm4*dikm
             			den(i)=den(i)+(x5(itypik)*dikm5+x4(itypik)*dikm4 &
                                           +x3(itypik)*dikm3)**2
             			for=-2.d0*(5.d0*a5(itypik)*dikm4+4.d0*a4(itypik)*dikm3 &
                                            +3.d0*a3(itypik)*dikm2)
             			eri=eri+a5(itypik)*dikm5+a4(itypik)*dikm4+a3(itypik)*dikm3
				prsr=(5.d0*a5(itypik)*dikm4+4.d0*a4(itypik)*dikm3+3.d0*a3(itypik)*dikm2)
 prsb=(x5(itypik)*dikm5+x4(itypik)*dikm4+x3(itypik)*dikm3)*(5.d0*x5(itypik)*dikm4+&
&4.d0*x4(itypik)*dikm3+3.d0*x3(itypik)*dikm2)
          		endif
			prsr=prsr*dik
			prsb=prsb*dik
			presri=presri+prsr
			presbi=presbi+prsb
			!stress tensor
			s00r=s00r+prsr*xik*xik/(dik*dik)
			s00b=s00b+prsb*xik*xik/(dik*dik)
			s01r=s01r+prsr*xik*yik/(dik*dik)
			s01b=s01b+prsb*xik*yik/(dik*dik)
			s02r=s02r+prsr*xik*zik/(dik*dik)
			s02b=s02b+prsb*xik*zik/(dik*dik)
			s11r=s11r+prsr*yik*yik/(dik*dik)
			s11b=s11b+prsb*yik*yik/(dik*dik)
			s12r=s12r+prsr*yik*zik/(dik*dik)
			s12b=s12b+prsb*yik*zik/(dik*dik)
			s22r=s22r+prsr*zik*zik/(dik*dik)
			s22b=s22b+prsb*zik*zik/(dik*dik)
		
          		forsudik=for/dik
          		frx(i)=frx(i)-forsudik*xik
          		fry(i)=fry(i)-forsudik*yik
          		frz(i)=frz(i)-forsudik*zik

   		20   continue
        	ebi=-dsqrt(den(i))
        	den(i)=-1.d0/ebi
        	eneri=ebi+eri
		ener_atom(i)=eneri
		presbi=den(i)*presbi
		s00b=s00b*den(i)
		s01b=s01b*den(i)
		s02b=s02b*den(i)
		s11b=s11b*den(i)
		s12b=s12b*den(i)
		s22b=s22b*den(i)
		s00=s00r+s00b
		s01=s01r+s01b
		s02=s02r+s02b
		s11=s11r+s11b
		s12=s12r+s12b
		s22=s22r+s22b
		call autopres(s00,s01,s02,s11,s12,s22,auto1,auto2,auto3)
		pres_atom(i)=presri+presbi
		stress_atom(i,1)=auto1
		stress_atom(i,2)=auto2
		stress_atom(i,3)=auto3
!		write (114,*) i, pres_atom(i),auto1,auto2,auto3
!              write (111,'(i6,3f18.8)') i, presri,presbi,pres_atom(i)
        	ener=ener+eneri
	ENDIF !exist neighbours	

   10 continue

   emetal=ener
   do 30 i=1,nat3d
        fbx=0.d0
        fby=0.d0
        fbz=0.d0

	IF(nvois(i) .ne. 0) THEN !exist neighbours
           do 40 j=1,nvois(i)
          	k=ivois(j,i)

          	if((itype(i).eq.1).and.(itype(k).eq.1)) then
             		itypik=1   ! A-A
          	else if((itype(i).eq.2).and.(itype(k).eq.2)) then
             		itypik=2   ! B-B
          	else
             		itypik=3   ! A-B
          	endif

          	dik0=dist(itypik)
          	psudik0=p(itypik)/dik0
          	dueapsudik0=2.d0*a(itypik)*psudik0
          	dueq=2.d0*q(itypik)
          	qsi2=qsi(itypik)*qsi(itypik)
          	qsi2qsudik0=qsi2*q(itypik)/dik0
!
          	xik=dispx(j,i)
          	yik=dispy(j,i)
          	zik=dispz(j,i)
          	dik=dsqrt(xik*xik+yik*yik+zik*zik)

          	if (dik.lt.cutoff_start(itypik)) then
             	espo=dik/dik0-1.d0
             	qexp=dexp(-dueq*espo)
             	f=qsi2qsudik0*qexp
          	else
             		dikm=dik-cutoff_end(itypik)
             		dikm2=dikm*dikm
             		dikm3=dikm2*dikm
             		dikm4=dikm3*dikm
             		dikm5=dikm4*dikm
             		f=-(x5(itypik)*dikm5+x4(itypik)*dikm4+x3(itypik)*dikm3)* &
              			(5.d0*x5(itypik)*dikm4+4.d0*x4(itypik)*dikm3 &
              			+3.d0*x3(itypik)*dikm2)
          	endif
          	denik=f*(den(i)+den(k))/dik
          	fbx=fbx+denik*xik
          	fby=fby+denik*yik
          	fbz=fbz+denik*zik
       	   40   continue
	ENDIF ! no neighbours
     	fx(i)=frx(i)+fbx  
        fy(i)=fry(i)+fby
        fz(i)=frz(i)+fbz
   30 continue

!       pressure and stress components are converted from eV/A^3 to GPa

     do i=1,nat3d
        pres_atom(i)=pres_atom(i)*4.d0/(3.d0*arete(1)**3)*1.60217662d2 
        stress_atom(i,1)=stress_atom(i,1)*4.d0/(3.d0*arete(1)**3)*1.60217662d2 
        stress_atom(i,2)=stress_atom(i,2)*4.d0/(3.d0*arete(1)**3)*1.60217662d2 
        stress_atom(i,3)=stress_atom(i,3)*4.d0/(3.d0*arete(1)**3)*1.60217662d2 
     enddo  

!if(substrate == 'si') then
!    call force_met_mgo
!else
!    e_met_mgo=0
    is_nan = .false.
!endif

!write(*,*)'force_rgl, ener= ',ener
!stop

END SUBROUTINE force_pres
