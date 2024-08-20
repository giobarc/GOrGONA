SUBROUTINE FORCE_RGL

USE PARAMETERS
USE FUNCTION_PAR
USE BH_MOD

implicit none

! Local variables
Real(8) :: den(nmax0),frx(nmax0),fry(nmax0),frz(nmax0)
Real(8) :: x(nmax0),y(nmax0),z(nmax0)
Real(8) :: fb(nat3d,nat3d)!,ds(nat3d,nat3d)
Real(8) :: ebi,eri,eneri,for,forsudik,denik
Real(8) :: fbx,fby,fbz,f
Real(8) :: xik,yik,zik,dik
Real(8) :: dik0,psudik0,dueapsudik0,dueq,qsi2,qsi2qsudik0,pexp,qexp,espo
Real(8) :: dikm,dikm2,dikm3,dikm4,dikm5
Real(8) :: pexpp,qexpp,qsiexpq
Integer :: i,j,k,itypik

 ener=0.d0

 do i=1,nat3d
    dfx(i)=fx(i)
    dfy(i)=fy(i)
    dfz(i)=fz(i)
    x(i)=xat(i)
    y(i)=yat(i)
    z(i)=zat(i)
    den(i)=0.d0                                                       
    frx(i)=0.d0                                                       
    fry(i)=0.d0                                                       
    frz(i)=0.d0                                                       
 enddo

do i=1,nat3d
!   write (78,*) i, nvois(i)
   ebi=0.d0                                                       
   eri=0.d0                                                          
   eneri=0.d0
   do j=1,nvois(i) 
      k=ivois(j,i)  
      if (k.gt.i) then
         if ((itype(i).eq.1).and.(itype(k).eq.1)) then          
            itypik=1   ! stesso metallo A          
         else if ((itype(i).eq.2).and.(itype(k).eq.2)) then             
            itypik=2   ! stesso metallo B          
         else             
            itypik=3  ! interazione A-B                
         endif      
         dik0=dist(itypik)
         psudik0=p(itypik)/dik0
         dueapsudik0=2.d0*a(itypik)*psudik0
         dueq=2.d0*q(itypik)
         qsi2=qsi(itypik)*qsi(itypik)
         qsi2qsudik0=qsi2*q(itypik)/dik0     
         xik=x(k)-x(i)
         yik=y(k)-y(i)
         zik=z(k)-z(i) 
         dik=dsqrt(xik*xik+yik*yik+zik*zik)
         if (dik.lt.cutoff_start(itypik)) then                         
            espo=1.d0-dik/dik0                                
            pexpp=dexp(p(itypik)*espo)
            qexpp=dexp(dueq*espo)
            qsiexpq=qsi2*qexpp
            for=pexpp*dueapsudik0
            eri=eri+2.d0*a(itypik)*pexpp
            fb(k,i) = (q(itypik)/dik0)*qsiexpq
         else 
            dikm=dik-cutoff_end(itypik)
            dikm2=dikm*dikm
            dikm3=dikm2*dikm
            dikm4=dikm3*dikm
            dikm5=dikm4*dikm  
            qsiexpq=x5(itypik)*dikm5+x4(itypik)*dikm4+x3(itypik)*dikm3
            for=-2.d0*(5.d0*a5(itypik)*dikm4+4.d0*a4(itypik)*dikm3+3.d0*a3(itypik)*dikm2)
            eri=eri+2.d0*(a5(itypik)*dikm5+a4(itypik)*dikm4+a3(itypik)*dikm3)       
            fb(k,i) = -((5.d0*x5(itypik)*dikm4+4.d0*x4(itypik)*dikm3+&
                    & 3.d0*x3(itypik)*dikm2))*qsiexpq
            qsiexpq = qsiexpq**2
         endif
         den(i)=qsiexpq+den(i)
         den(k)=qsiexpq+den(k)
         forsudik=for/dik
         frx(i)=frx(i)-forsudik*xik                                   
         fry(i)=fry(i)-forsudik*yik                                    
         frz(i)=frz(i)-forsudik*zik  
         frx(k)=frx(k)+forsudik*xik                                   
         fry(k)=fry(k)+forsudik*yik                                    
         frz(k)=frz(k)+forsudik*zik  
      endif !su k>i   
   enddo !su j    
   ebi=-dsqrt(den(i))                                                
   den(i)=-1.d0/ebi                                                  
   eneri=eri+ebi
   ener=ener+eneri
   energ(i)=eneri
   !write (*,*) i,energ(i)           
 enddo !su i

! write (*,*) ener
! stop
 do i=1,nat3d
    if (nvois(i).gt.0) then !exist neighbours
       do j=1,nvois(i)
          k=ivois(j,i)
          if (k.gt.i)then
             xik=x(k)-x(i)
             yik=y(k)-y(i)
             zik=z(k)-z(i) 
             dik=dsqrt(xik*xik+yik*yik+zik*zik)
             denik=fb(k,i)*(den(i)+den(k))/dik
             frx(i)=frx(i)+denik*xik
             fry(i)=fry(i)+denik*yik
             frz(i)=frz(i)+denik*zik
             frx(k)=frx(k)-denik*xik
             fry(k)=fry(k)-denik*yik
             frz(k)=frz(k)-denik*zik
          endif!k>i
       enddo! su j
     endif ! on nvois(i)
     fx(i)=frx(i)
     fy(i)=fry(i)
     fz(i)=frz(i)
    
 enddo !su i
 
 e_met_mgo=0
 is_nan = .false.

END SUBROUTINE FORCE_RGL
