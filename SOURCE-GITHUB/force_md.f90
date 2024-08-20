subroutine force_rgl
use iso_c_binding!in realta` e` uno statement per fortran 2003 ma se compila non mi lamento
USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster
USE STRUCTURES

Implicit None
Integer :: i,j,k,icolor,ierr
Integer :: itypik
Real*8 :: den(natom),frx(natom),fry(natom),frz(natom)
Real*8 :: dikm,dikm2,dikm3,dikm4,dikm5
Real*8 :: dik0,dik,xik,yik,zik,xuar,yvar,zwar
Real*8 :: ebi,eneri,eri,f,for,forsudik,denik
Real(8) :: fb(natom,natom)!,ds(natom,natom)
Real*8 :: espo,pexpp,qexpp,qsiexpq
Real*8 :: psudik0,dueapsudik0,dueq,qsi2,qsi2qsudik0


!write(*,*)'dentro a force i parametri valgono:'
!write(*,*)'p,        q,         a,        qsi'
!write(*,'(4f10.6)')p(1),q(1),a(1),qsi(1)
!write(*,'(4f10.6)')p(2),q(2),a(2),qsi(2)
!write(*,'(4f10.6)')p(3),q(3),a(3),qsi(3)

ener=0.d0                                                        



do i=1,natom
   dfx(i)=fx(i)
   dfy(i)=fy(i)
   dfz(i)=fz(i)
   den(i)=0.d0                                                       
   frx(i)=0.d0                                                       
   fry(i)=0.d0                                                       
   frz(i)=0.d0                                                       
enddo

do i=1,natom
   ebi=0.d0                                                       
   eri=0.d0                                                          
   eneri=0.d0
   !se l-atomo i non ha vicini, pongo la forza = 0
   do j=1,nvois(i) 
      k=ivois(j,i)  
!      write (*,*) "prima di k>i, i,k=",i,k
      if (k.gt.i) then
      if((itype(i).eq.1).and.(itype(k).eq.1)) then          
         itypik=1   ! stesso metallo A          
      else if((itype(i).eq.2).and.(itype(k).eq.2)) then             
         itypik=2   ! stesso metallo B          
      else             
         itypik=3  ! interazione A-B                
      endif      
      dik0=nn(itypik)
      psudik0=p(itypik)/dik0
      dueapsudik0=2.d0*a(itypik)*psudik0
      dueq=2.d0*q(itypik)
      qsi2=qsi(itypik)*qsi(itypik)
      qsi2qsudik0=qsi2*q(itypik)/dik0     
      xik=x(k)+u(k)-x(i)-u(i)
      yik=y(k)+v(k)-y(i)-v(i)
      zik=z(k)+w(k)-z(i)-w(i) 
!      xik=distanzex(i,k)                                     
!      yik=distanzey(i,k)                                           
!      zik=distanzez(i,k) 
!      ds(i,k)=dsqrt(xik*xik+yik*yik+zik*zik)           
      dik=dsqrt(xik*xik+yik*yik+zik*zik)
!      write (44,'(2i6,f14.9)') i,k,dik
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
!      write (44,'(2i6,2e18.8)') i,k,dik0,for

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
   
   ebi=dsqrt(den(i))                                                
   den(i)=1.d0/ebi                                                  
   eneri=eri-ebi
   ener=ener+eneri                                                   
!        write (33,'(i6,4e17.8)') i,ebi,eri,eneri,ener
!write (49,'(i6,4e17.8)') i, frx(i), fry(i), frz(i),den(i)

enddo !su i  
  do 30 i=1,natom
     IF(nvois(i) .gt. 0) THEN !exist neighbours
        do 40 j=1,nvois(i)
           k=ivois(j,i)
           if(k.gt.i)then
!              write (*,*)i,k,ds(i,k)
              xik=x(k)+u(k)-x(i)-u(i)
              yik=y(k)+v(k)-y(i)-v(i)
              zik=z(k)+w(k)-z(i)-w(i) 
!              xik=distanzex(i,k)                                     
!              yik=distanzey(i,k)                                           
!              zik=distanzez(i,k) 
!      ds(i,k)=dsqrt(xik*xik+yik*yik+zik*zik)           
      dik=dsqrt(xik*xik+yik*yik+zik*zik)
              denik=fb(k,i)*(den(i)+den(k))/dik
              frx(i)=frx(i)+denik*xik
              fry(i)=fry(i)+denik*yik
              frz(i)=frz(i)+denik*zik
              frx(k)=frx(k)-denik*xik
              fry(k)=fry(k)-denik*yik
              frz(k)=frz(k)-denik*zik
           endif!k>i
40      enddo
     ENDIF ! no neighbours
     fx(i)=frx(i)
     fy(i)=fry(i)
     fz(i)=frz(i)
!write (55,'(i6,3e18.10)') i, 2.89d0*dsqrt(2.d0)*fx(i),& 
! 2.89d0*dsqrt(2.d0)*fy(i), 2.89d0*dsqrt(2.d0)*fz(i)
30 enddo
!stop

EndSubroutine force_rgl
