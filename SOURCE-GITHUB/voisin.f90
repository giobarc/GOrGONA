subroutine voisin

USE PARAMETERS  !uso il modulo di definizione dei parametri
USE CLUSTER     !uso il modulo dove defin1:natomco variabili e parametri cluster
USE STRUCTURES

Implicit None
Integer :: i,jv,j,ikjk
Real*8 :: dij2,xij,yij,zij
Real*8 :: cutoff_end2(3)
 do i=1,3
    cutoff_end2(i)=cutoff_end(i)*cutoff_end(i)! !voisin.f90 works with Angstrom coordinates
 enddo
do i=1,natom
   nvois(i)=0
enddo

do i=1,natom-1                                                   
   do jv=1,nv4(i)                                                   
      j=iv4(jv,i)   
                      
      ikjk=3
      if((itype(i).eq.1).and.(itype(j).eq.1)) then
        ikjk=1
      else if ((itype(i).eq.2).and.(itype(j).eq.2)) then
        ikjk=2
      endif

!      distanzex(i,j)=x(j)+u(j)-x(i)-u(i)                     
!      distanzey(i,j)=y(j)+v(j)-y(i)-v(i)        
!      distanzez(i,j)=z(j)+w(j)-z(i)-w(i)
      xij=x(j)+u(j)-x(i)-u(i)
      yij=y(j)+v(j)-y(i)-v(i)
      zij=z(j)+w(j)-z(i)-w(i) 

        !periodic boundary conditions
!        if(boundary_cond .eq. 'si')then
!           if (abs(xij) .gt. x_boxedge/2.) then
!              if (xij .gt. 0.)then
!                 xij=xij-x_boxedge
!              else
!                 xij=xij+x_boxedge
!              endif
!           endif
!           if (abs(yij) .gt. y_boxedge/2.) then
!              if (yij .gt. 0.)then
!                 yij=yij-y_boxedge
!              else
!                 yij=yij+y_boxedge
!              endif
!           endif
!           if(zboundary_cond .eq. 'si')then
!             if (abs(zij) .gt. z_boxedge/2.) then
!                if (zij .gt. 0.)then
!                   zij=zij-z_boxedge
!                else
!                   zij=zij+z_boxedge
!                endif
!             endif
!           endif

!      distanzex(i,j)=xij
!      distanzey(i,j)=yij
!      distanzez(i,j)=zij

!        endif
      
      dij2=xij*xij+yij*yij+zij*zij

       
      !neighbor list update for force_rgl.f90
      if (dij2.lt.cutoff_end2(ikjk)) then
         nvois(i)=nvois(i)+1
         nvois(j)=nvois(j)+1
         ivois(nvois(i),i)=j
         ivois(nvois(j),j)=i

      endif
      
enddo      !fine ciclo su jv
   
!   if (nvois(i).ge.nmax) then                                         
!      write (*,'(i4,1x,i4,a18)'),i,nvois(i),'troppi NN in voisin'
!      stop                                                          
!   endif
enddo       !fine ciclo su i

!check
!write(*,*)'second check loop in voisin:'
!do i=10,15
!   write(*,*)i,' ha ',nvois(i),' vicini'
!   write(*,*)'lista:'
!   do j=1,nvois(i)
!      write(*,*)ivois(j,i)
!  enddo
!enddo



!simmetrizzo la matrice delle distanze 
!(solo le componenti aggiornate da voisin):
!do i=1,natom-1                                                   
!   do jv=1,nv4(i)                                                   
!      j=iv4(jv,i) 
!      distanzex(j,i)=-distanzex(i,j)
!      distanzey(j,i)=-distanzey(i,j)
!      distanzez(j,i)=-distanzez(i,j)
!   enddo
!enddo

EndSubroutine voisin
