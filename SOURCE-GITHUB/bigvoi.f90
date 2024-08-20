subroutine bigvoi

USE FUNCTION_PAR, Only: nv1,nmax0,nvmax0,nvois,ivois,&
                        nat3d,arete,itype,cutz,cutoff_end
USE BH_MOD, Only: x => xat,y => yat, z => zat

  	implicit none
!	logical, intent(in) ::first
	integer*4 ::i,j,ikjk
  	Real(8) :: dij2,xij,yij,zij
        Real(8) :: cutoff_end2(3)

      do i=1,3  
          cutoff_end2(i)=cutoff_end(i)*cutoff_end(i)
      enddo
      do i=1,nat3d
        nvois(i)=0
 !       if(first) nv1(i)=0
      enddo
      


      do i=1,nat3d-1
         do j=i+1,nat3d
           ikjk=3
           if((itype(i).eq.1).and.(itype(j).eq.1)) then
             ikjk=1
           else if ((itype(i).eq.2).and.(itype(j).eq.2)) then
             ikjk=2
           endif
           xij=x(j)-x(i)
           yij=y(j)-y(i)
           zij=z(j)-z(i)
           dij2=xij*xij+yij*yij+zij*zij

!	   if(first .and. (dij2 .lt. cutz(itype(i),itype(j))))then
!      			nv1(i)=nv1(i)+1
!      			nv1(j)=nv1(j)+1
!	   endif
           if (dij2.lt.cutoff_end2(ikjk)) then
              	nvois(i) = nvois(i)+1
              	ivois(nvois(i),i) = j
              	nvois(j) = nvois(j)+1
              	ivois(nvois(j),j) = i
           endif
        enddo

        if (nvois(i).ge.nvmax0) then
           write (*,530) i,nvois(i)
           stop
        endif
	

      enddo
! stop
  530 format (' troppi vicini in bigvoi  ',2i5)

      return
      end
