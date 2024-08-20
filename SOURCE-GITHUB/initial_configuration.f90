SUBROUTINE INITIAL_CONFIGURATION

USE PARAMETERS
USE BH_MOD

implicit none
!local variables
Integer :: ind
Real(8) :: rho,theta,phi

do ind=1,nspec1
   spec(ind)=nomemet(1)
enddo

do ind=nspec1+1,nat3d
   spec(ind)=nomemet(2)
enddo

   do ind=1,nat3d
      call ranmar(zzz)
      xin(ind)=x_min+(x_max-x_min)*zzz   
      call ranmar(zzz)
      yin(ind)=y_min+(y_max-y_min)*zzz   
      call ranmar(zzz)
      zin(ind)=z_min+(z_max-z_min)*zzz   
   enddo

   open(22,file='initial_random.xyz',status='unknown')
   write (22,*) nat3d
   write (22,*) nomemet(1),' ',nomemet(2)

   do ind=1,nat3d
      write (22,'(a2,3f15.9)') spec(ind),xin(ind),yin(ind),zin(ind)
!      write (*,'(a2,3f15.9)') spec(ind),xin(ind),yin(ind),zin(ind)
   enddo
   close(22)
!write(*,*)'esco da initial config'

END SUBROUTINE INITIAL_CONFIGURATION
