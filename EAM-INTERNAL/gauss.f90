subroutine gauss
USE PARAMETERS
USE BH_MOD

Implicit None

Real*8 :: gset1,gset2,gset3,gset4,r                                          
Real*8 :: g4
Integer :: rip

rip=0

do 
   call ranmar(zzz)
   gset1=2.d0*zzz-1.d0
   call ranmar(zzz)
   gset2=2.d0*zzz-1.d0
   r=gset1*gset1+gset2*gset2 
   if (r<1) then     
      g1=gset1*sqrt(-2.d0*log(r)/r)
      g2=gset2*sqrt(-2.d0*log(r)/r)        
      rip=1        
   else        
      rip=0         
   endif
   if (rip==1) exit        
enddo

do          
   call ranmar(zzz)
   gset3=2.d0*zzz-1.d0
   call ranmar(zzz)
   gset4=2.d0*zzz-1.d0
   r=gset3*gset3+gset4*gset4
   if (r<1) then           
      g3=gset3*sqrt(-2.d0*log(r)/r)
      g4=gset4*sqrt(-2.d0*log(r)/r)       
      rip=2        
   else        
      rip=1         
   endif
   if (rip==2) exit       
enddo

End Subroutine GAUSS
