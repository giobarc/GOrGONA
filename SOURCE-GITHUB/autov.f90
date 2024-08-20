subroutine autopres(s00,s01,s02,s11,s12,s22,auto1,auto2,auto3)
use function_par
use bh_mod
use histo
implicit none

Real*8 :: s00,s01,s02,s11,s12,s22
Real*8 :: auto1,auto2,auto3
Real*8 :: bb,cc,dd,pp,qq
Real*8 :: alpha,theta1,theta2,theta3,lambda1,lambda2,lambda3,qqsualpha3!,pigreco

!pigreco=dacos(-1.d0)
!write (*,*) pigreco


  bb = -s00-s11-s22
  cc = s00*s11+s00*s22+s11*s22-s01*s01-s02*s02-s12*s12
  dd = s00*s12*s12+s11*s02*s02+s22*s01*s01-s00*s11*s22-2.d0*s02*s01*s12
  pp = cc - bb*bb/3.d0
  qq = 2.d0*bb*bb*bb/27.d0 - bb*cc/3.d0 +dd

  alpha = dsqrt(-pp/3.d0)
  qqsualpha3=-qq/(2.d0*alpha*alpha*alpha)
  if (qqsualpha3.gt.1.d0) qqsualpha3=1.d0
  if (qqsualpha3.lt.1.d0) qqsualpha3=-1.d0

  theta1 = dacos(qqsualpha3)/3.d0
  theta2 = dacos(qqsualpha3)/3.d0 + 2.d0*pigreco/3.d0
  theta3 = dacos(qqsualpha3)/3.d0 + 4.d0*pigreco/3.d0


  lambda1 = 2.d0*alpha*dcos(theta1) - bb/3.d0
  lambda2 = 2.d0*alpha*dcos(theta2) - bb/3.d0
  lambda3 = 2.d0*alpha*dcos(theta3) - bb/3.d0
!  write (117,*) alpha,qqsualpha3,lambda1,lambda2,lambda3
    
! delirio di istruzioni condizionate per ordinare tre numeri, così tanto per divertirsi
! auto1,auto2,auto3 in increasing order  
if (lambda1.lt.lambda2) then
   if (lambda1.lt.lambda3) then
       auto1=lambda1
       if (lambda2.lt.lambda3) then
           auto2=lambda2
           auto3=lambda3
       else
           auto2=lambda3
           auto3=lambda2
       endif
   else
       auto2=lambda1
       if (lambda2.lt.lambda3) then
           auto1=lambda2
           auto3=lambda3
       else
           auto1=lambda3
           auto3=lambda2
       endif
   endif
else 
   if (lambda1.lt.lambda3) then
       auto2=lambda1
       if (lambda2.lt.lambda3) then
           auto1=lambda2
           auto3=lambda3
       else
           auto1=lambda3
           auto3=lambda2
       endif
    else
       auto3=lambda1
       if (lambda2.lt.lambda3) then
           auto1=lambda2
           auto2=lambda3
       else
           auto1=lambda3
           auto2=lambda2
       endif
    endif
endif




end subroutine autopres
