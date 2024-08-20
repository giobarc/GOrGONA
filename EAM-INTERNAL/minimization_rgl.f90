
! DRIVER PER L'USO DI 'setulb'; ROUTINE DI MINIMIZZAZIONE DA MCS ANL. GOV.
! pacchetto di nome 'Lbfgsb.2.1.tar.gz' scaricabile da
! http://www.ece.nwu.edu/~nocedal/lbfgsb.html
!=========================================================================

Subroutine minimization_rgl(number_of_parameters,n_iter)

USE PARAMETERS
USE BH_MOD
USE FUNCTION_PAR

Implicit None

Integer, Intent(In) :: number_of_parameters
! Needed by the minimization routine setulb()
Character(60) :: task, csave
Logical :: lsave(4)
Integer :: n, m, iprint, nbd(number_of_parameters), &
           iwa(3*number_of_parameters)
Integer :: isave(44)
Real(8) :: f, factr, pgtol
Real(8) :: dsave(29)
Real(8), Allocatable :: wa(:)
Real(8) :: l(number_of_parameters)
Real(8) :: u(number_of_parameters)
Real(8) :: g(number_of_parameters)
Real(8) :: xdble(number_of_parameters)
!needed by the bigvoi routine
Logical :: first !Output for bigvoi
                 !If first=.true., bigvoi updates the number of neighbours
                 !for the metal-oxide interaction
! local 
Integer, Parameter :: itmax=2000 !The maximum number of calls to setulb
Real(8), Parameter :: epsi_lbfgsb=1.e-50_8 !Close to zero, used to check that 
                                           !setulb does not return the same energy
                                           !as in the previous step
Integer :: ind !counter
Integer :: n_z(nat3d) !It stores the number of neighbours during minimization
Integer :: MAX_ITER = 50 !The maximum number of neighbours updates
Integer :: n_iter !The number of neighbours updates
Real(8) :: gradiente(number_of_parameters)
Real(8) :: rtol,funzval_old,df,funzval
Integer :: iter !It counts the number of calls to Setulb()

!first = .true.
iprint = -1 ! We suppress the default output.
! We suppress both code-supplied stopping tests because the
! user is providing his own stopping criteria.
factr=1.0d0 !1.0d+7
pgtol=0.d0 !1.0d-5

! We now specify nbd which defines the bounds on the variables:
!                l   specifies the lower bounds,
!                u   specifies the upper bounds.
!     nbd(i)=0 if x(i) is unbounded,
!            1 if x(i) has only a lower bound,
!            2 if x(i) has both lower and upper bounds, and
!            3 if x(i) has only an upper bound.

!l(1:number_of_parameters)=0.000000000000001d0
!u(1:number_of_parameters)=0.999999999999999d0
nbd(1:number_of_parameters)=0 !unbounded minimization

! settaggio delle dimensioni del problema e delle 'limited memory corrections'
! N.B. valori raccomandati dagli autori per m: 3 <= m <= 20
n=number_of_parameters
m=7 ! m is the number of limited memory corrections

! allocations/initializations
Allocate(wa(2*m*n+4*n+12*m*m+12*m))
iter=0
iwa=0
wa=0.d0
isave=0
dsave=0.d0
lsave=.False.

! associazione del puntatore con le coordinate del punto in ingresso
do ind=1,nat3d
   xdble(3*ind-2)=xat(ind)
   xdble(3*ind-1)=yat(ind)
   xdble(3*ind)=zat(ind)
enddo

!serve per la routine di minimizzazione che vuole argomenti dbl-prc
funzval=ener
f=funzval

do ind=1,nat3d
   gradiente(3*ind-2)=-fx(ind)
   gradiente(3*ind-1)=-fy(ind)
   gradiente(3*ind)=-fz(ind)
enddo

g(1:number_of_parameters)=gradiente(1:number_of_parameters)
funzval_old=funzval

!We start the iteration by initializing task.
task = 'START'
!------- the beginning of the loop ----------
call bigvoi 

111  continue

iter=iter+1

call setulb(n,m,xdble,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,&
           csave,lsave,isave,dsave)

do ind=1,nat3d
!   if(xdble(3*ind)<=0*arete(1) .and. substrate=='si' )then
!      !here we check that no atoms are under the surface of the oxide support
!      is_nan = .true.
!      Deallocate(wa)
!      Return
!   endif
   xat(ind) = xdble(3*ind-2)
   yat(ind) = xdble(3*ind-1)
   zat(ind) = xdble(3*ind)
enddo

call bigvoi 
!if(first)then
!   do ind=1,nat3d
!      n_z(ind)=nv1(ind)
!   enddo
!   first=.false.
!endif

call force_rgl

!ATTENZIONE: in force_rgl is_nan e' posto .false.
!come puo' realizzarsi la condizione sotto?
!if(is_nan) then 
!  Deallocate(wa)
!  Return
!endif



IF (task(1:2) .eq. 'FG') THEN
         funzval=ener
         f=funzval
         do ind=1,nat3d
            gradiente(3*ind-2)=-fx(ind)
            gradiente(3*ind-1)=-fy(ind)
            gradiente(3*ind)=-fz(ind)
         enddo

         g(1:number_of_parameters)=gradiente(1:number_of_parameters)
         !   se il numero massimo di iterazioni e' raggiunto il migliore individuo
         !   trovato sino ad ora viene restituito in P0 e funzval
    	 If(iter.Ge.itmax) Then
!       	    Deallocate(wa)
            open (22,file='errors_bh.out',status='unknown',position='append')
            write(22,*)'WARNING:'
            !write(22,*)'At step ',ind_mc
            write(22,*)'the minimization exits with iter > itmax '
            close(22)
       	    Return
     	EndIf ! chiude If(rtol.Lt.ftol) Then
        !        go back to the minimization routine.
	 Goto 111
ENDIF ! chiude if (task(1:2) .eq. 'FG') then

IF (task(1:5) .eq. 'NEW_X') THEN
    ! the minimization routine has returned with a new iterate.
    ! At this point have the opportunity of stopping the iteration 
    ! or observing the values of certain parameters
    rtol=2._8*abs(funzval_old-funzval)/(abs(funzval_old)+abs(funzval))
    df=abs(funzval_old-funzval)
    funzval_old=funzval
    !   se la tolleranza e' soddisfatta o il numero massimo
    !   di iterazioni e' raggiunto il migliore individuo
    !   trovato sino ad ora viene restituito in P0 e funzval
    isucc=10
    If((iter.Ge.itmax .or.  n_iter > MAX_ITER ).and.(df.gt.epsi_lbfgsb)) Then
!       Deallocate(wa)
       open (22,file='errors_bh.out',status='unknown',position='append')
       write(22,*)'WARNING at step: ',ind_tot, ' the minimization exits with iter > itmax '       
       close(22)
       isucc=11
       !write (92,*) iter
       Return
    EndIf
	
    If((rtol.Lt.ftoll).and.(df.gt.epsi_lbfgsb)) Then

!	if( .not.first )then
!		call bigvoi

!		do ind=1,nat3d
!			if( (n_z(ind) .ne. nv1(ind)))then
!				first=.true.
!			endif
!			n_z(ind)=nv1(ind)
!		enddo
	 	call force_rgl
		if(is_nan) then
!			Deallocate(wa)
			Return
		endif
         	funzval=ener
         	f=funzval
         	do ind=1,nat3d
            		gradiente(3*ind-2)=-fx(ind)
            		gradiente(3*ind-1)=-fy(ind)
            		gradiente(3*ind)=-fz(ind)
         	enddo
    		g(1:number_of_parameters)=gradiente(1:number_of_parameters)

!		if(first)then
!			first=.false.
!			n_iter=n_iter+1
!			write(*,*)n_iter
!			goto 111 
!		endif
!	endif
       !write (92,*) iter

!	Deallocate(wa)
       	Return
     EndIf ! chiude If(rtol.Lt.ftoll) Then
ENDIF ! chiude if (task(1:5) .eq. 'NEW_X') then   

	Goto 111

End Subroutine minimization_rgl


