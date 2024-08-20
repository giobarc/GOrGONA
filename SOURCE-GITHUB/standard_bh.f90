Subroutine standard_bh(sbh_ener,sbh_enermin,sbh_ncall_output,sbh_label_outcoord,sbh_iwalk,sbh_nrestart,&
                       & steps_sincelast_gm,sbh_temperatura)
  
  USE PARAMETERS, ONLY: maxpar,cbol
  USE BH_MOD, ONLY: ind_tot,nat3d,xat,yat,zat,nglo,restart,xpos,ypos,zpos,xin,yin,zin,ener_restart
  USE FUNCTION_PAR, ONLY: arete,boltz_fact,ener

  
  Implicit none
  
  Integer :: i,j,sbh_ncall_output,sbh_iwalk,sbh_nrestart,prova
  Integer :: sbh_output_instruction
  Integer, parameter :: sbh_nmc=200 !the number of mc steps of the standard basin-hopping
  Integer :: steps_sincelast_gm
  Real(8) :: sbh_x(maxpar),sbh_y(maxpar),sbh_z(maxpar)
  Real(8) :: sbh_ener,sbh_enerold,sbh_deltaen
  Real(8) :: sbh_enermin
  Real(8) :: sbh_temperatura
  Real(8) :: rrr
  Logical :: sbh_label_outcoord
! DEBUG:
  integer :: n_init  


  sbh_output_instruction=0


  do j=1,nat3d
     sbh_x(j)=xpos(j,sbh_iwalk)
     sbh_y(j)=ypos(j,sbh_iwalk)
     sbh_z(j)=zpos(j,sbh_iwalk)
  enddo
  
  i=0
  do while (i .lt. sbh_nmc)
  i=i+1
     sbh_enerold=sbh_ener
          
     ind_tot=ind_tot+1
     
     do j=1,nat3d
        xin(j)=sbh_x(j)
        yin(j)=sbh_y(j)
        zin(j)=sbh_z(j)
     enddo
     
     call move_shell(ind_tot)
     
   !write(*,*)'2'  
     do j=1,nat3d
        xat(j)=xin(j)/arete(1)
        yat(j)=yin(j)/arete(1)
        zat(j)=zin(j)/arete(1)
     enddo

     prova=nat3d*3
     call minimization_rgl(prova,n_init)
     do j=1,nat3d
        xin(j)=xat(j)*arete(1)
        yin(j)=yat(j)*arete(1)
        zin(j)=zat(j)*arete(1)
     enddo
     sbh_ener=ener
     sbh_deltaen=sbh_ener-sbh_enerold
          
     !deciding whether accepting or not the move:
     if((-sbh_deltaen/(cbol*sbh_temperatura)) .gt. 0.)then
        boltz_fact=1.
     elseif((-sbh_deltaen/(cbol*sbh_temperatura)) .lt. 0)then	        
        boltz_fact=exp(-sbh_deltaen/(cbol*sbh_temperatura))
     elseif((-sbh_deltaen/(cbol*sbh_temperatura)) .eq. 0)then	
        write(*,*)'esponente di boltzmann nullo'
        boltz_fact=1.
     endif
     !write(*,*)'6'
     call ranmar(rrr)
     !write(*,*)'7'
     if(rrr .le. boltz_fact)then
	!the move is accepted
	!coordinates are updated:
        do j=1,nat3d
           sbh_x(j)=xin(j)
           sbh_y(j)=yin(j)
           sbh_z(j)=zin(j)
        enddo
     else
     sbh_ener=sbh_enerold
     endif
     !write(*,*)'8'
     !if the current lowest minimum has been found:
     if (sbh_ener.lt.sbh_enermin) then
        nglo=nglo+1
        sbh_ncall_output=sbh_ncall_output+1
        sbh_label_outcoord=.true.
        sbh_output_instruction=2
        call output_glo(sbh_iwalk,sbh_output_instruction)
        sbh_enermin=sbh_ener
	i=1 !standard_bh has to cycle for sbh_nmc more steps
     endif !su ener.lt enermin
     !write(*,*)'9'
     !writing output at every MC step:
     sbh_ncall_output=sbh_ncall_output+1
     sbh_label_outcoord=.false.
     sbh_output_instruction=3
     call output_glo(sbh_iwalk,sbh_output_instruction)
     !write(*,*)'10'
     !OPTION RESTART
     if(restart .eq. 'si')then
        if(sbh_ener .lt. ener_restart)then
           sbh_nrestart=sbh_nrestart+1
           write(*,*)'Step=',ind_tot
           write(*,*)'Ener=',sbh_ener
           write(*,*)'ener_restart=',ener_restart
           write(*,*)'Il programma riparte per la',sbh_nrestart,'-esima volta'	 
           !scrittura nel file successi.out
           if(sbh_nrestart .eq. 1)then
	      open(300,file='successi.out',status='unknown')
	      write(300,*)'#MCstep, ener, ripartenza num.'
	      call flush(300)
	      close(300)
           endif
           open(300,file='successi.out',status='unknown',position='append')
           write(300,'(i10,f16.8,i5)')ind_tot,sbh_ener,sbh_nrestart
           call flush(300)
           close(300)
           !optimization starts again from nwalkers random configurations:
           call initialization   
	   steps_sincelast_gm=0
	   return	     
        endif
     endif
     
  enddo
!write(*,*)'11'
  do j=1,nat3d
     xpos(j,sbh_iwalk)=sbh_x(j)
     ypos(j,sbh_iwalk)=sbh_y(j)
     zpos(j,sbh_iwalk)=sbh_z(j)
  enddo
!write(*,*)'12'
  
End Subroutine standard_bh
