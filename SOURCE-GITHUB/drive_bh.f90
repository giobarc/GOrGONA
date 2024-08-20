SUBROUTINE DRIVE_BH(n_err,n_iter_max)
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  USE HISTO
  
  implicit none

  !local variables
  integer :: ind,ind2,number_of_parameters,ii,indpres
  integer :: output_instruction
             ! = 1 first call to the output subroutine, before MC loop starts - file headings
             ! = 2 current lowest energy minimum has been found - its coordinates are saved
             ! = 3 energy and order parameter at the current MC step are recorded
             ! = 4 every 1000 steps the percentage of accepted/rejected moves are recorded
  integer :: ncall_output
  integer :: iwk,iwalk
  integer :: lt1neigh,lt3neigh
  integer :: n_neigh_old,n_neigh_old2
  real(8) :: enermin
  real(8) :: epsihisto
  logical :: label_output_coordinates,label_output_coordinates2,sbh,low_walker
  !modifica temperatura
  real(8) :: temperatura,sbh_temperatura,lowenergy,highenergy,eps
  real(8) :: pr1,pr2,pr3,pr4,pr5,pr6,pr7,pr8,pr9,pr10
  integer :: sbh_neighbours(15),lowenergy_walker,highenergy_walker,steps_sincelast_gm
  integer :: nrestart
  integer :: sbh_freq,sbh_start
  integer :: nlista,n_pass_all(2),n_iter(2)

  logical :: isnan, minimization
  
! refresh random num. generator:
  character*8 :: date
  character*10 :: time
  integer :: rnd_num1,rnd_num2

! DEBUG:
  Integer:: n_err, n_debug, n_iter_max

  write(*,*)'IN DRIVE CUTOFF',cutoff_start,cutoff_end   
  !These assignments regards the call to standard basin hopping
  sbh=.false. !if .false., standard basin-hopping is never performed
  sbh_start=100 !if sbh=.true., standard basin_hopping is performed starting from sbh_start steps
  sbh_freq=3000 !if no new minima have been found but the last sbh call was sbh_freq steps ago, sbh starts again
                !and it is applied to the lowest/highest energy walker
  sbh_temperatura=1500.
  low_walker=.false. !the lowest energy walker (true) or the highest energy walker (false) performs sbh
  n_err = 0
  histo_sig555(:)=0
  
  output_instruction=0


!  setting of random numbers generator
  if(refresh == 'no') then
        call rmarin(1234,5678)
  else
	call date_and_time(date,time)
	time(7:9)=time(8:10)
	rnd_num1=1000*(iachar(time(6:6))-iachar('0'))+100*(iachar(time(7:7))-iachar('0'))+10*&
		(iachar(time(8:8))-iachar('0'))+iachar(time(9:9))-iachar('0')
	if(rnd_num1<1000)rnd_num1=rnd_num1+1000
	rnd_num2=1000*(iachar(time(7:7))-iachar('0'))+100*(iachar(time(8:8))-iachar('0'))+10*&
		(iachar(time(9:9))-iachar('0'))+iachar(time(10:10))-iachar('0')
	if(rnd_num2<1000)rnd_num2=rnd_num2+1000
	write(*,*) rnd_num1,' ' ,rnd_num2
        call rmarin(rnd_num1,rnd_num2)
  endif
  

  !printing to terminal:
  write (*,*) 'imet(1),imet(2)=', imet(1),imet(2)
  write (*,*) 'p= ',(p(ind),ind=1,3)
  write (*,*) 'q= ',(q(ind),ind=1,3)
  write (*,*) 'a= ',(a(ind),ind=1,3)
  write (*,*) 'qsi= ',(qsi(ind),ind=1,3)
  !write (*,*) ' string_global=',string_global

  
  !conversion of coordinates, initialization and first minimization
  number_of_parameters=3*nat3d

  first_minimization = .true.
  !cluster(s) is(are) initialized:
  call initialization 
!  first_minimization = .false. !it is setted .false. aftre the first minima is found
  !first call to output_glo.f90
  ncall_output=1
  iwalk=1
  output_instruction=1
  call output_glo(iwalk,output_instruction) ! writing headers of output files
  nglo=1
  output_instruction=2
  call output_glo(iwalk,output_instruction) ! writing headers of output files

  ! counters initializations
   !counter for successive lowest energy found
  ind_tot=0 !counter for the total number of PC steps
  noutput=1 !number of output files
  enermin=0 !the lowest energy encountered
  nrestart=0 !the times ener_restart has been reached
  n_accepted_moves=0 !the number of accepted moves
  n_rejected_moves=0 !the number of rejected moves
  n_true_moves=0 !the number of accepted moves towards a different minimum
  nlista=0
  sbh_neighbours(:)=0
  steps_sincelast_gm=0
  !probmove1,probmove2,probmove3,ecc... have been assigned by input file
  !the probability array used to choose the kind of move is assigned here:
  pr1=probmove1 ! zzz < pr1 --> move bonds
  pr2=pr1+probmove2 ! pr1 < zzz < pr2 --> move ball
  pr3=pr2+probmove3 ! pr2 < zzz < pr3 --> move shell
  pr4=pr3+probmove4 ! pr3 < zzz < pr4 --> move shake
  pr5=pr4+probmove5 ! pr4 < zzz < pr5 --> move highenergyatoms
  pr6=pr5+probmove6 ! pr5 < zzz < pr6 --> move single
  pr7=pr6+probmove7 ! pr6 < zzz < pr7 --> move brownian
  pr8=pr7+probmove8 ! pr7 < zzz < pr8 --> move brownian surf
  pr9=pr8+probmove9 ! pr8 < zzz < pr9 --> move mirror
  pr10=pr9+probmove10 ! pr9 < zzz < pr10 --> move exchange
  if (pr10.lt.0.9999999999d0) then
      write(*,*) pr10,"problem in the sum of probabilities of moves"
      stop
  endif

  do iwk=1,nwalkers
     contwalker(iwk)=1
  enddo

  do ind_mc=1,n_mc
  !write(*,*)ind_mc
     	ind_tot=ind_tot+1 !counter of total number of MC steps
     	steps_sincelast_gm=steps_sincelast_gm+1     
     	if (ind_tot.ge.n_mc) return
     	!is it the right moment to perform some standard bh steps?
     	if((sbh) .and. (ind_mc .gt. sbh_start))then
!	  if((mod(ind_tot,sbh_freq) .eq. 0) .or. (label_output_coordinates2))then 
     	  if(mod(ind_tot,sbh_freq) .eq. 0)then 
		if(mod(ind_tot,sbh_freq) .eq. 0)then
        	   !looking for the walker with the lowest energy
        	   lowenergy=enerwalker(1)
        	   lowenergy_walker=1
        	   do ind=1,nwalkers
        	      if(enerwalker(ind) .lt. lowenergy)then
        	         lowenergy=enerwalker(ind)
        	         lowenergy_walker=ind
        	      endif
        	   enddo
		   highenergy=enerwalker(1)
        	   highenergy_walker=1
        	   do ind=1,nwalkers
        	      if(enerwalker(ind) .gt. highenergy)then
        	         highenergy=enerwalker(ind)
        	         highenergy_walker=ind
        	      endif
        	   enddo
        	   !...the program performes some steps of standard bh at sbh_temperature
		   if(low_walker)then
        	        call standard_bh(ener,enermin,ncall_output,label_output_coordinates, &
			  lowenergy_walker,nrestart,steps_sincelast_gm,sbh_temperatura)
		   elseif(.not.low_walker)then
		        write(*,*)'chiamo sbh, sul walker',highenergy_walker
		        call standard_bh(ener,enermin,ncall_output,label_output_coordinates,&
			  highenergy_walker,nrestart,steps_sincelast_gm,sbh_temperatura)
		   endif
		elseif(label_output_coordinates2)then
        	   call standard_bh(ener,enermin,ncall_output,label_output_coordinates,iwalk,nrestart,&
        	        & steps_sincelast_gm,sbh_temperatura)
        	   label_output_coordinates2=.false.
		endif
     	  endif
     	endif !sbh
     	!sbh_neighbours(iwalk)=n_neigh_old
     
     	!choosing the walker to be moved:
     	SELECT CASE (walkers_turnover)
     	CASE(1) ! walkers move as 1,2...,nwalkers,1,2...,nwalkers,1,2...
        	if(iwalk .ne. nwalkers)iwalk=iwalk+1
        	if(iwalk .eq. nwalkers)iwalk=1
     	CASE(2) ! walkers random exctraction
        	call ranmar(zzz)
        	iwalk=int(nwalkers*zzz)+1
        	iwalk=min(iwalk,nwalkers)
     	CASE(3)
	if(mod(ind_tot,freq_walker_exchange) .eq. 0)then
           if(iwalk .ne. nwalkers)then
              iwalk=iwalk+1
           else
              iwalk=1
           endif
	endif
	
     	END SELECT
          
     	call cont_neigh(iwalk,n_neigh_old) !checking whether the iwk walker has any neighbors

     	contwalker(iwalk)=contwalker(iwalk)+1 !MC steps counter for walker iwalk
 
     	ener_old=enerwalker(iwalk)
!        write (*,*) "ener_old=",ener_old
        esub_old=esubwalker(iwalk)
!        write (*,*) "esub_old=",esub_old
     	walker_parameter_old(:)=walker_parameter(iwalk,:)
     
     	do ind=1,nat3d
     		xin(ind)=xpos(ind,iwalk)
        	yin(ind)=ypos(ind,iwalk)
        	zin(ind)=zpos(ind,iwalk)
                spec(ind)=specpos(ind,iwalk)
                !write(*,*)xin(ind),yin(ind),zin(ind)
     	enddo


	!temperature

	    temperatura=0.d0


	minimization=.false.
	eps=1.0*arete(1) !the lowest value of z-coord. in xin
	do while (.not.minimization)	
        !zzz determina la scelta della mossa.
 	       call ranmar(zzz)
        	if(zzz.le.pr1) then
        	   string_move='bonds'
        	   call move_bonds(eps)
                   temperatura=temp_bonds
        	elseif((zzz.gt.pr1).and.(zzz.le.pr2)) then
        	   string_move='ball'
        	   call move_ball(eps)	!(ind_tot)
                   temperatura=temp_ball
		elseif((zzz.gt.pr2).and.(zzz.le.pr3)) then
        	   string_move='shell'
        	   call move_shell(eps)	!(ind_tot)
                   temperatura=temp_shell
        	elseif((zzz.gt.pr3).and.(zzz.le.pr4)) then
        	   string_move='shake'
        	   call move_shake(eps)
                   temperatura=temp_shake
		elseif((zzz.gt.pr4).and.(zzz.le.pr5)) then
        	   string_move='highe'
        	   call move_highenergyatoms(eps)
                   temperatura=temp_highener
        	elseif((zzz.gt.pr5).and.(zzz.le.pr6)) then
        	   string_move='single'
        	   call move_single
                   temperatura=temp_single
        	elseif((zzz.gt.pr6).and.(zzz.le.pr7)) then
        	   string_move='brow'
        	   call move_brownian
                   temperatura=temp_brow_accept
         	elseif((zzz.gt.pr7).and.(zzz.le.pr8)) then
        	   string_move='brsu'
        	   call move_brownian_surf
                   temperatura=temp_browsurf_accept
         	elseif((zzz.gt.pr8).and.(zzz.le.pr9)) then
        	   string_move='mirr'
        	   call move_mirror
                   temperatura=temp_mirror
          	elseif((zzz.gt.pr9).and.(zzz.le.pr10)) then
        	   string_move='exch'
        	   call move_exchange(temperatura) 
        	   ! here the variable temperatura is set inside the move_exchange routine; it depends on the exchange type
         	endif

		do ind=1,nat3d
        	   xsave(ind)=xat(ind)
        	   ysave(ind)=yat(ind)
        	   zsave(ind)=zat(ind)
        	   xat(ind)=xin(ind)/arete(1)
        	   yat(ind)=yin(ind)/arete(1)
		   zat(ind)=zin(ind)/arete(1)
        	enddo
     		call cont_neigh(iwalk,n_neigh_old2)
		n_iter(2)=0
		n_iter(1)=0

		call minimization_rgl(number_of_parameters,n_iter(2))
!                write (*,*) "after minimization, ind_mc=",ind_mc, e_met_mgo
 		if((ener .lt. 0) .and. (.not. isnan(ener)) .and. (.not. is_nan))then
			 minimization=.true.
			if(n_iter_max<n_iter(2))n_iter_max=n_iter(2)
			if(n_iter_max<n_iter(1))n_iter_max=n_iter(1)
		else
			if(isnan(ener)) n_err=n_err + 1

			do ind=1,nat3d
		        	xat(ind)=xsave(ind)
		        	yat(ind)=ysave(ind)
		        	zat(ind)=zsave(ind)
		        enddo

			if(isnan(ener)) write (*,*) 'do not accept: E=nan: ',n_debug

		endif
	enddo
	
	do ind=1,nat3d
	   xin(ind)=xat(ind)*arete(1)
           yin(ind)=yat(ind)*arete(1)
           zin(ind)=zat(ind)*arete(1)
        enddo
	
	call comparison(lt1neigh,lt3neigh)!cicle again if 
        !there are atoms with less than 3 nearest neighbors
	!i parametri d'ordine strutturali vengono assegnati da struttura.f90:
	call struttura(xin,yin,zin)
        !il parametro d'ordine invece e' assegnato qui in drive.f90:
	!if the energy is the order parameter, its value is not assigned by struttura:
!        if(order_parameter(1) .eq. 5)then
!          parameter_value(1)=ener
!	elseif(order_parameter(2) .eq. 5)then
!          parameter_value(2)=ener	
!        endif
        call force_pres
!        open(44,file="pressures.xyz",status='unknown')
!        write (44,*) nat3d
!        write (44,*) nomemet(1), ' ',nomemet(2), '  ', ener
!	do indpres=1,nat3d
!	   write (44,'(a2,3x,4f18.8)') spec(indpres),xin(indpres),yin(indpres),zin(indpres),&
!	   & pres_atom(indpres)
!	enddo
!	close(44)
!S	stop
        !the candidate has energy and order parameter...:
	walker_parameter(iwalk,:)=parameter_value(:)
        enerwalker(iwalk)=ener
        esubwalker(iwalk)=e_met_mgo
!                write (*,*) "after esubwalker, ind_mc=",ind_mc, e_met_mgo
        
	!****************
        !SCELTA ALGORITMO
	!****************
        epsihisto=0.000001 !used by mew2 algorithm when perc=perc_old
        
        SELECT CASE (choose_algorithm)
	
	CASE ('_pew_')
           !**********************************************************************
	   !ALGORITMO PEW
	   !**********************************************************************
           !questo algoritmo fa correre n walkers contemporaneamente e li respinge
           !di una quantita' fissata se hanno un parametro d'ordine vicino
           if((n_neigh_old .gt. 0).and.(string_move .ne. 'exch'))then
              deltaen=ener-ener_old-repwalkers
           else
              deltaen=ener-ener_old
           endif
	CASE ('_pew1')
           !**********************************************************************
	   !ALGORITMO PEW1
	   !**********************************************************************
           !questo algoritmo fa correre n walkers contemporaneamente e li respinge
           !di una quantita' fissata se hanno un parametro d'ordine vicino
           if((n_neigh_old .gt. 0).and.(maxenw.eq.1)) then
              deltaen=ener-ener_old-repwalkers
           else
              deltaen=ener-ener_old
           endif

        CASE ('_mew1')
	   !**********************************************************************
	   !ALGORITMO MEW1
	   !**********************************************************************
           !questo algoritmo attribuisce un vantaggio fissato (repwalkers) all'uscita da 
           !configurazioni con una barra piu' alta di quella di arrivo
	   
	   !looking for the height of the bars corresponding to the old and new order parameters:
           do ind=0,nbars(1)
	   do ind2=ind,nbars(2)
              hi_l(1)=hist_step(1)*ind+start_histo(1)
              hi_s(1)=hi_l(1)+hist_step(1)
	      hi_l(2)=hist_step(2)*ind+start_histo(2)
              hi_s(2)=hi_l(2)+hist_step(2)
              if (((walker_parameter(iwalk,1).ge.hi_l(1)).and.(walker_parameter(iwalk,1).lt.hi_s(1))).and.&
	          &(walker_parameter(iwalk,2).ge.hi_l(2)).and.(walker_parameter(iwalk,2).lt.hi_s(2))) then
                 perc=perc_histo(ind+1,ind2+1) !bar height for the new minimum
              endif
              if (((walker_parameter_old(1).ge.hi_l(1)).and.(walker_parameter_old(1).lt.hi_s(1))).and.&
	          &(walker_parameter_old(2).ge.hi_l(2)).and.(walker_parameter_old(2).lt.hi_s(2))) then
                 perc_old=perc_histo(ind+1,ind2+1) !bar height for the old minimum
              endif
	   enddo
           enddo
           
           if(perc_old .gt. perc) then
              deltaen=ener-ener_old-repwalkers
           else
              deltaen=ener-ener_old
           endif
        CASE('_mew2')
           !**********************************************************************
	   !ALGORITMO MEW2
	   !**********************************************************************
           !questo algoritmo attribuisce un vantaggio fissato (repwalkers) all'uscita da 
           !configurazioni con una barra piu' alta di quella di arrivo, e uno svantaggio
           !pari al vantaggio se e' vero il contrario
           !looking for the height of the bars corresponding to the old and new order parameter:
           do ind=0,nbars(1)
	   do ind2=ind,nbars(2)
              hi_l(1)=hist_step(1)*ind+start_histo(1)
              hi_s(1)=hi_l(1)+hist_step(1)
	      hi_l(2)=hist_step(2)*ind+start_histo(2)
              hi_s(2)=hi_l(2)+hist_step(2)
              if (((walker_parameter(iwalk,1).ge.hi_l(1)).and.(walker_parameter(iwalk,1).lt.hi_s(1))).and.&
	          &(walker_parameter(iwalk,2).ge.hi_l(2)).and.(walker_parameter(iwalk,2).lt.hi_s(2))) then
                 perc=perc_histo(ind+1,ind2+1) !bar height for the new minimum
              endif
              if (((walker_parameter_old(1).ge.hi_l(1)).and.(walker_parameter_old(1).lt.hi_s(1))).and.&
	          &(walker_parameter_old(2).ge.hi_l(2)).and.(walker_parameter_old(2).lt.hi_s(2))) then
                 perc_old=perc_histo(ind+1,ind2+1) !bar height for the old minimum
              endif
	   enddo
           enddo
           
           if(perc_old .gt. (perc+epsihisto))then
              deltaen=ener-ener_old-repwalkers
           elseif(perc_old .lt. (perc-epsihisto))then
              deltaen=ener-ener_old+repwalkers
           else
              deltaen=ener-ener_old
           endif
        CASE('histo')
	   !**********************************************************************
	   !ALGORITMO HISTO
	   !**********************************************************************
           !questo algoritmo attribuisce un vantaggio/svantaggio all'uscita pesato sulla
           !differenza tra l'altezza delle barre
           !looking for the height of the bars corresponding to the old and new order parameter:

           do ind=0,nbars(1)-1
	   do ind2=0,nbars(2)-1
              hi_l(1)=hist_step(1)*ind+start_histo(1)
              hi_s(1)=hi_l(1)+hist_step(1)
	      hi_l(2)=hist_step(2)*ind2+start_histo(2)
              hi_s(2)=hi_l(2)+hist_step(2)
	      
              if (((walker_parameter(iwalk,1).ge.hi_l(1)).and.(walker_parameter(iwalk,1).lt.hi_s(1))).and.&
	          &(walker_parameter(iwalk,2).ge.hi_l(2)).and.(walker_parameter(iwalk,2).lt.hi_s(2))) then
		  !write(*,*)'walker:',walker_parameter(iwalk,1),hi_l(1),hi_s(1)
                 perc=perc_histo(ind+1,ind2+1) !bar height for the new minimum
!		if (ind_mc .gt. 100)then
!	   	write(*,*)'perc: ',perc
!	   	endif
              endif
              if (((walker_parameter_old(1).ge.hi_l(1)).and.(walker_parameter_old(1).lt.hi_s(1))).and.&
	          &(walker_parameter_old(2).ge.hi_l(2)).and.(walker_parameter_old(2).lt.hi_s(2))) then
		  !write(*,*)'old walker:',walker_parameter(iwalk,1),hi_l(1),hi_s(1)
                 perc_old=perc_histo(ind+1,ind2+1) !bar height for the old minimum
!		 if (ind_mc .gt. 100)then
!	   	write(*,*)'perc_old: ',perc_old
!	   	endif
              endif
	   enddo
           enddo
           
	   !write(*,*)walker_parameter(iwalk,1)
	   !write(*,*)perc,perc_old,weight,weight*(perc-perc_old)
!	   if (ind_mc .gt. 100)then
!	   write(*,*)'ener-ener_old',ener-ener_old
!	   endif
           deltaen=ener-ener_old+weight*(perc-perc_old)
!	   if (ind_mc .gt. 100)then
!	   write(*,*)'deltaen',deltaen
!	   endif
	   
	CASE('_hexc')
	   !**********************************************************************
	   !ALGORITMO HISTO_EXCAPE
	   !**********************************************************************
           !questo algoritmo attribuisce solo un vantaggio all'uscita pesato sulla
           !differenza tra l'altezza delle barre
           !looking for the height of the bars corresponding to the old and new order parameter:
           do ind=0,nbars(1)
	   do ind2=ind,nbars(2)
              hi_l(1)=hist_step(1)*ind+start_histo(1)
              hi_s(1)=hi_l(1)+hist_step(1)
	      hi_l(2)=hist_step(2)*ind+start_histo(2)
              hi_s(2)=hi_l(2)+hist_step(2)
              if (((walker_parameter(iwalk,1).ge.hi_l(1)).and.(walker_parameter(iwalk,1).lt.hi_s(1))).or.&
	          &(walker_parameter(iwalk,2).ge.hi_l(2)).and.(walker_parameter(iwalk,2).lt.hi_s(2))) then
                 perc=perc_histo(ind+1,ind2+1) !bar height for the new minimum
              endif
              if (((walker_parameter_old(1).ge.hi_l(1)).and.(walker_parameter_old(1).lt.hi_s(1))).or.&
	          &(walker_parameter_old(2).ge.hi_l(2)).and.(walker_parameter_old(2).lt.hi_s(2))) then
                 perc_old=perc_histo(ind+1,ind2+1) !bar height for the old minimum
              endif
	   enddo
           enddo
	   
	   if(perc_old .gt. perc) then
              deltaen=ener-ener_old+weight*(perc-perc_old)
	   else
	      deltaen=ener-ener_old
	   endif
	   
	  CASE('misto')
	   !**********************************************************************
	   !ALGORITMO MISTO 
	   !**********************************************************************
           !questo algoritmo e' una formula mista di histo e pew.
           !Se un walker ha vicini, si scalda secondo quanto previsto dall'agoritmo pew (repwalkers)
           !Se invece non ha alcun vicino, si scalda secondo quanto previsto dall'algoritmo
	   !histo (weight).
           do ind=0,nbars(1)
	   do ind2=ind,nbars(2)
              hi_l(1)=hist_step(1)*ind+start_histo(1)
              hi_s(1)=hi_l(1)+hist_step(1)
	      hi_l(2)=hist_step(2)*ind+start_histo(2)
              hi_s(2)=hi_l(2)+hist_step(2)
              if (((walker_parameter(iwalk,1).ge.hi_l(1)).and.(walker_parameter(iwalk,1).lt.hi_s(1))).and.&
	          &(walker_parameter(iwalk,2).ge.hi_l(2)).and.(walker_parameter(iwalk,2).lt.hi_s(2))) then
                 perc=perc_histo(ind+1,ind2+1) !bar height for the new minimum
              endif
              if (((walker_parameter_old(1).ge.hi_l(1)).and.(walker_parameter_old(1).lt.hi_s(1))).and.&
	          &(walker_parameter_old(2).ge.hi_l(2)).and.(walker_parameter_old(2).lt.hi_s(2))) then
                 perc_old=perc_histo(ind+1,ind2+1) !bar height for the old minimum
              endif
	   enddo
           enddo
	              
	   if(n_neigh_old .gt. 0)then ! il walker ha vicini, agisce pew
              deltaen=ener-ener_old-repwalkers
           else                       ! il walker non ha vicini, agisce histo
              deltaen=ener-ener_old+weight*(perc-perc_old)
           endif

        END SELECT
	!keep the energy of the trial configuration
        ener_att=ener


	!deciding whether accepting or not the move:
	if((-deltaen/(cbol*temperatura)) .gt. 0.)then
           boltz_fact=1.
	elseif((-deltaen/(cbol*temperatura)) .lt. 0)then	        
           boltz_fact=exp(-deltaen/(cbol*temperatura))
	elseif((-deltaen/(cbol*temperatura)) .eq. 0)then	
           boltz_fact=1.
        endif
	
        call ranmar(zzz)
        accflag=0
        if ((zzz.le.boltz_fact).and.(lt1neigh.eq.0)) then !if the move is accepted:
           !write(*,*)'move is accepted at step',ind_tot
	   !updating of moves counters:
	   accflag=1
	   n_accepted_moves=n_accepted_moves+1
	   if(dabs(ener-ener_old) .gt. 0.005d0)then
	      n_true_moves=n_true_moves+1
	      if(lista_minimi)call minima_collecting(ener,walker_parameter(iwalk,1),ind_mc)
	   endif
	   !coordinates are updated:
           do ind=1,nat3d
              xpos(ind,iwalk)=xin(ind)
              ypos(ind,iwalk)=yin(ind)
              zpos(ind,iwalk)=zin(ind)
              specpos(ind,iwalk)=spec(ind)
           enddo
           
	   !histogram is updated:
           call calcola_histo
           
	   !if the current lowest minimum has been found:
	   if (ener.lt.enermin) then
!	   if(ener .lt. 0)then
!write(*,*)'call output index=',ind_tot,' n_mc=',n_mc,' nglo=',nglo
	      nglo=nglo+1
	      
	      ncall_output=ncall_output+1
	      label_output_coordinates=.true.
	      label_output_coordinates2=.true. !this flag is used by standard_bh
              output_instruction=2
	      call output_glo(iwalk,output_instruction)
	      enermin=ener
           endif !su ener.lt enermin
	   	   
        else !if the move is NOT accepted:
           
	   !write(*,*)ind_mc, 'rifiuto'
	   n_rejected_moves=n_rejected_moves+1
           ener=ener_old
           e_met_mgo=esub_old
	   enerwalker(iwalk)=ener_old
           esubwalker(iwalk)=esub_old
	   walker_parameter(iwalk,1)=walker_parameter_old(1)
	   walker_parameter(iwalk,2)=walker_parameter_old(2)
           
	   call calcola_histo !histogram is updated with the old data
           
        endif ! accepting / not accepting
        
	!writing output at every MC step:
	ncall_output=ncall_output+1
	label_output_coordinates=.false.
        output_instruction=3
!        write (*,*) "ind_mc=", ind_mc, ener, e_met_mgo
	call output_glo(iwalk,output_instruction)
        !Every 1000 steps, percentages of acceptance/rejection of moves are recorded:
        If (mod(ind_tot,1000) .eq. 0) Then
            output_instruction=4
            call output_glo(iwalk,output_instruction)
        Endif
	
	!OPTION RESTART
        if(restart .eq. 'si')then
           if(ener .lt. ener_restart)then
	      nrestart=nrestart+1
	      write(*,*)'Step=',ind_tot
	      write(*,*)'Ener=',ener
	      write(*,*)'ener_restart=',ener_restart
	      write(*,*)'Il programma riparte per la',nrestart,'-esima volta'	 
              
	      !scrittura nel file successi.out
	      if(nrestart .eq. 1)then
	        open(300,file='successi.out',status='unknown')
	        write(300,*)'#MCstep, ener, ripartenza num.'
	        call flush(300)
	        close(300)
	      endif
	      open(300,file='successi.out',status='unknown',position='append')
	      write(300,'(i10,f16.8,i5)')ind_tot,ener,nrestart
	      write(300,*)'number of disallowed steps:',n_err
	      call flush(300)
	      close(300)
              !optimization starts again from nwalkers random conenfigurations:
	      first_minimization = .true.
              call initialization   
!	      first_minimization = .false.
	      steps_sincelast_gm=0
           endif
        endif

  enddo ! su ind_mc
  
  return
  
END SUBROUTINE DRIVE_BH

