SUBROUTINE READ_BH

  USE PARAMETERS
  USE BH_MOD
  USE HISTO
  USE FUNCTION_PAR, only : p,q,a,qsi,ecoh,rat,mass,cutoff,cutoff_start,cutoff_end,&
                           ddd,amgo,cutz,znmax,arete,dist,&
                           a5,a4,a3,x5,x4,x3
  IMPLICIT NONE

  !Local variables
  Integer :: ind,i,j,k
  Real(8) :: ar,br,cr,ab,bb,cb,dik0,u2,v2,zr,zrp,za,zap,frac
  Real(8) :: rmet1,rmet2,cont_big,cont_small
  Character(2) :: metal(2)
  Character(20):: param_file
 
  open (22,file='errors_bh.out',status='unknown')
  write(22,*)'ERRORS AND WARNINGS'
  close(22)

  write (*,*) 'START READING FILE input_bh.in'
  open(19,file='input_bh.in',status='old')

  ! choice between seeded and unseeded search
  read(19,*)!******Where do we start from?
  read(19,*) choice_seeding
  IF (choice_seeding.eq.'seeded') THEN !SEEDED RUN
     write (*,*) 'START READING FILE seed.in'
     open (21,file='seed.in',status='old')
     read (21,*) nat3d
     if (nat3d.gt.maxpar) then
        write (*,*) ' the cluster is too big'
        open(22,file='errors_bh.out',status='unknown',position='append')
        write (22,*) ' the cluster is too big: nat3d > maxpar'
        close(22)
        stop
     endif
     read (21,*) nomemet(1),nomemet(2)
     nspec1=0
     nspec2=0
     do ind=1,nat3d
        read (21,*) spec(ind),xin(ind),yin(ind),zin(ind)
        if (spec(ind).eq.nomemet(1)) then
           nspec1=nspec1+1
        else if (spec(ind).eq.nomemet(2)) then
           nspec2=nspec2+1
        else
           write (*,*) 'the species of atom ',ind, ' does not match'
           open (22,file='errors_bh.out',status='unknown',position='append')
           write (22,*) 'the species of atom ',ind, ' does not match'
           close(22)
           stop
        endif
     enddo
     close(21)
     write (*,*) 'END READING FILE seed.in'
  ELSE ! UNSEEDED RUN
     write (*,*) 'START READING FILE noseed_rgl.in'
     open(21,file='noseed_rgl.in',status='old')
     read (21,*) nat3d
     read (21,*) nomemet(1),nomemet(2)
     write(*,*)'in noseed species are:',nomemet(1),nomemet(2)
     read (21,*) nspec1,nspec2 
     if (nspec1+nspec2.ne.nat3d) then
        write (*,*) ' bad atoms count in noseed.in'
        open (22,file='errors_bh.out',status='unknown',position='append')
        write (22,*) ' bad atoms count in noseed.in'
        close(22)
        stop
     endif

     !Assign species to atoms
     do ind=1,nspec1
        spec(ind)=nomemet(1)
     enddo
     do ind=nspec1+1,nat3d
        spec(ind)=nomemet(2)
     enddo
     
	read (21,*)
        read (21,*) x_min,x_max
        read (21,*) y_min,y_max
        read (21,*) z_min,z_max
      close(21)
     write (*,*) 'END READING FILE noseed_rgl.in'
  ENDIF
  read (19,*) refresh !refresh or not refresh the seed for random generators
  read(19,*)!******Interaction potentials
  read (19,*) cutoff !cutoff or no cutoff for the metal-metal interaction
  read (19,*) param_file

  !STEPS AND TEMPERATURE:
  read(19,*)!******Run settings
  read (19,*) n_mc ! number of Monte Carlo steps

  !CHOICE AND TUNING OF THE MOVE
  read(19,*)!******Moves
  read (19,*) probmove1,temp_bonds ! probability of choosing the 'bonds' move
  read (19,*) probmove2,temp_ball ! probability of choosing the 'ball' move
  read (19,*) probmove3,shell_thickness,temp_shell ! probability of choosing the 'shell' move
  read (19,*) probmove4,rho_inf_sk,rho_sup_sk,temp_shake ! probability of choosing the 'shake' move
  read (19,*) probmove5,rho_inf_he,rho_sup_he,temp_highener ! probability of choosing the 'highenergy' move
  read (19,*) probmove6,rho_inf_sg,rho_sup_sg,temp_single
  read (19,*) probmove7,npas_bro,temp_brow,temp_brow_accept ! probability of choosing the 'brownian' move
  read (19,*) probmove8,npas_brosurf,temp_browsurf,temp_browsurf_accept ! probability of choosing the 'brownian on the surface' move
  read (19,*) probmove9,zcut,deltaz,temp_mirror
  read (19,*) probmove10 ! probability of choosing the 'exchange' move
  !sub-probabilities for the various kinds of exchange moves:
  read (19,*) prob_exch_all,temp_all ! probability of choosing the 'exchange all' move
  read (19,*) prob_exch_sur,temp_sur ! probability of choosing the 'exchange surface' move
  read (19,*) prob_exch_csh,bignn_csh,smallnn_csh,temp_csh ! probability of choosing the 'exchange core-shell' move
  read (19,*) prob_exch_int,bignn_int,smallnn_int,temp_int ! probability of choosing the 'exchange interface' move
  read (19,*) prob_exch_sep,exp_sep,temp_sep,small_inside ! probability of choosing the 'exchange separation' move
  read (19,*) prob_exch_mix,exp_mix,temp_mix ! probability of choosing the 'exchange mixing' move
  read (19,*) prob_exch_mult,temp_mult,exp_mult
  read (19,*) prob_exch_gru

  read(19,*)!*******Order parameters
  !CHOICE OF THE ORDER PARAMETER
  read (19,*) order_parameter(1),order_parameter(2)
  read(19,*)!*******Choice of the GO algorithm and of the number of walkers
  read (19,'(a5)') choose_algorithm !algorithm
  read (19,*) nwalkers ! number of walkers
  read(19,*)!*******Rules for selecting the walkers  
  read (19,*) walkers_turnover ! turn-over algorithm to be applied to walkers turnover=1 is sequential, turnover=2 is random
  read (19,*) freq_walker_exchange ! frequency of walker turnover, if turnover=3
  read(19,*)!******* Parameters for the interactions between walkers
  read (19,*) repwalkers ! walkers' repulsive interaction
  read (19,*) widthwalkers(1),widthwalkers(2) ! width of the repulsive interactions
  read(19,*)!******* Parameters for the histograms 
  read (19,*) weight ! weight of the histogram
  read (19,*) start_histo(1) !starting value of the first order parameter:
  read (19,*) hist_step(1) !step of the first histogram
  read (19,*) nbars(1) ! number of intervals of the first histogram
  read (19,*) start_histo(2)  !starting value of the second order parameter:
  read (19,*) hist_step(2) !step of the second histogram
  read (19,*) nbars(2) ! number of intervals of the second histogram
  read (19,*) renew_histo ! refresh frequency of the histogram
  read(19,*)!******Output instructions
  read (19,*)output_walkers !want different output for different walkers?
  read(19,*)lista_minimi !want a list of the minima encountered during optimization?
  read(19,*)param_min,param_max !order parameter range where collecting minima
  read(19,*)parameter_partitioning !order parameter step to distinguish minima  
  read(19,*)!******Option for restarting the simulation    
  read (19,'(a2)')restart !restart
  read(19,*)ener_restart !energy threshold for restart [eV]
  close(19)
  write (*,*) 'END READING FILE input_bh.in'

  write (*,*) 'START READING FILE ', param_file
  open(19,file=param_file,status='old')
  read(19,*)
  read(19,*)
  read(19,*)metal(1), metal(2)
            !check consistency with metals declared in seed.in
            if ((metal(1) .ne. nomemet(1)) .and. (metal(2) .ne. nomemet(2)))then
               write(*,*)'####################################################'
               write(*,*)'ERROR:'
               write(*,*)'Metals declared in seed.in or noseed_rgl.in'
               write(*,*)'MUST BE THE SAME AND IN THE SAME ORDER' 
               write(*,*)'as in the parametrization file.'
               write(*,*)'####################################################'
               open(22,file='errors_bh.out',status='unknown',position='append')
               write(22,*)'####################################################'
               write(22,*)'ERROR:'
               write(22,*)'Metals declared in seed.in or noseed_rgl.in'
               write(22,*)'MUST BE THE SAME AND IN THE SAME ORDER' 
               write(22,*)'as in the parametrization file.'
               write(22,*)'####################################################'
               close(22)
               stop
            endif
            !end check consistency
  read(19,*)
  read(19,*)
  read(19,*)p(1),p(2),p(3)
  read(19,*)q(1),q(2),q(3)
  read(19,*)a(1),a(2),a(3)
  read(19,*)qsi(1),qsi(2),qsi(3)
  read(19,*)
  read(19,*)
  read(19,*)ecoh(1),ecoh(2)
  read(19,*)rat(1),rat(2)
  read(19,*)mass(1),mass(2)
  read(19,*)
  read(19,*)
  read(19,*)cutoff_start(1),cutoff_end(1)
  read(19,*)cutoff_start(2),cutoff_end(2)
  read(19,*)cutoff_start(3),cutoff_end(3)
  read(19,*)
  read(19,*)
  read(19,*)(((ddd(1,i,j,k),i=1,3),j=1,3),k=1,3)
  read(19,*)
  read(19,*)
  read(19,*)(((ddd(2,i,j,k),i=1,3),j=1,3),k=1,3)
  read(19,*)
  read(19,*)
  read(19,*)amgo(1),amgo(2)
  read(19,*)
  read(19,*)
  read(19,*)
  read(19,*)cutz(1,1),cutz(2,2),cutz(1,2),cutz(2,1)
  read(19,*)znmax(1),znmax(2)
  write (*,*) 'END READING FILE ', param_file
  close(19)
  ! flag for move_mirror
  zcutflag=0
  if (zcut.gt.1.d0) zcutflag=1
  !Unit conversions

  arete(1)=rat(1)*dsqrt(8.d0)
  arete(2)=rat(2)*dsqrt(8.d0)
  arete(3)=(arete(1)+arete(2))/2.0d0

  !nn are the nearest neighbours distances in A
  nn(1)=arete(1)/dsqrt(2.d0)
  nn(2)=arete(2)/dsqrt(2.d0)
  nn(3)=arete(3)/dsqrt(2.d0)

  !dist are the nearest neighbours distances in arete(1) units
  dist(1)=1.d0/dsqrt(2.d0)
  dist(2)=nn(2)/arete(1)
  dist(3)=nn(3)/arete(1)

  !cutoff_end and cutoff_start are converted into arete(1) units
  do i=1,3
  if(cutoff .eq. 'si')then
     cutoff_start(i)=cutoff_start(i)/arete(1)
     cutoff_end(i)=cutoff_end(i)/arete(1)
  elseif(cutoff .eq. 'no')then
     cutoff_start(i)=2000.d0
     cutoff_end(i)=2000.d0
  else
     open (22,file='errors_bh.out',status='unknown',position='append')
     write(22,*)'ERROR:'
     write(22,*)'cutoff must be =`si` or =`no`'
     write(22,*)'PROGRAM EXIT'
     close(22)
     stop
  endif
  enddo
  !cutz is converted into arete(1) units
  cutz(1,1)=cutz(1,1)/arete(1)
  cutz(1,2)=cutz(1,2)/arete(1)
  cutz(2,1)=cutz(2,1)/arete(1)
  cutz(2,2)=cutz(2,2)/arete(1)
  !cutz will be compared to square distances in bigvoi.f90:
  cutz(1,1)=cutz(1,1)*cutz(1,1)
  cutz(2,2)=cutz(2,2)*cutz(2,2)
  cutz(1,2)=cutz(1,2)*cutz(1,2)
  cutz(2,1)=cutz(2,1)*cutz(2,1)

  !Cutoff parameters a5,a4,a3,x5,x4,x3
  do i=1,3 
        dik0=dist(i)
        ar=-a(i)*dexp(-p(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**3 
        br=-(p(i)/dik0)*a(i)*dexp(-p(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**2
        cr=-((p(i)/dik0)**2) &
             *a(i)*dexp(-p(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))
        ab=-qsi(i)*dexp(-q(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**3
        bb=-(q(i)/dik0)*qsi(i)*dexp(-q(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))**2 
        cb=-((q(i)/dik0)**2) &
             *qsi(i)*dexp(-q(i)*(cutoff_start(i)/dik0-1.d0))/(cutoff_end(i)-cutoff_start(i))
        x5(i)=(12.d0*ab-6.d0*bb+cb)/(2.d0*(cutoff_end(i)-cutoff_start(i))**2)
        x4(i)=(15.d0*ab-7.d0*bb+cb)/(cutoff_end(i)-cutoff_start(i))
        x3(i)=(20.d0*ab-8.d0*bb+cb)/2.d0
        a5(i)=(12.d0*ar-6.d0*br+cr)/(2.d0*(cutoff_end(i)-cutoff_start(i))**2)
        a4(i)=(15.d0*ar-7.d0*br+cr)/(cutoff_end(i)-cutoff_start(i))
        a3(i)=(20.d0*ar-8.d0*br+cr)/2.d0
     enddo

  ! If two species are present, two lists of atoms are created: list_big and list_small, 
  !containing the list of indexes corresponding to the atoms with the big and small radius, respectively.
  IF(nomemet(1) .ne. nomemet(2))THEN
    If(rat(1) .lt. rat(2))Then
       nsmall=nspec1
       nbig=nspec2
    ElseIf(rat(1) .gt. rat(2))Then
       nsmall=nspec2
       nbig=nspec1
    EndIf
    !write(*,*)'Test: nbig=',nbig,'nsmall=',nsmall
    
    !allocate(list_big(nbig))
    !allocate(list_small(nsmall))
    
    cont_small=0
    cont_big=0
    do ind=1,nat3d
       if((spec(ind).eq.nomemet(1)) .and. (rat(1) .lt. rat(2)))then!the ind atom is small
         cont_small=cont_small+1
         list_small(cont_small)=ind 
       elseif((spec(ind).eq.nomemet(2)) .and. (rat(2) .lt. rat(1)))then!the ind atom is small
         cont_small=cont_small+1
         list_small(cont_small)=ind
       elseif((spec(ind).eq.nomemet(1)) .and. (rat(1) .gt. rat(2)))then!the ind atom is big
         cont_big=cont_big+1
         list_big(cont_big)=ind
       elseif((spec(ind).eq.nomemet(2)) .and. (rat(2) .gt. rat(1)))then!the ind atom is big
         cont_big=cont_big+1
         list_big(cont_big)=ind
       endif
    enddo
!  write(*,*)list_big(:)
!  write(*,*)'*****'
!  write(*,*)list_small(:)
  ENDIF


  ! writes everything in output_bh.out file
  open (22,file='output_bh.out',status='unknown')
  write (*,*)'###############################################################################'
  write (*,*)'RUN CONFIGURATION'
  write (*,*)'###############################################################################'
  write (22,*)'###############################################################################'
  write (22,*)'RUN CONFIGURATION'
  write (22,*)'###############################################################################'
  write (*,*) ' type of optimization: ',choice_seeding
  write (22,*) ' type of optimization: ',choice_seeding
  write (*,*) ' refresh random numbers: ',refresh
  write (22,*) ' refresh random numbers: ',refresh
  write (*,*) ' cut-off: ',cutoff
  write (22,*) ' cutoff: ',cutoff
  write (*,*) ' parametrization: ',param_file
  write (22,*) ' parametrization: ',param_file
  write (*,*)'##############################################################################'
  write (*,*)'MC STEPS'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'MC STEPS'
  write (22,*)'##############################################################################'
  write (*,*) ' number of MC steps: ',n_mc
  write (22,*) ' number of MC steps: ',n_mc
  write (*,*)'##############################################################################'
  write (*,*)'MOVES'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'MOVES'
  write (22,*)'##############################################################################'
  write (*,*) 'probability of choosing move n.1: bonds ',probmove1
  write (22,*) 'probability of choosing move n.1: bonds ',probmove1
  write (*,*) 'probability of choosing move n.2: ball ',probmove2
  write (22,*) 'probability of choosing move n.2: ball ',probmove2
  write (*,*) 'probability of choosing move n.3: shell ',probmove3
  write (22,*) 'probability of choosing move n.3: shell ',probmove3
  write (*,*) 'probability of choosing move n.4: shake ',probmove4
  write (22,*) 'probability of choosing move n.4: shake ',probmove4
  write (*,*) 'probability of choosing move n.5: high energy ',probmove5
  write (22,*) 'probability of choosing move n.5: high energy ',probmove5
  write (*,*) 'probability of choosing move n.6: single ',probmove6
  write (22,*) 'probability of choosing move n.6: single ',probmove6
  write (*,*) 'probability of choosing move n.7: brownian ',probmove7
  write (22,*) 'probability of choosing move n.7: brownian ',probmove7
  write (*,*) 'probability of choosing move n.8: brownian surf ',probmove8
  write (22,*) 'probability of choosing move n.8: brownian surf ',probmove8
  write (*,*) 'probability of choosing move n.9: mirror ',probmove9
  write (22,*) 'probability of choosing move n.9: mirror ',probmove9
  write (*,*) 'probability of choosing move n.10: exchange ',probmove10
  write (22,*) 'probability of choosing move n.10: exchange ',probmove10
   write (*,*) 'probability of choosing the exchange move "all": ',prob_exch_all
  write (22,*) 'probability of choosing the exchange move "all"',prob_exch_all
  write (*,*) 'probability of choosing the exchange move "sur": ',prob_exch_sur
  write (22,*) 'probability of choosing the exchange move "sur": ',prob_exch_sur
  write (*,*) 'probability of choosing the exchange move "csh": ',prob_exch_csh
  write (22,*) 'probability of choosing the exchange move "csh": ',prob_exch_csh
  write (*,*) 'probability of choosing the exchange move "int": ',prob_exch_int
  write (22,*) 'probability of choosing the exchange move "int": ',prob_exch_int
  write (*,*) 'probability of choosing the exchange move "sep": ',prob_exch_sep
  write (22,*) 'probability of choosing the exchange move "sep": ',prob_exch_sep
  write (*,*) 'probability of choosing the exchange move "mix": ',prob_exch_mix
  write (22,*) 'probability of choosing the exchange move "mix": ',prob_exch_mix
  write (*,*) 'probability of choosing the exchange move "GRU": ',prob_exch_gru
  write (22,*) 'probability of choosing the exchange move "GRU": ',prob_exch_gru
  write (*,*)'##############################################################################'
  write (*,*)'ORDER PARAMETERS'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'ORDER PARAMETERS'
  write (22,*)'##############################################################################'
  write (*,*) ' order parameters:',order_parameter(1),order_parameter(2)
  write (22,*) ' order parameters:',order_parameter(1),order_parameter(2)
  write (*,*)'##############################################################################'
  write (*,*)'CHOICE OF THE ALGORITHM AND OF THE NUMBER OF WALKERS'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'CHOICE OF THE ALGORITHM AND OF THE NUMBER OF WALKERS'
  write (22,*)'##############################################################################'
  write (*,*) ' algorithm: ',choose_algorithm
  write (22,*) ' algorithm: ',choose_algorithm
  write (*,*) ' number of walkers: ',nwalkers
  write (22,*) ' number of walkers: ',nwalkers
  write (*,*)'##############################################################################'
  write (*,*)'RULES FOR SELECTING THE WALKERS'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'RULES FOR SELECTING THE WALKERS'
  write (22,*)'##############################################################################'
  write (*,*) ' turn-over algorithm to be applied to walkers: ',walkers_turnover
  write (22,*) ' turn-over algorithm to be applied to walkers: ',walkers_turnover
  write (*,*) ' frequency of walker turnover, if turnover=3: ',freq_walker_exchange
  write (22,*) ' frequency of walker turnover, if turnover=3: ',freq_walker_exchange
  write (*,*)'##############################################################################'
  write (*,*)'PARAMETERS FOR THE INTERACTIONS BETWEEN WALKERS'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'PARAMETERS FOR THE INTERACTIONS BETWEEN WALKERS'
  write (22,*)'##############################################################################'
  write (*,*) ' walkers` repulsive interaction: ',repwalkers
  write (22,*) ' walkers` repulsive interaction: ',repwalkers
  write (*,*) ' widths of the repulsive interactions: ', widthwalkers(1),widthwalkers(2)
  write (22,*) ' widths of the repulsive interactions: ', widthwalkers(1),widthwalkers(2)
  write (*,*)'##############################################################################'
  write (*,*)'PARAMETERS FOR THE HISTOGRAM'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'PARAMETERS FOR THE HISTOGRAM'
  write (22,*)'##############################################################################'
  write (*,*) 'weight of the histogram: ',weight
  write (22,*) 'weight of the histogram: ',weight
  write (*,*) ' starting value of the first order parameter: ',start_histo(1)
  write (22,*) ' starting value of the first order parameter: ',start_histo(1)
  write (*,*) ' step of the first histogram: ',hist_step(1)
  write (22,*) ' step of the first histogram: ',hist_step(1)
  write (*,*) ' number of intervals of the first histogram: ',nbars(1)
  write (22,*) ' number of intervals of the first histogram: ',nbars(1)
  write (*,*) ' starting value of the second order parameter: ',start_histo(2)
  write (22,*) ' starting value of the second order parameter: ',start_histo(2)
  write (*,*) ' step of the second histogram: ',hist_step(2)
  write (22,*) ' step of the second histogram: ',hist_step(2)
  write (*,*) ' number of intervals of the second histogram: ',nbars(2)
  write (22,*) ' number of intervals of the second histogram: ',nbars(2)
  write (*,*) ' refresh frequency of the histogram: ',renew_histo
  write (22,*) ' refresh frequency of the histogram: ',renew_histo
  write (*,*)'##############################################################################'
  write (*,*)'OUTPUT'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'OUTPUT'
  write (22,*)'##############################################################################'
  write (*,*) ' want different output for different walkers? ',output_walkers
  write (22,*) ' want different output for different walkers? ',output_walkers
  write (*,*) ' want a list of the minima encountered during optimization?',lista_minimi
  write (22,*) ' want a list of the minima encountered during optimization?',lista_minimi
  write (*,*) ' order parameter range where collecting minima: ',param_min,param_max 
  write (22,*) ' order parameter range where collecting minima: ',param_min,param_max 
  write (*,*) ' order parameter step to distinguish minima: ',parameter_partitioning
  write (22,*) ' order parameter step to distinguish minima: ',parameter_partitioning
  write (*,*)'##############################################################################'
  write (*,*)'OPTION FOR RESTARTING THE SIMULATION'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'OPTION FOR RESTARTING THE SIMULATION'
  write (22,*)'##############################################################################'
  write (*,*) ' restart: ',restart
  write (22,*) ' restart: ',restart
  write (*,*) ' energy threshold for restart [eV]: ',ener_restart
  write (22,*) ' energy threshold for restart [eV]: ',ener_restart
  write (*,*)'##############################################################################'
  write (*,*)'CLUSTER'
  write (*,*)'##############################################################################'
  write (22,*)'##############################################################################'
  write (22,*)'CLUSTER'
  write (22,*)'##############################################################################'
  write (*,*) ' number of atoms in the cluster = ', nat3d
  write (22,*) ' number of atoms in the cluster = ', nat3d
  if ((nspec1 .lt. nat3d) .and. (nspec2 .lt. nat3d)) then
     write (*,*) ' bimetallic cluster of ', nomemet(1),'-',nomemet(2)
     write (*,*) ' number of ', nomemet(1), ' atoms = ', nspec1
     write (*,*) ' number of ', nomemet(2), ' atoms = ', nspec2
     write (22,*) ' bimetallic cluster of ', nomemet(1),'-',nomemet(2)
     write (22,*) ' number of ', nomemet(1), ' atoms = ', nspec1
     write (22,*) ' number of ', nomemet(2), ' atoms = ', nspec2
  elseif ((nspec1 .eq. nat3d) .and. (nspec2 .eq. 0))then
     write (*,*) ' monometallic cluster of ', nomemet(1)
     write (*,*) ' number of ', nomemet(1), ' atoms = ', nspec1
     write (22,*) ' monometallic cluster of ', nomemet(1)
     write (22,*) ' number of ', nomemet(1), ' atoms = ', nspec1
  elseif ((nspec2 .eq. nat3d) .and. (nspec1 .eq. 0))then
     write (*,*) ' monometallic cluster of ', nomemet(2)
     write (*,*) ' number of ', nomemet(2), ' atoms = ', nspec2
     write (22,*) ' monometallic cluster of ', nomemet(2)
     write (22,*) ' number of ', nomemet(2), ' atoms = ', nspec2
  endif
  if (choice_seeding.eq.'seeded') then
     write (*,*) 
     write (*,*) 'species and initial coordinates of the seed'
     write (22,*) 
     write (22,*) 'species and initial coordinates of the seed'
     do ind=1,nat3d
        write (*,'(a2,3f15.9)') spec(ind),xin(ind),yin(ind),zin(ind)
        write (22,'(a2,3f15.9)') spec(ind),xin(ind),yin(ind),zin(ind)
     enddo
  else
     write (*,*) 
     write (*,*) ' limits of the box'
        write (*,*) 'x_min,x_max=',x_min,x_max
        write (*,*) 'y_min,y_max=',y_min,y_max
        write (*,*) 'z_min,z_max=',z_min,z_max
     write (22,*) 
     write (22,*) ' limits of the box'
        write (22,*) 'x_min,x_max=',x_min,x_max
        write (22,*) 'y_min,y_max=',y_min,y_max
        write (22,*) 'z_min,z_max=',z_min,z_max
  endif
  close(22)  
  
END SUBROUTINE READ_BH
