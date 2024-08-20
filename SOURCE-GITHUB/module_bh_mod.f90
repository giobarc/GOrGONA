MODULE BH_MOD
  
  USE PARAMETERS
  
  Implicit None
  SAVE
  
  Integer, Parameter :: nmaxwalker=30
  Integer :: nat3d    ! number of atoms in the cluster
  Integer :: nspec1,nspec2 ! number of atoms of the different species
  Integer :: nbig, nsmall ! number of atoms with big and small radius
  !Integer, dimension(:), allocatable :: list_big, list_small
  Integer :: list_big(maxpar),list_small(maxpar)
  Integer :: number_of_steps ! number of steps in the simulation
  !step counter for walkers
  Integer :: contwalker(nmaxwalker)
  Integer :: npas ! total number of quenching steps
  Integer :: nfric ! number of steps in quenching with friction
  Integer :: n_mc, ind_mc ! number of Monte Carlo steps and running index
  Integer :: ind_tot ! counter for Monte Carlo steps (does not depend upon the number of walkers)
  Integer :: noutput ! number of output files
  Integer :: maxenw ! check wheter the chosen walker has maximum energy among its neighbours
  Integer :: nglo ! number of  putative gm (only for escape algorithm)
  Integer :: nmaxfail ! maximum number of refused moves
  Integer :: isucc ! flag for the success of local minimization
  Integer :: n_accepted_moves,n_rejected_moves,n_true_moves,accflag
  Integer :: zcutflag
!  Character*2 :: nomemet1,nomemet2 ! symbols of the elements 1 and 2
  character*2 :: refresh ! refres random number generator (si/no)
  Character*2 :: nomemet(2)
  Character*2 :: spec(maxpar) ! species of the atoms of the cluster
  Character*12 :: choice_seeding ! chooses between seeded and unseeded search
  Character*5  :: string_move ! chooses the kind of move
  Character*2  :: output_walkers
  Character*2  :: restart
  Character*4 :: list_type  
  Character*2 :: specpos(maxpar,nmaxwalker) 
  Real(8) :: ener_restart
  Real(8) :: ener_atom(maxpar)
  ! global minimum
  !matrices for storing coordinates and energies of the last refused minima
  Real(8) :: ratv(9)
  Real(8) :: xclus(maxpar,maxfail)
  Real(8) :: yclus(maxpar,maxfail)
  Real(8) :: zclus(maxpar,maxfail)
  Real(8) :: enerclus(maxfail),percclus(maxfail)
  ! running coordinates in angstrom
  Real(8) :: xin(maxpar),yin(maxpar),zin(maxpar)
  ! running ccordinates in the minimization routines (in arete units)
  Real(8) :: xat(maxpar),yat(maxpar),zat(maxpar)
  ! final output for each escape minimization
  Real(8) :: xsave(maxpar),ysave(maxpar),zsave(maxpar)
  ! coordinates of the walkers
  Real(8) :: xpos(maxpar,nmaxwalker),ypos(maxpar,nmaxwalker),zpos(maxpar,nmaxwalker)
  ! energies of the walkers
  Real(8) :: enerwalker(nmaxwalker)
  Real(8) :: esubwalker(nmaxwalker)
  ! order parameters ofthe walkers:
  Real(8) :: walker_parameter(nmaxwalker,2),walker_parameter_old(2)
  ! limits of the box for cartesian coordinates
  Real(8) :: x_min,x_max,y_min,y_max,z_min,z_max
  ! this is the variable for calling ranmar
  Real(8) :: zzz
  !temperatures for each kind of move:
  Real(8) :: temp_bonds,temp_ball,temp_shell,temp_shake,temp_highener,temp_single
  Real(8) :: shell_thickness
  Real(8) :: temp_brow,temp_brow_accept,temp_browsurf,temp_browsurf_accept,temp_mirror
  Real(8) :: temp_all,temp_sur,temp_csh,temp_int,temp_sep,temp_mix,temp_mult
  !probability of hyperquenching
  Real(8) :: pquench
  ! temperature for choosing among the refused clusters
  Real(8) :: escapetemp
  !minimum and maximum radii for atom displacement in subroutine move (old type)
  Real(8) :: rho_inf_sk,rho_sup_sk,rho_inf_he,rho_sup_he,rho_inf_sg,rho_sup_sg
  !parameter for mirror move
  Real(8) :: zcut,deltaz
  !parameters for multiple exchanges
  Real(8) :: exp_mult
  !probability to choose different kinds of moves
  Real(8) :: probmove1,probmove2,probmove3,probmove4,probmove5
  Real(8) :: probmove6,probmove7,probmove8,probmove9,probmove10
  !probabilities of the different kinds of exchange moves:
  Real(8) :: prob_exch_all,prob_exch_sur,prob_exch_csh,prob_exch_int, prob_exch_sep,prob_exch_mix
  Real(8) :: prob_exch_mult, prob_exch_gru
  Real(8) :: exp_sep,exp_mix ! parameter for the exchange move SEP and MIX
  Integer  :: bignn_csh,smallnn_csh,bignn_int,smallnn_int ! parameters for the exchange move CSH
  Integer :: npas_bro,npas_brosurf
  Logical :: small_inside
  !Character(3) :: exchange_type
  Integer :: vicini(maxpar,maxpar)
  Integer :: num_vicini(maxpar),num_vicini_1l(maxpar)
  Integer :: quantivicini(maxpar),quantiviAg(maxpar),quantivi_eq(maxpar),quantivi_neq(maxpar)
  Integer :: pv_matrix(maxpar,maxpar)
  Real*8  :: dist_matrix(maxpar,maxpar)
  Real*8 :: param_min,param_max,parameter_tolerance,parameter_interval
  Character*7 :: gruppo
  Character*3 :: caratt
  Real*8  :: nn(3)
  ! number of walkers, kind of walkers turn-over, repulsive strength of walkers, width of repulsive interaction
  Integer :: nwalkers,walkers_turnover,freq_walker_exchange,parameter_partitioning
  Real*8 :: repwalkers
  Real*8 :: widthwalkers(2)
  Real*8 :: legami_specie1,legami_specie2,legami_specie3,legami_totali
  Character(5) :: choose_algorithm
  Real*8 ::g1,g2,g3
  Logical :: lista_minimi
  Logical :: pvett((maxpar*(maxpar+1))/2)
  Logical :: first_minimization

END MODULE BH_MOD
