Subroutine cna
  
  USE PARAMETERS
  USE HISTO
  USE BH_MOD
  
  Implicit None
  
  Integer, Parameter :: nmax_coppie=100000
  Integer :: i,j,n,m,k
  Integer :: nna(nmax_coppie),nnb(nmax_coppie),nnc(nmax_coppie)
  Integer :: vicini_comuni(nmax_coppie,50),num_coppia
  Integer :: fact,vettore_da_perm(nat3d),qprova(nat3d)
  Integer :: lunghezza_max
  Integer :: primo_elemcoppia(nmax_coppie),secondo_elemcoppia(nmax_coppie)
  Integer, dimension(:,:), allocatable :: permut
  Integer, dimension(:), allocatable :: lunghezza_catena
  Logical :: PV_picc(13,13)
  Integer :: somma_sig_ccn,comp=0,a,b
  !Real :: sig555,sig322,sig433,sig321,sig331,sig422,sig421,sig333,sig554,sig666 
  !Real :: sig533,sig211,sig544,sig200,sig311,sig300
  Real :: tot_perc
  Logical :: first
  
  !Per ogni coppia di primi vicini (individuo le
  !coppie di pv esplorando solo la meta' sup della matrice) devo
  !individuare i vicini comuni usando le info contenute in vicini
  PV_picc(:,:)=.false.
  num_coppia=0
!  open(55,file='dettagli_cna.out',status='unknown')
!  open(66,file='signatures.out',status='unknown')
  
!  write(*,*)'comm_neigh'
  do i=1,nat3d-1
     do j=i+1,nat3d
     comp=(i-1)*nat3d-(((i-2)*(i-1))/2)+j-i+1
     !write(*,*)comp
        if (pvett(comp)) then !ho individuato una coppia
           num_coppia=num_coppia+1
	   !write(*,*)num_coppia
	   !write(*,*)'coppia',num_coppia
           !write(55,*)'la coppia',num_coppia,'e` formata da:',i,j
	   primo_elemcoppia(num_coppia)=i
	   secondo_elemcoppia(num_coppia)=j
           !call flush(55)
           nna(num_coppia)=0 !inizializzo il numero di vicini comuni della coppia
           nnb(num_coppia)=0
           !confronto i loro primi vicini per vedere se ne hanno in comune:
           do n=1,num_vicini(i)
              do m=1,num_vicini(j)
                 if(vicini(i,n) .eq. vicini(j,m) )then !un vicino comune
                    nna(num_coppia)=nna(num_coppia)+1 !incremento il numero di vicini comuni
                    vicini_comuni(num_coppia,nna(num_coppia))=vicini(i,n)
                    !conto il numero di legami ai pv esistenti tra i vicini comuni
                    if (nna(num_coppia) .gt. 1)then ! se ho almeno 2 pv comuni
                       do k=1,nna(num_coppia)-1 !cfr tra il vicino nna(num_coppia) e i precedenti
		       a=vicini_comuni(num_coppia,nna(num_coppia)-k)
		       b=vicini_comuni(num_coppia,nna(num_coppia))
		       comp=(a-1)*nat3d-(((a-2)*(a-1))/2)+b-a+1
                          if ( pvett(comp))then
                             nnb(num_coppia)=nnb(num_coppia)+1 !ho trovato un legame ai pv tra i pv comuni
                             !salvo nella matrice legami il legame pv, contandoli 2 volte
                             !legami(num_coppia,(2*nnb(num_coppia))-1)=vicini_comuni(num_coppia,nna(num_coppia))
                             !legami(num_coppia,2*nnb(num_coppia))=vicini_comuni(num_coppia,nna(num_coppia)-k)
                          endif
                       enddo !ho controllato se il nuovo pv comune ha legami ai pv con gli altri
                    endif
                 endif
              enddo
           enddo !ho trovato tutti i vicini comuni
        endif
     enddo
     
  enddo !ho controllato tutte le possibili coppie
!  write (*,*) "num_coppia=", num_coppia 
  do i=1,num_coppia !per ogni coppia devo valutare la lunghezzamax della catena
     !calcolo il fattoriale del numero di vicini comuni della coppia:
     fact=1
     do j=1,nna(i)
        fact=fact*j
     enddo
     !ora conosco le dimensioni della matrice permut e posso allocarla:
     !write(*,*)fact,nna(i)
     allocate(permut(fact,nna(i)))
     allocate(lunghezza_catena(fact))
     !ora riempio permut con le permutazioni delle coponenti del vettore vicini_comuni(i)
     first=.true.
     vettore_da_perm(1:nna(i))=vicini_comuni(i,1:nna(i))
     do m=1,fact
        permut(m,1:nna(i))=vettore_da_perm(1:nna(i))
        lunghezza_catena(m)=0
        !per ogni permutazione dei vicini comuni, calcolo la lunghezza della catena:
        if(nna(i) .eq. 2)then !il caso e' banale: se i vicini comuni sono vicini tra loro, la catena e' lunga 1
	   if(permut(m,2) .lt. permut(m,1))comp=(permut(m,2)-1)*nat3d-(((permut(m,2)-2)*(permut(m,2)-1))/2)+permut(m,1)-permut(m,2)+1
	   if(permut(m,2) .gt. permut(m,1))comp=(permut(m,1)-1)*nat3d-(((permut(m,1)-2)*(permut(m,1)-1))/2)+permut(m,2)-permut(m,1)+1
           if (pvett(comp)) then
              lunghezza_catena(m)=1            
           endif
        elseif (nna(i) .gt. 2)then   
           n=2
           do while (n .le. nna(i))
              if(permut(m,n) .lt. permut(m,n-1))comp=(permut(m,n)-1)*nat3d- &
	    	    (((permut(m,n)-2)*(permut(m,n)-1))/2)+permut(m,n-1)-permut(m,n)+1
	      if(permut(m,n) .gt. permut(m,n-1))comp=(permut(m,n-1)-1)*nat3d- & 
	    	    (((permut(m,n-1)-2)*(permut(m,n-1)-1))/2)+permut(m,n)-permut(m,n-1)+1
              if (pvett(comp))then
                 lunghezza_catena(m)=lunghezza_catena(m)+1
		 if(permut(m,nna(i)) .lt. permut(m,1))comp=(permut(m,nna(i))-1)*nat3d-(((permut(m,nna(i))-2)*(permut(m,nna(i))-1))/2)+&
		 &permut(m,1)-permut(m,nna(i))+1
		 if(permut(m,nna(i)) .gt. permut(m,1))comp=(permut(m,1)-1)*nat3d-(((permut(m,1)-2)*(permut(m,1)-1))/2)+&
		 &permut(m,nna(i))-permut(m,1)+1
                 if((lunghezza_catena(m) .eq. (nna(i)-1)) .and. (pvett(comp))) then
                    lunghezza_catena(m)=lunghezza_catena(m)+1
                 endif
              elseif(.not.pvett(comp))then !appena trovo un elemento che non e' vicino del suo
	      !precedente, esco dal loop e passo ad analizzare la permutazione successiva
                 n=nna(i)
              endif
              n=n+1
           enddo
        
        
        call BestLex(vettore_da_perm(1:nna(i)),nna(i),qprova(1:nna(i)),first)
       
        endif
        
        
     enddo
     
     !trovo il massimo del vettore lunghezza_catena
     lunghezza_max=lunghezza_catena(1)
     do m=2,fact
        if(lunghezza_catena(m) .gt. lunghezza_max) then
           lunghezza_max=lunghezza_catena(m)
        endif
     enddo
     
     nnc(i)=lunghezza_max
     
     deallocate(permut)
     deallocate(lunghezza_catena)
  enddo
  
  !scrivo i dettagli della cna nel file dettagli_cna.out
  
!  do i=1,num_coppia
!     write(55,*)'la coppia',i,'ha',nna(i),'pv comuni'
!     write(55,*)'la coppia',i,'ha',nnb(i),'legami a pv tra i pvcomuni'
!     write(55,*)'la coppia',i,'ha una catena lunga',nnc(i)
!  enddo
  
!  write(55,*)'****************  SIGNATURES (5,5,5)'
!  do i=1,num_coppia
!  if(nna(i)+nnb(i)+nnc(i) .eq. 15)then
!  write(55,*)primo_elemcoppia(i),secondo_elemcoppia(i)
!  endif
!  enddo
  
!  write(55,*)'****************  SIGNATURES (4,2,1)'
!  do i=1,num_coppia
!  if((nna(i)+nnb(i)+nnc(i) .eq. 7) .and. (nna(i) .eq. 4))then
!  write(55,*)primo_elemcoppia(i),secondo_elemcoppia(i)
!  endif
!  enddo
  
!  write(55,*)'****************  SIGNATURES (4,2,2)'
!  do i=1,num_coppia
!  if((nna(i)+nnb(i)+nnc(i) .eq. 8) .and. (nna(i) .eq. 2*nnb(i)))then
!  write(55,*)primo_elemcoppia(i),secondo_elemcoppia(i)
!  endif
!  enddo
  
!  write(55,*)'****************  SIGNATURES (3,1,1)'
!  do i=1,num_coppia
!  if((nna(i) .eq. 3) .and. (nnb(i)*nnc(i) .eq. 1))then
!  write(55,*)primo_elemcoppia(i),secondo_elemcoppia(i)
!  endif
!  enddo
  
  !scrivo l'output che contiene le signature, signatures.out
  sig555=0
  sig544=0
  sig533=0
  sig433=0
  sig421=0
  sig422=0
  sig322=0
  sig321=0
  sig331=0
  sig333=0
  sig554=0
  sig666=0
  sig211=0
  sig200=0
  sig311=0
  sig300=0
  
  do i=1,num_coppia
     if ((nna(i).eq.5) .and. (nnb(i).eq.5) .and. (nnc(i).eq.5)) then
        sig555=sig555+1
     elseif ((nna(i).eq.4) .and. (nnb(i).eq.2) .and. (nnc(i).eq.1)) then
        sig421=sig421+1
     elseif ((nna(i).eq.5) .and. (nnb(i).eq.4) .and. (nnc(i).eq.4)) then
        sig544=sig544+1
     elseif ((nna(i).eq.5) .and. (nnb(i).eq.3) .and. (nnc(i).eq.3)) then
        sig533=sig533+1
     elseif ((nna(i).eq.6) .and. (nnb(i).eq.6) .and. (nnc(i).eq.6)) then
        sig666=sig666+1
     elseif ((nna(i).eq.4) .and. (nnb(i).eq.2) .and. (nnc(i).eq.2)) then
        sig422=sig422+1
     elseif ((nna(i).eq.3) .and. (nnb(i).eq.3) .and. (nnc(i).eq.1)) then
        sig331=sig331+1 
     elseif ((nna(i).eq.3) .and. (nnb(i).eq.3) .and. (nnc(i).eq.3)) then
        sig333=sig333+1
     elseif ((nna(i).eq.3) .and. (nnb(i).eq.2) .and. (nnc(i).eq.1)) then
        sig321=sig321+1
     elseif ((nna(i).eq.3) .and. (nnb(i).eq.2) .and. (nnc(i).eq.2)) then
        sig322=sig322+1
     elseif ((nna(i).eq.3) .and. (nnb(i).eq.1) .and. (nnc(i).eq.1)) then
        sig311=sig311+1
     elseif ((nna(i).eq.4) .and. (nnb(i).eq.3) .and. (nnc(i).eq.3)) then
        sig433=sig433+1
     elseif ((nna(i).eq.5) .and. (nnb(i).eq.5) .and. (nnc(i).eq.4)) then
        sig554=sig554+1
     elseif ((nna(i).eq.2) .and. (nnb(i).eq.1) .and. (nnc(i).eq.1)) then
        sig211=sig211+1
     elseif ((nna(i).eq.2) .and. (nnb(i).eq.0) .and. (nnc(i).eq.0)) then
        sig200=sig200+1
     elseif ((nna(i).eq.3) .and. (nnb(i).eq.0) .and. (nnc(i).eq.0)) then
        sig300=sig300+1
     endif
  enddo
  
  !individuo le strutture chained-common-neighbour, ccn, cioe' quelle che hanno
  !SOLO signatures del tipo (n,n,n) e (n,n-1,n-1)
  
  somma_sig_ccn=sig666+sig555+sig322+sig433+sig333+sig211+sig544
  !write(*,*)'num_coppia',num_coppia
!  if (somma_sig_ccn .eq. num_coppia)then
!     write(*,*)'la struttura E` una CCN, chained-common-neighbours'
!     write(66,*)'la struttura E` una CCN, chained-common-neighbours'
!  else
!     write(*,*)'la struttura NON E` una CCN, chained-common-neighbours'
!     write(66,*)'la struttura NON E` una CCN, chained-common-neighbours'
!  endif
  
  !SCRITTURA OUTPUT
  !normalizzo rispetto al numero di coppie
  sig666=sig666/num_coppia
  sig555=sig555/num_coppia
  sig544=sig544/num_coppia
  sig533=sig533/num_coppia
  sig322=sig322/num_coppia
  sig433=sig433/num_coppia
  sig333=sig333/num_coppia
  sig554=sig554/num_coppia
  sig321=sig321/num_coppia
  sig311=sig311/num_coppia
  sig331=sig331/num_coppia
  sig422=sig422/num_coppia
  sig421=sig421/num_coppia
  sig211=sig211/num_coppia
  sig200=sig200/num_coppia
  sig300=sig300/num_coppia
!      !PER ANALISI SEQUENZE CRESCITA
!      open(67,file='sign.out',status='unknown')
!      write(67,*)'555,421,422'
!      write(67,*)sig555,sig421,sig422
!      close(67)
!  write(66,*)'sig666',sig666
!  write(66,*)'sig555',sig555
!  write(66,*)'sig554',sig554
!  write(66,*)'sig544',sig544
!  write(66,*)'sig533',sig533
!  write(66,*)'sig433',sig433
!  write(66,*)'sig422',sig422
!  write(66,*)'sig421',sig421
!  write(66,*)'sig333',sig333
!  write(66,*)'sig331',sig331
!  write(66,*)'sig322',sig322
!  write(66,*)'sig321',sig321
!  write(66,*)'sig311',sig311
!  write(66,*)'sig300',sig300
!  write(66,*)'sig211',sig211
!  write(66,*)'sig200',sig200
  
  tot_perc=sig666+sig555+sig554+sig544+sig533+sig433+sig422+sig421+sig333+ &
       & sig331+sig322+sig321+sig311+sig300+sig211+sig200
!  write(66,*)'percentuale totale',tot_perc
  
!  close(55)
!  close(66)
  
  return
  
!  !individuo le strutture poli-Ih perfette: ho una struttura perfetta poli-Ih se e'
!  !verificata questa condizione: ciascuno dei suoi atomi devono avere almeno una delle
!  !due seguenti caratteristiche:
!  !1-avere 12 vicini con i quali costituire un perfetto Ih13
!  !2-avere almeno un vicino che realizza la richiesta 1
  
  !inizializzo:
!  atomi_ok_poliIh=0
  
  !prima conto gli atomi che sono al centro di Ih perfetti
!  do i=1,nat3d
!     if (num_vicini(i) .eq. 12) then
!        !write(*,*)'l`atomo',i,'ha 12 vicini'
!        !faccio la cna del clusterino piccolo, costituito da i e dai suoi 12 vicini
!        !devono venire 12 sig 555 e tre sig 322. In questo caso, l'atomo i si trova
!        !al centro di un Ih13.
!        !costruisco l'Ih:
!        Ih(1)=i
!        do k=2,13
!           Ih(k)=vicini(i,k-1)
!        enddo
!        !costruisco la matrice PV_picc dei primi vicini dell'Ih:
!        PV_picc(:,:)=.true. !inizializzo a 1, la diagonale non sara' piu' toccata
!        do k=1,13
!           do j=k+1,13
!	   comp=(k-1)*nat3d-(((k-2)*(k-1))/2)+j-k+1
!              PV_picc(k,j)=pvett(comp)
!              PV_picc(j,k)=pvett(comp)
!           enddo
!        enddo
!        !la passo alla sub che ripete la cna su questo piccolo cluster
!        !      write(*,*)'sto per passare alla cna piccola il clusterino formato da:'
!        !      write(*,*)Ih(1),Ih(2),Ih(3),Ih(4),Ih(5),Ih(6),Ih(7),Ih(8),Ih(9),Ih(10),Ih(11),Ih(12),Ih(13)
!        call cna_piccola(PV_picc,Ih_ok,Ih)
!        if (Ih_ok) then
!           write(*,*)'ho trovato un Ih centrato nell`atomo',i
!           write(66,*)'ho trovato un Ih centrato nell`atomo',i
!           atomi_ok_poliIh=atomi_ok_poliIh+1
!           !scrivo in un vettore gli atomi che sono il centro di un Ih13:
!           atomi_centroIh(atomi_ok_poliIh)=i
!        endif
!     endif
!  enddo
  
  
!  !se l'atomo i non ha 12 vicini, verifico che almeno uno dei suoi vicini
!  !abbia 12 vicini.
!  do i=1,nat3d
!     if (num_vicini(i) .ne. 12) then
!        !controlla se almeno uno dei suoi vicini, k, appartiene alla categoria precedente      
!        controllo=.true.
!        do k=1,num_vicini(i)
!           do j=1,atomi_ok_poliIh
!              if (vicini(i,k) .eq. atomi_centroIh(j)) then
!                 atomi_ok_poliIh=atomi_ok_poliIh+1
!                 controllo=.false.
!                 exit
!              endif
!           enddo
!           if (.not. controllo) exit
!        enddo
!     endif
!  enddo
  
  
!  if (atomi_ok_poliIh .eq. nat3d) then
!     write(*,*)'la struttura E` un perfetto poliIh'
!     write(66,*)'la struttura E` un perfetto poliIh'
!  endif
!  if (atomi_ok_poliIh .ne. nat3d) then
!     write(*,*)'la struttura NON E` un perfetto poliIh'
!     write(66,*)'la struttura NON E` un perfetto poliIh'
!  endif
  
!  !SCRITTURA OUTPUT
!  !normalizzo rispetto al numero di coppie
!  sig666=sig666/num_coppia
!  sig555=sig555/num_coppia
!  sig544=sig544/num_coppia
!  sig533=sig533/num_coppia
!  sig322=sig322/num_coppia
!  sig433=sig433/num_coppia
!  sig333=sig333/num_coppia
!  sig554=sig554/num_coppia
!  sig321=sig321/num_coppia
!  sig311=sig311/num_coppia
!  sig331=sig331/num_coppia
!  sig422=sig422/num_coppia
!  sig421=sig421/num_coppia
!  sig211=sig211/num_coppia
!  sig200=sig200/num_coppia
!  sig300=sig300/num_coppia
!  write(66,*)'sig666',sig666
!  write(66,*)'sig555',sig555
!  write(66,*)'sig554',sig554
!  write(66,*)'sig544',sig544
!  write(66,*)'sig533',sig533
!  write(66,*)'sig433',sig433
!  write(66,*)'sig422',sig422
!  write(66,*)'sig421',sig421
!  write(66,*)'sig333',sig333
!  write(66,*)'sig331',sig331
!  write(66,*)'sig322',sig322
!  write(66,*)'sig321',sig321
!  write(66,*)'sig311',sig311
!  write(66,*)'sig300',sig300
!  write(66,*)'sig211',sig211
!  write(66,*)'sig200',sig200
!  
!  tot_perc=sig666+sig555+sig554+sig544+sig533+sig433+sig422+sig421+sig333+ &
!       & sig331+sig322+sig321+sig311+sig300+sig211+sig200
!  write(66,*)'percentuale totale',tot_perc
!  
!  close(55)
!  close(66)
  
End Subroutine cna


                                                   
