 Subroutine struttura(xcord,ycord,zcord)

  !***********************
  ! This subroutine perform a general analysis of cluster structure.
  !1)nearest neighbors analysis
  !2)distinction between pure and mixed bonds (bimetallic clusters)
  !3)call to cna analysis
  !***********************

  USE PARAMETERS, ONLY: maxpar
  USE FUNCTION_PAR, ONLY: itype,dist,arete !imet1,imet2
  USE BH_MOD, ONLY: vicini,num_vicini,num_vicini_1l,pv_matrix,dist_matrix,&
       gruppo,caratt,nn,nat3d,&
       legami_specie1,legami_specie2,legami_specie3,legami_totali,ind_tot,&
       pvett,ratv
  USE HISTO, ONLY: order_parameter,parameter_value,legamimisti,sig555,sig422,frac_contact

  Implicit None
  Integer :: i,j,comp
  Real*8 :: xcord(maxpar),ycord(maxpar),zcord(maxpar),atomi_bassi,num_medio_vicini_1l,atomi_1l
  Real*8 :: snn(3)

  !riempio la matrice PV:
  pv_matrix(1:nat3d,1:nat3d)=0 !=1 se i,j sono primi vicini, =0 altrimenti
  vicini(1:nat3d,1:nat3d)=0 !vicini(i,j) e' il j-esimo vicino di i
  num_vicini(1:nat3d)=0 !num_vicini(i) e' il numero di vicini di i
  num_vicini_1l(1:nat3d)=0 !num_vicini_1l(i) e' il numero di vicini degli atomi nel 1o layer (z<3.)
  pvett(:)=.false.

    !first neighbours are counted with a 20% tolerance on the nn distance
    snn(1)=nn(1)+0.2d0*nn(1)
    snn(2)=nn(2)+0.2d0*nn(2)
    snn(3)=nn(3)+0.2d0*nn(3)

  !Riempio la matrice delle distanze:
    Do i=1,nat3d
        Do j=i+1,nat3d
    	    dist_matrix(i,j)=sqrt((xcord(i)-xcord(j))*(xcord(i)-xcord(j))+&
                &(ycord(i)-ycord(j))*(ycord(i)-ycord(j))+&
                &(zcord(i)-zcord(j))*(zcord(i)-zcord(j)))
    	    dist_matrix(j,i)=dist_matrix(i,j)
        EndDo !ciclo su j
    EndDo !ciclo su i

   atomi_bassi=0.d0
   do i =1,nat3d
      if (zcord(i).le.3.d0) then
         atomi_bassi=atomi_bassi+1.d0
      endif
   enddo
   frac_contact=atomi_bassi/dfloat(nat3d)


  ! COMINCIO L'ANALISI DEI PRIMI VICINI
  ! individuo i primi vicini
  legami_specie1=0.
  legami_specie2=0.
  legami_specie3=0.
  
  !il seguente loop fa tre cose:
  !1- riempie la matrice pv_matrix(i,j) con 1 se ie j sono primi vicini, 0 altrimenti
  !2- conta il numero di vicini di ciascun atomo e lo scrive nella variabile num_vicini
  !3- riempie la matrice vicini in modo t.c. vicini(i,n)=n-esimo vicino dell'atomo i
     Do i=1,nat3d
        Do j=i+1,nat3d
        comp=(i-1)*nat3d-(((i-2)*(i-1))/2)+j-i+1
           if ((itype(i).eq.1).and.(itype(j).eq.1)) then
              if (dist_matrix(i,j) .lt. snn(1)) then
                 pv_matrix(j,i)=1
                 pv_matrix(i,j)=1
   	      pvett(comp)=.true.
                 legami_specie1=legami_specie1+1
                 num_vicini(i)=num_vicini(i)+1
                 num_vicini(j)=num_vicini(j)+1
                 if(zcord(i) .lt. 3.)then
                   if(zcord(j) .lt. 3)then
                      num_vicini_1l(i)=num_vicini_1l(i)+1
                      num_vicini_1l(j)=num_vicini_1l(j)+1
                   endif
                 endif
                 vicini(i,num_vicini(i))=j
                 vicini(j,num_vicini(j))=i
              endif
           else if ((itype(i).eq.2).and.(itype(j).eq.2)) then
              if (dist_matrix(i,j) .lt. snn(2)) then
                 pv_matrix(j,i)=1
                 pv_matrix(i,j)=1
   	      pvett(comp)=.true.
                 legami_specie2=legami_specie2+1
                 num_vicini(i)=num_vicini(i)+1
                 num_vicini(j)=num_vicini(j)+1
                 if(zcord(i) .lt. 3.)then
                   if(zcord(j) .lt. 3)then
                      num_vicini_1l(i)=num_vicini_1l(i)+1
                      num_vicini_1l(j)=num_vicini_1l(j)+1
                   endif
                 endif
                 vicini(i,num_vicini(i))=j
                 vicini(j,num_vicini(j))=i
              endif
           else
              if (dist_matrix(i,j) .lt. snn(3)) then
                 pv_matrix(j,i)=1
                 pv_matrix(i,j)=1
                 if(zcord(i) .lt. 3.)then
                   if(zcord(j) .lt. 3)then
                      num_vicini_1l(i)=num_vicini_1l(i)+1
                      num_vicini_1l(j)=num_vicini_1l(j)+1
                   endif
                 endif
   	      pvett(comp)=.true.
                 legami_specie3=legami_specie3+1
                 num_vicini(i)=num_vicini(i)+1
                 num_vicini(j)=num_vicini(j)+1
                 vicini(i,num_vicini(i))=j
                 vicini(j,num_vicini(j))=i
              endif
           endif
        EndDo
     EndDo
  legami_totali=legami_specie1+legami_specie2+legami_specie3

  !contoo gli atomi del primo layer che hanno >= 4 vicini nel primo layer, e medio il numero di vicini nel 1l di questa categoria
  atomi_1l=0. !atomi nel primo layer con piu' di 4 vicini nel primo layer
  num_medio_vicini_1l=0. !il loro numero medio di vicini 

  Do i=1,nat3d
     !if((zcord(i) .lt. 3.) .and. (num_vicini_1l(i) .ge. 4))then
     if(zcord(i) .lt. 3)then
       atomi_1l=atomi_1l+1.
       num_medio_vicini_1l=num_medio_vicini_1l+num_vicini_1l(i)
     endif
  EndDo

  if(atomi_1l .ne. 0)then 
     num_medio_vicini_1l=num_medio_vicini_1l/atomi_1l
  else
     num_medio_vicini_1l=0.
  endif


  !Parametro d'ordine
  !1 = legami misti
  !2 = legami puri tipo1
  !3 = legami puri tipo2
  !4 = signature (5,5,5)
  !5 = energy
  !6 = fraction of atoms in the first layer (z<3)
  !7 = average nn number of the first-layer atoms
  if(order_parameter(1) .eq. 1)then
    !il parametro d'ordine "legamimisti" e' la percentuale dei legami misti rispetto al numero totale di legami:
    legamimisti=legami_specie3/legami_totali
    parameter_value(1)=legamimisti
  elseif(order_parameter(1) .eq. 2)then
    parameter_value(1)=legami_specie1/legami_totali
  elseif(order_parameter(1) .eq. 3)then
    parameter_value(1)=legami_specie2/legami_totali
  elseif(order_parameter(1) .eq. 4)then
!    write (*,*) "chiamo cna"
    Call cna !cna assign the parameter value as the value of signature (5,5,5)
!    write (*,*) "esco da cna ", sig555
    parameter_value(1)=sig555
  elseif(order_parameter(1) .eq. 5)then
!    write (*,*) "chiamo cna"
    Call cna !cna assign the parameter value as the value of signature (4,2,2)
!    write (*,*) "esco da cna ", sig422
    parameter_value(1)=sig422
  elseif(order_parameter(1) .eq. 6)then
    !parameter_value(1)=num_medio_vicini_1l
    parameter_value(1)=frac_contact
    !write(*,*)'p_ord:',parameter_value(1)
  elseif(order_parameter(1) .eq. 7)then
    parameter_value(1)=num_medio_vicini_1l
  endif

  if(order_parameter(2) .eq. 1)then
    !il parametro d'ordine "legamimisti" e' la percentuale dei legami misti rispetto al numero totale di legami:
    legamimisti=legami_specie3/legami_totali
    parameter_value(2)=legamimisti
  elseif(order_parameter(2) .eq. 2)then
    parameter_value(2)=legami_specie1/legami_totali
  elseif(order_parameter(2) .eq. 3)then
    parameter_value(2)=legami_specie2/legami_totali
  elseif(order_parameter(2) .eq. 4)then
  Call cna !cna assign the parameter value as the value of signature (5,5,5)
    parameter_value(2)=sig555
  elseif(order_parameter(2) .eq. 5)then
  Call cna !cna assign the parameter value as the value of signature (4,2,2)
    parameter_value(2)=sig422
  elseif(order_parameter(2) .eq. 6)then
    !parameter_value(2)=num_medio_vicini_1l
    parameter_value(2)=frac_contact
  elseif(order_parameter(2) .eq. 7)then
    parameter_value(2)=num_medio_vicini_1l

  endif




  return

End Subroutine struttura
