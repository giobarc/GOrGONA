        subroutine grouping 

        USE PARAMETERS
        USE FUNCTION_PAR
        USE BH_MOD
        USE HISTO

        implicit none 

        integer, Parameter :: nmax_coppie=100000
        integer, dimension(:,:), allocatable:: permut
        integer, dimension(:), allocatable:: lunghezza_catena
        integer,dimension(maxpar):: occ1,occ2
        integer,dimension(maxpar):: coord,coordb,coords,coordbN,coordb1
        integer,dimension(maxpar):: orb_ini,orb_fin
        integer,dimension(maxpar):: orb_ini1,orb_fin1
        integer,dimension(maxpar):: orb_inib,orb_finb
        integer,dimension(maxpar):: string
        integer,dimension(maxpar):: vettore_da_perm,qprova
        integer,dimension(maxpar):: sigb,sigs
        integer:: vicini_comuni(nmax_coppie,50)
        integer::i,ii,j,k,nat,kk,kkk,nmet1,nmet1old,nmet2,nat1,n,m
        integer:: norb,norb1,norbb,min_norbb,min_norbb1
        integer:: lunghezza_max,fact
        integer:: somma_sig_ccn,comp,b,al
        integer:: nbulk,index_bulk
        integer:: n_ini,n_fin,n_step
        integer:: l,ll,lll
        integer:: l1,ll1,lll1
        integer:: l14,ll14,lll14,l15,ll15,lll15,l16,ll16,lll16
        integer:: nsteps,deltamet1,deltamet2
        integer:: n1,n2,n3,n4,n5,n1old,n2old,n3old,n4old,n5old
        integer:: indmin,smax
        integer:: ind_good,index_l
        integer:: indpres
        real(kind=8),dimension(maxpar):: x,y,z,vx,vy,vz,n_v,n_vb,n_vs
        real(kind=8),dimension(maxpar):: xx,yy,zz 
        real(kind=8),dimension(maxpar):: n_vbN,sigbN,pres_atombN
        real(kind=8),dimension(maxpar):: pres_atomb,pres_atoms 
        real(kind=8),dimension(maxpar,3):: c,c1,cb,cs,cbn,cbn1,cbn2,cbn3
        real(kind=8),dimension(maxpar,3):: csave
        real(kind=8),dimension(maxpar,3):: cnew
        real(kind=8),dimension(maxpar,3):: cb1
        real(kind=8),dimension(maxpar,3):: cn,cn1,cn2,cn3
        real(kind=8),dimension(maxpar,3):: css,css1,css2,css3
        real(kind=8),dimension(maxpar):: ene,att,rep,dist_z,distN
        real(kind=8),dimension(maxpar):: distanza 
        real(kind=8),dimension(maxpar):: eneb,distb,enebN,distbN
        real(kind=8),dimension(maxpar):: distsN,enee,enes,dists
        real(kind=8),dimension(maxpar):: distb1,eneb1
        real(kind=8),dimension(3):: snn
        real(kind=8):: xtot,ytot,ztot,tot_perc,cutoff1,n_max
        real(kind=8):: dij,ene_met1,ene_met2,enerMin
        real(kind=8):: r1,A1,csi1,pp1,q1,r2,A2,csi2,pp2,q2
        real(kind=8):: r3,AA3,csi3,pp3,q3
        real(kind=8):: p1,p2,p3,p4,p5,random,enerold,p21,p22,p31,p32,pch
        real(kind=8):: p41,p42,p51,p52,parameterMC,eneratt
        real(kind=8):: prob,rndm1,rndm2,rndm3,rndm4,rndm5
        real(kind=8):: deltae,expdeltae,rndm
        real(kind=8):: alfa,beta,gama
        real(kind=8):: prova1,prova2,enetot
        real(kind=8):: tmp4x2,tmp4y2,tmp4z2
        real(kind=8):: tmp4x3,tmp4y3,tmp4z3
        real(kind=8):: tmp4x4,tmp4y4,tmp4z4
        real(kind=8):: tmp5x2,tmp5y2,tmp5z2
        real(kind=8):: tmp5x3,tmp5y3,tmp5z3
        real(kind=8):: tmp5x4,tmp5y4,tmp5z4
        real(kind=8):: tmp5x5,tmp5y5,tmp5z5
        real(kind=8):: tmp6x2,tmp6y2,tmp6z2
        real(kind=8):: tmp6x3,tmp6y3,tmp6z3
        real(kind=8):: tmp6x4,tmp6y4,tmp6z4
        real(kind=8):: tmp6x5,tmp6y5,tmp6z5
        real(kind=8):: tmp6x6,tmp6y6,tmp6z6
        real(kind=8):: min024,min034,min044
        real(kind=8):: min025,min035,min045,min055
        real(kind=8):: min026,min036,min046,min056,min066
        real(kind=8):: dist024,dist034,dist044
        real(kind=8):: dist025,dist035,dist045,dist055
        real(kind=8):: dist026,dist036,dist046,dist056,dist066
        real(kind=8):: sumQ4,min_sumQ4
        real(kind=8):: sumQ5,min_sumQ5
        real(kind=8):: sumQ6,min_sumQ6
        real(kind=8):: tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp8,tmp9,tmp10
        character(len=2):: tmp7
        real(kind=8):: emin,dmin,pmin
        character(len=2),dimension(maxpar):: chtype
        character(len=2),dimension(maxpar):: specb,specb1,specs,specbN
        character(len=2),dimension(maxpar):: specsave
        character(len=100):: glo_file,res_file
        logical :: PV_picc(13,13)
        logical:: first
        Real(8) :: eneri 
        Integer :: num_coppia 
        Integer :: nna(nmax_coppie),nnb(nmax_coppie),nnc(nmax_coppie)
        Integer :: sig(maxpar)
        Integer :: primo_elemcoppia(nmax_coppie)
        Integer :: secondo_elemcoppia(nmax_coppie)

        !Read input parameters
        open(1,file='input_gru.dat',status='old')
        read(1,*) p1      
        read(1,*) p2      
        read(1,*) p3      
        read(1,*) p4      
        read(1,*) p5      
        read(1,*) parameterMC
        read(1,*) nsteps
        read(1,*) p21,p22
        read(1,*) p31,p32
        read(1,*) p41,p42
        read(1,*) p51,p52
        close(1)

        nat=nat3d
 
!        open (99,file='out_GRU.dat',status='replace')
!        write(res_file,'(A3)')'RES'
! 
!        write(99,*) 'Coordinates have been read...'
!        write(99,*) 'Atoms:',nat
!        write(99,*) 'The following Gouping criteria will be used:'
!        write(99,*) 'CASE 0 Energy and distance from most symm atom'
!        write(99,*) 'CASE 1 Energy and zcoord'
!        write(99,*) 'CASE 2 Energy and distance from z-axis'
!        write(99,*) 'CASE 3 Combination of 1 and 2'
!        write(99,*) 'CASE 4 Like 0 with sign(x)'
!        write(99,*) 'CASE 5 Like 0 with sign(y)'
!        write(99,*) 'CASE 6 Like 0 with sign(z)'
!        write(99,*) 'Specify p1 (threshold for energy):'
!        write(99,*) 'p1 from user:',p1
!        write(99,*) 'Specify p2 (threshold for distance):'
!        write(99,*) 'p2 from user:',p2
!        write(99,*) 'Specify p3 (threshold for z-coord):'
!        write(99,*) 'p3 from user:',p3
!        write(99,*) 'Specify p4 (threshold for distance from z-axis):'
!        write(99,*) 'p4 from user:',p4
!        write(99,*) 'Specify p5 (threshold for pressure):'
!        write(99,*) 'p4 from user:',p5

        do i=1,nat
          c(i,1)=xin(i)
          c(i,2)=yin(i)
          c(i,3)=zin(i)
        enddo

        index_l=0
 110    continue
 
        do i=1,nat
           xin(i)=c(i,1)
           yin(i)=c(i,2)
           zin(i)=c(i,3)
           xat(i)=xin(i)/arete(1)
           yat(i)=yin(i)/arete(1)
           zat(i)=zin(i)/arete(1)
           if (spec(i).eq.nomemet(1)) then !assigning chemical species
              itype(i)=1
           else
              itype(i)=2
           endif
        enddo

        !call force_rgl
        call minimization_rgl(3*nat,0)

        enetot=0.d0

        do i=1,nat
          enetot=enetot+energ(i)
          ene(i)=energ(i)
        enddo

!        write(99,*) 'Starting Energy is', enetot, 'eV'
         
        call force_pres

!        open (3,file='data4',status='replace')
!        do i=1,nat
!         write(3,'(A6,5F10.4,I6,3F10.4,I6,F10.4)') spec(i),c(i,1),c(i,2)&
!     &,c(i,3),ene(i),pres_atom(i)
!        enddo
!        close (3,status='keep')

        ! Calculation of coordination and surface vector
        do i=1,nat
          coord(i)=0
          vx(i)=0.d0
          vy(i)=0.d0
          vz(i)=0.d0
          do j=1,nat
            if (i.ne.j) then
              dij=dsqrt((c(i,1)-c(j,1))**2+(c(i,2)-c(j,2))**2+(c(i,3)-c(&
     &j,3))**2)
              distanza(i)=dsqrt(c(i,1)**2+(c(i,2))**2+(c(i,3))**2) !Distance
              dist_z(i)=dsqrt(c(i,1)**2+(c(i,2))**2)
              if (dij.le.3.d0) then
                coord(i)=coord(i)+1     ! Coordination
                vx(i)=vx(i)+c(i,1)-c(j,1)
                vy(i)=vy(i)+c(i,2)-c(j,2)
                vz(i)=vz(i)+c(i,3)-c(j,3)
              endif
            endif
          enddo
          n_v(i)=vx(i)**2+vy(i)**2+vz(i)**2
        enddo
 
        !riempio la matrice PV:
        pv_matrix(1:nat,1:nat)=0 !=1 se i,j sono primi vicini, =0 altrimenti
        vicini(1:nat,1:nat)=0 !vicini(i,j) e' il j-esimo vicino di i
        num_vicini(1:nat)=0 !num_vicini(i) e' il numero di vicini di i
        pvett(:)=.false.
 
        !first neighbours are counted with a 20% tolerance on the nn
        !distance
        snn(1)=nn(1)+0.2d0*nn(1)
        snn(2)=nn(2)+0.2d0*nn(2)
        snn(3)=nn(3)+0.2d0*nn(3)
 
        !Riempio la matrice delle distanze:
        do i=1,nat
          do j=i+1,nat
            dist_matrix(i,j)=dsqrt((c(i,1)-c(j,1))**2+(c(i,2)-c(j,2))**2&
     &+(c(i,3)-c(j,3))**2)
            dist_matrix(j,i)=dist_matrix(i,j)
          enddo !ciclo su j
        enddo !ciclo su i
 
        ! COMINCIO L'ANALISI DEI PRIMI VICINI
        ! individuo i primi vicini
        legami_specie1=0.d0
        legami_specie2=0.d0
        legami_specie3=0.d0
 
        !il seguente loop fa tre cose:
        !1- riempie la matrice pv_matrix(i,j) con 1 se ie j sono primi vicini, 0 altrimenti
        !2- conta il numero di vicini di ciascun atomo e lo scrive nella variabile num_vicini
        !3- riempie la matrice vicini in modo t.c. vicini(i,n)=n-esimo vicino dell'atomo i
 
        do i=1,nat
          do j=i+1,nat
            comp=(i-1)*nat-(((i-2)*(i-1))/2)+j-i+1
           if ((spec(i).eq.nomemet(1)).and.(spec(j).eq.nomemet(1))) then
              if (dist_matrix(i,j).lt.snn(1)) then
                pv_matrix(j,i)=1
                pv_matrix(i,j)=1
                pvett(comp)=.true.
                legami_specie1=legami_specie1+1
                num_vicini(i)=num_vicini(i)+1
                num_vicini(j)=num_vicini(j)+1
                vicini(i,num_vicini(i))=j
                vicini(j,num_vicini(j))=i
              endif
      else if ((spec(i).eq.nomemet(2)).and.(spec(j).eq.nomemet(2))) then
              if (dist_matrix(i,j).lt.snn(2)) then
                pv_matrix(j,i)=1
                pv_matrix(i,j)=1
                pvett(comp)=.true.
                legami_specie2=legami_specie2+1
                num_vicini(i)=num_vicini(i)+1
                num_vicini(j)=num_vicini(j)+1
                vicini(i,num_vicini(i))=j
                vicini(j,num_vicini(j))=i
              endif
            else
              if (dist_matrix(i,j).lt.snn(3)) then
                pv_matrix(j,i)=1
                pv_matrix(i,j)=1
                pvett(comp)=.true.
                legami_specie3=legami_specie3+1
                num_vicini(i)=num_vicini(i)+1
                num_vicini(j)=num_vicini(j)+1
                vicini(i,num_vicini(i))=j
                vicini(j,num_vicini(j))=i
              endif
            endif
          enddo
        enddo
        legami_totali=legami_specie1+legami_specie2+legami_specie3
 
        num_coppia=0
        PV_picc(:,:)=.false.
 
        do i=1,nat-1
          do j=i+1,nat
            comp=(i-1)*nat-(((i-2)*(i-1))/2)+j-i+1
            if (pvett(comp)) then !ho individuato una coppia
              num_coppia=num_coppia+1
              primo_elemcoppia(num_coppia)=i
              secondo_elemcoppia(num_coppia)=j
              nna(num_coppia)=0 !inizializzo il numero di vicini comuni della coppia
              nnb(num_coppia)=0 !confronto i loro primi vicini per vedere se ne hanno in comune:
              do n=1,num_vicini(i)
                do m=1,num_vicini(j)
                  if(vicini(i,n).eq.vicini(j,m) )then !un vicino comune
                    nna(num_coppia)=nna(num_coppia)+1 !incremento il numero di vicini comuni
                   vicini_comuni(num_coppia,nna(num_coppia))=vicini(i,n)!conto il numero di legami ai pv esistenti tra i vicini comuni
                      if (nna(num_coppia).gt.1)then ! se ho almeno 2 pv comuni
                        do k=1,nna(num_coppia)-1 !cfr tra il vicino nna(num_coppia) e i precedenti
                        al=vicini_comuni(num_coppia,nna(num_coppia)-k)
                        b=vicini_comuni(num_coppia,nna(num_coppia))
                        comp=(al-1)*nat-(((al-2)*(al-1))/2)+b-al+1
                          if (pvett(comp)) then
                            nnb(num_coppia)=nnb(num_coppia)+1 !ho trovato un legame ai pv tra i pv comuni salvo nella matrice legami il legame pv, contandoli 2 volte
                          endif
                       enddo !ho controllato se il nuovo pv comune ha legami ai pv con gli altri
                    endif
                  endif
                enddo
              enddo !ho trovato tutti i vicini comuni
            endif
          enddo
        enddo !ho controllato tutte le possibili coppie
 
        do i=1,num_coppia !per ogni coppia devo valutare la lunghezzamax della catena calcolo il fattoriale del numero di vicini comuni della coppia:
        fact=1
        do j=1,nna(i)
          fact=fact*j
        enddo
        !ora conosco le dimensioni della matrice permut e posso allocarla:
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
            if(permut(m,2).lt.permut(m,1))comp=(permut(m,2)-1)*nat-     &
     &(((permut(m,2)-2)*(permut(m,2)-1))/2)+permut(m,1)-permut(m,2)+1
            if(permut(m,2).gt.permut(m,1))comp=(permut(m,1)-1)*nat-     &
     &(((permut(m,1)-2)*(permut(m,1)-1))/2)+permut(m,2)-permut(m,1)+1
            if (pvett(comp)) then
              lunghezza_catena(m)=1
            endif
          elseif (nna(i) .gt. 2)then
            n=2
            do while (n .le. nna(i))
              if(permut(m,n).lt.permut(m,n-1))comp=(permut(m,n)-1)*nat- &
     &(((permut(m,n)-2)*(permut(m,n)-1))/2)+permut(m,n-1)-permut(m,n)+1
              if(permut(m,n).gt.permut(m,n-1))comp=(permut(m,n-1)-1)*nat&
     &-(((permut(m,n-1)-2)*(permut(m,n-1)-1))/2)+permut(m,n)-permut(m,n-&
     &1)+1
              if (pvett(comp))then
                lunghezza_catena(m)=lunghezza_catena(m)+1
                if(permut(m,nna(i)).lt.permut(m,1))comp=(permut(m,nna(i)&
     &)-1)*nat-(((permut(m,nna(i))-2)*(permut(m,nna(i))-1))/2)+permut(m,&
     &1)-permut(m,nna(i))+1
                if(permut(m,nna(i)).gt.permut(m,1))comp=(permut(m,1)-1)*&
     &nat-(((permut(m,1)-2)*(permut(m,1)-1))/2)+permut(m,nna(i))-permut(&
     &m,1)+1
         if((lunghezza_catena(m).eq.(nna(i)-1)).and.(pvett(comp))) then
                  lunghezza_catena(m)=lunghezza_catena(m)+1
                endif
              elseif(.not.pvett(comp))then !appena trovo un elemento che non e' vicino del suo precedente, esco dal loop e passo ad analizzare la permutazione successiva
                n=nna(i)
              endif
              n=n+1
            enddo
            call BestLex(vettore_da_perm(1:nna(i)),nna(i),qprova(1:nna(i&
     &)),first)
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

        do i=1,nat
          sig(i)=0
        enddo
        do i=1,num_coppia
          if (nna(i).eq.5) then
           sig(primo_elemcoppia(i))=sig(primo_elemcoppia(i))+1
           sig(secondo_elemcoppia(i))=sig(secondo_elemcoppia(i))+1
          endif
          if (nnb(i).eq.5) then
           sig(primo_elemcoppia(i))=sig(primo_elemcoppia(i))+1
           sig(secondo_elemcoppia(i))=sig(secondo_elemcoppia(i))+1
          endif
          if (nnc(i).eq.5) then
           sig(primo_elemcoppia(i))=sig(primo_elemcoppia(i))+1
           sig(secondo_elemcoppia(i))=sig(secondo_elemcoppia(i))+1
          endif
        enddo
 
        ! Ordering the atoms by surf_vector 
        do i=1,nat
          n_max=100.d0
          do j=i,nat
            if (n_v(j).lt.n_max) then
              n_max=n_v(j)
              indmin=j
            endif
          enddo
          tmp1=c(i,1)
          tmp2=c(i,2)
          tmp3=c(i,3)
          tmp4=ene(i)
          tmp5=distanza(i)
          tmp6=coord(i)
          tmp7=spec(i)
          tmp8=n_v(i)
          tmp9=sig(i)
          tmp10=pres_atom(i)
          c(i,1)=c(indmin,1)
          c(i,2)=c(indmin,2)
          c(i,3)=c(indmin,3)
          coord(i)=coord(indmin)
          ene(i)=ene(indmin)
          distanza(i)=distanza(indmin)
          spec(i)=spec(indmin)
          n_v(i)=n_v(indmin)
          sig(i)=sig(indmin)
          pres_atom(i)=pres_atom(indmin)
          c(indmin,1)=tmp1
          c(indmin,2)=tmp2
          c(indmin,3)=tmp3
          ene(indmin)=tmp4
          distanza(indmin)=tmp5
          coord(indmin)=tmp6
          spec(indmin)=tmp7
          n_v(indmin)=tmp8
          sig(indmin)=tmp9
          pres_atom(indmin)=tmp10
        enddo
 
        ! Calculation of the orbits based on coord 
        norb=1
        orb_ini(norb)=1
        do i=2,nat
          if (n_v(i).gt.10) then                    ! This parametr can be tuned
            if (coord(i).lt.10) then
              orb_fin(norb)=i-1
              norb=norb+1
              orb_ini(norb)=i
            endif
          endif
 11       continue
        enddo
        orb_fin(norb)=nat
 
        ! Definition of bulk and surface quantities 
        do i=1,nat
          if (i.le.orb_fin(1)) then
            specb(i)=spec(i)
            cb(i,1)=c(i,1)
            cb(i,2)=c(i,2)
            cb(i,3)=c(i,3)
            coordb(i)=coord(i)
            eneb(i)=ene(i)
            distb(i)=distanza(i)
            n_vb(i)=n_v(i)
            sigb(i)=sig(i)
            pres_atomb(i)=pres_atom(i)
          else
            specs(i)=spec(i)
            cs(i,1)=c(i,1)
            cs(i,2)=c(i,2)
            cs(i,3)=c(i,3)
            coords(i)=coord(i)
            enes(i)=ene(i)
            dists(i)=distanza(i)
            n_vs(i)=n_v(i)
            sigs(i)=sig(i)
            pres_atoms(i)=pres_atom(i)
          endif
        enddo
 
        nbulk=orb_fin(1)
 
!        write(99,*) 'Number of bulk atoms is', nbulk
 
        min_norbb=100
 
        ! Definition of the central atom; it works on bulk atoms 
        do k=1,nbulk
          do j=1,nbulk
            cbN(j,1)=cb(j,1)-c(k,1)
            cbN(j,2)=cb(j,2)-c(k,2)
            cbN(j,3)=cb(j,3)-c(k,3)
            enebN(j)=ene(j)
            coordbN(j)=coord(j)
            distbN(j)=dsqrt((cb(j,1)-cb(k,1))**2+(cb(j,2)-cb(k,2))**2+  &
     &      (cb(j,3)-cb(k,3))**2)
            specbN(j)=specb(j)
            n_vbN(j)=n_vb(j)
            sigbN(j)=sigb(j)
            pres_atombN(j)=pres_atomb(j)
          enddo    
          ! Ordering the atoms by distance 
          do i=1,nbulk
            dmin=100.d0
            do j=i,nbulk
              if (distbN(j).lt.dmin) then
                dmin=distbN(j)
                indmin=j
              endif
            enddo
            tmp1=cbN(i,1)
            tmp2=cbN(i,2)
            tmp3=cbN(i,3)
            tmp4=enebN(i)
            tmp5=distbN(i)
            tmp6=coordbN(i)
            tmp7=specbN(i)
            tmp8=n_vbN(i)
            tmp9=sigbN(i)
            tmp10=pres_atombN(i)
            cbN(i,1)=cbN(indmin,1)
            cbN(i,2)=cbN(indmin,2)
            cbN(i,3)=cbN(indmin,3)
            coordbN(i)=coordbN(indmin)
            enebN(i)=enebN(indmin)
            distbN(i)=distbN(indmin)
            specbN(i)=specbN(indmin)
            n_vbN(i)=n_vbN(indmin)
            sigbN(i)=sigbN(indmin)
            pres_atombN(i)=pres_atombN(indmin)
            cbN(indmin,1)=tmp1
            cbN(indmin,2)=tmp2
            cbN(indmin,3)=tmp3
            enebN(indmin)=tmp4
            distbN(indmin)=tmp5
            coordbN(indmin)=tmp6
            specbN(indmin)=tmp7
            n_vbN(indmin)=tmp8
            sigbN(indmin)=tmp9
            pres_atombN(indmin)=tmp10
          enddo
          norbb=1
          orb_inib(norbb)=1
          do i=2,nbulk
            if ((abs(distbN(i)-distbN(i-1))).gt.p2) then
              orb_finb(norbb)=i-1
              norbb=norbb+1
              orb_inib(norbb)=i
            endif
          enddo
          orb_finb(norbb)=nbulk
          if (norbb.lt.min_norbb) then
            min_norbb=norbb
            index_bulk=k
          endif
        enddo
 
        ! Centering of the structure
        do j=1,nbulk
          cb(j,1)=cb(j,1)-c(index_bulk,1)
          cb(j,2)=cb(j,2)-c(index_bulk,2)
          cb(j,3)=cb(j,3)-c(index_bulk,3)
          distb(j)=dsqrt(cb(j,1)**2+(cb(j,2))**2+(cb(j,3))**2) !Distance
        enddo

!        write(99,*) 'Looking for maximum symmetry on atoms'
!        write(99,*) 'Minimum norbb:',min_norbb, index_bulk
!        write(99,*) 'Looking for maximum symmetry on orientation'
        
        min_sumQ4=20.d0
        min_sumQ5=20.d0
        min_sumQ6=20.d0
 
        ! Finding main orientations
        do l=0,35
          call flush(99)
           !write(99,*) 'Index l is',l
          alfa=l*5*3.1416d0/180.d0
          do ll=0,35
            beta=ll*5*3.1416d0/180.d0
            do lll=0,35
              gama=lll*5*3.1416d0/180.d0
              do j=1,nbulk
                cbn1(j,1)=(+cb(j,1)*cos(alfa)+cb(j,2)*sin(alfa))
                cbn1(j,2)=(-cb(j,1)*sin(alfa)+cb(j,2)*cos(alfa))
                cbn1(j,3)= +cb(j,3)
                cbn2(j,1)= +cbn1(j,1)
                cbn2(j,2)=(+cbn1(j,2)*cos(beta)+cbn1(j,3)*sin(beta))
                cbn2(j,3)=(-cbn1(j,2)*sin(beta)+cbn1(j,3)*cos(beta))
                cbn3(j,1)=(+cbn2(j,3)*cos(gama)+cbn2(j,1)*sin(gama))
                cbn3(j,2)= +cbn2(j,2)
                cbn3(j,3)=(-cbn2(j,3)*sin(gama)+cbn2(j,1)*cos(gama))
              enddo
              do i=1,nbulk
                if (distb(i).lt.3.d0) then
                if ((cbn3(i,1).ge.0.2d0).or.(cbn3(i,2).ge.0.2d0)) then
                  tmp4x2=+cbn3(i,2)
                  tmp4y2=-cbn3(i,1)
                  tmp4z2=+cbn3(i,3)
                  tmp4x3=-cbn3(i,1)
                  tmp4y3=-cbn3(i,2)
                  tmp4z3=+cbn3(i,3)
                  tmp4x4=-cbn3(i,2)
                  tmp4y4=+cbn3(i,1)
                  tmp4z4=+cbn3(i,3)
                  !!!
                  tmp5x2=+cbn3(i,1)*0.309d0+cbn3(i,2)*0.951d0
                  tmp5y2=-cbn3(i,1)*0.951d0+cbn3(i,2)*0.309d0
                  tmp5z2=+cbn3(i,3)
                  tmp5x3=-cbn3(i,1)*0.809d0+cbn3(i,2)*0.588d0
                  tmp5y3=-cbn3(i,1)*0.588d0-cbn3(i,2)*0.809d0
                  tmp5z3=+cbn3(i,3)
                  tmp5x4=-cbn3(i,1)*0.809d0-cbn3(i,2)*0.588d0
                  tmp5y4=+cbn3(i,1)*0.588d0-cbn3(i,2)*0.809d0
                  tmp5z4=+cbn3(i,3)
                  tmp5x5=+cbn3(i,1)*0.309d0-cbn3(i,2)*0.951d0
                  tmp5y5=+cbn3(i,1)*0.951d0+cbn3(i,2)*0.309d0
                  tmp5z5=+cbn3(i,3)
                  !!!
                  tmp6x2=+cbn3(i,1)*0.500d0+cbn3(i,2)*0.866d0
                  tmp6y2=-cbn3(i,1)*0.866d0+cbn3(i,2)*0.500d0
                  tmp6z2=+cbn3(i,3)
                  tmp6x3=-cbn3(i,1)*0.500d0+cbn3(i,2)*0.866d0
                  tmp6y3=-cbn3(i,1)*0.866d0-cbn3(i,2)*0.500d0
                  tmp6z3=+cbn3(i,3)
                  tmp6x4=-cbn3(i,1)
                  tmp6y4=-cbn3(i,2)
                  tmp6z4=+cbn3(i,3)
                  tmp6x5=-cbn3(i,1)*0.500d0-cbn3(i,2)*0.866d0
                  tmp6y5=+cbn3(i,1)*0.866d0-cbn3(i,2)*0.500d0
                  tmp6z5=+cbn3(i,3)
                  tmp6x6=+cbn3(i,1)*0.500d0-cbn3(i,2)*0.866d0
                  tmp6y6=+cbn3(i,1)*0.866d0+cbn3(i,2)*0.500d0
                  tmp6z6=+cbn3(i,3)
                  exit
                endif
                endif
              enddo
              min024=10.d0
              min025=10.d0
              min026=10.d0
              do i=1,nbulk
                dist024=dsqrt((cbn3(i,1)-tmp4x2)**2+(cbn3(i,2)-tmp4y2)  &
     &**2+(cbn3(i,3)-tmp4z2)**2) 
                if (dist024.le.min024) then
                  min024=dist024
                endif
                dist025=dsqrt((cbn3(i,1)-tmp5x2)**2+(cbn3(i,2)-tmp5y2)  &
     &**2+(cbn3(i,3)-tmp5z2)**2) 
                if (dist025.le.min025) then
                  min025=dist025
                endif
                dist026=dsqrt((cbn3(i,1)-tmp6x2)**2+(cbn3(i,2)-tmp6y2)  &
     &**2+(cbn3(i,3)-tmp6z2)**2) 
                if (dist026.le.min026) then
                  min026=dist026
                endif
              enddo
              min034=10.d0
              min035=10.d0
              min036=10.d0
              do i=1,nbulk
                dist034=dsqrt((cbn3(i,1)-tmp4x3)**2+(cbn3(i,2)-tmp4y3)  &
     &**2+(cbn3(i,3)-tmp4z3)**2) 
                if (dist034.le.min034) then
                  min034=dist034
                endif
                dist035=dsqrt((cbn3(i,1)-tmp5x3)**2+(cbn3(i,2)-tmp5y3)  &
     &**2+(cbn3(i,3)-tmp5z3)**2) 
                if (dist035.le.min035) then
                  min035=dist035
                endif
                dist036=dsqrt((cbn3(i,1)-tmp6x3)**2+(cbn3(i,2)-tmp6y3)  &
     &**2+(cbn3(i,3)-tmp6z3)**2) 
                if (dist036.le.min036) then
                  min036=dist036
                endif
              enddo
              min044=10.d0
              min045=10.d0
              min046=10.d0
              do i=1,nbulk
                dist044=dsqrt((cbn3(i,1)-tmp4x4)**2+(cbn3(i,2)-tmp4y4)  &
     &**2+(cbn3(i,3)-tmp4z4)**2) 
                if (dist044.le.min044) then
                  min044=dist044
                endif
                dist045=dsqrt((cbn3(i,1)-tmp5x4)**2+(cbn3(i,2)-tmp5y4)  &
     &**2+(cbn3(i,3)-tmp5z4)**2) 
                if (dist045.le.min045) then
                  min045=dist045
                endif
                dist046=dsqrt((cbn3(i,1)-tmp6x4)**2+(cbn3(i,2)-tmp6y4)  &
     &**2+(cbn3(i,3)-tmp6z4)**2) 
                if (dist046.le.min046) then
                  min046=dist046
                endif
              enddo
              min055=10.d0
              min056=10.d0
              do i=1,nbulk
                dist055=dsqrt((cbn3(i,1)-tmp5x5)**2+(cbn3(i,2)-tmp5y5)  &
     &**2+(cbn3(i,3)-tmp5z5)**2) 
                if (dist055.le.min055) then
                  min055=dist055
                endif
                dist056=dsqrt((cbn3(i,1)-tmp6x5)**2+(cbn3(i,2)-tmp6y5)  &
     &**2+(cbn3(i,3)-tmp6z5)**2) 
                if (dist056.le.min056) then
                  min056=dist056
                endif
              enddo
              min066=10.d0
              do i=1,nbulk
                dist066=dsqrt((cbn3(i,1)-tmp6x6)**2+(cbn3(i,2)-tmp6y6)  &
     &**2+(cbn3(i,3)-tmp6z6)**2) 
                if (dist066.le.min066) then
                  min066=dist066
                endif
              enddo
              sumQ4=0.33d0*(min024+min034+min044)
              sumQ5=0.25d0*(min025+min035+min045+min055)
              sumQ6=0.20d0*(min026+min036+min046+min056+min066)
              if (sumQ4.le.min_sumQ4) then
                min_sumQ4=sumQ4
                l14=l
                ll14=ll
                lll14=lll
              endif                  
              if (sumQ5.le.min_sumQ5) then
                min_sumQ5=sumQ5
                l15=l
                ll15=ll
                lll15=lll
              endif                  
              if (sumQ6.le.min_sumQ6) then
                min_sumQ6=sumQ6
                l16=l
                ll16=ll
                lll16=lll
              endif                  
            enddo  !lll
          enddo  ! ll
        enddo  !  l
!        write(99,*) 'Minimum sum Q4:',min_sumQ4,l14,ll14,lll14
!        write(99,*) 'Minimum sum Q5:',min_sumQ5,l15,ll15,lll15
!        write(99,*) 'Minimum sum Q6:',min_sumQ6,l16,ll16,lll16
        if (min_sumQ4.lt.min_sumQ5) then
          if (min_sumQ4.lt.min_sumQ6) then
            l1=l14
            ll1=ll14
            lll1=lll14
          else
            if (min_sumQ5.lt.min_sumQ6) then 
              l1=l15
              ll1=ll15
              lll1=lll15
            else
              l1=l16
              ll1=ll16
              lll1=lll16
            endif
          endif
        else
          if (min_sumQ5.lt.min_sumQ6) then
            l1=l15
            ll1=ll15
            lll1=lll15
          else
            if (min_sumQ4.lt.min_sumQ6) then 
              l1=l14
              ll1=ll14
              lll1=lll14
            else
              l1=l16
              ll1=ll16
              lll1=lll16
            endif
          endif
        endif 
          
!        write(99,'(A20,3I10)') 'I chose',l1,ll1,lll1
!        write(99,*) 'The structure has been successfully re-centered and&
!     & re-oriented'
        
        ! Ridefinition of the center of mass and distance 
        do j=1,nat
          c1(j,1)=c(j,1)-c(index_bulk,1)
          c1(j,2)=c(j,2)-c(index_bulk,2)
          c1(j,3)=c(j,3)-c(index_bulk,3)
        enddo
 
        do j=1,nat
          c(j,1)=c1(j,1)
          c(j,2)=c1(j,2)
          c(j,3)=c1(j,3)
        enddo

        ! Rotation according to the proper angle 
        alfa=l1*5*3.1416d0/180.d0
        beta=ll1*5*3.1416d0/180.d0
        gama=lll1*5*3.1416d0/180.d0
        do j=1,nat
          cn1(j,1)=(+c(j,1)*cos(alfa)+c(j,2)*sin(alfa))
          cn1(j,2)=(-c(j,1)*sin(alfa)+c(j,2)*cos(alfa))
          cn1(j,3)= +c(j,3)
          cn2(j,1)= +cn1(j,1)
          cn2(j,2)=(+cn1(j,2)*cos(beta)+cn1(j,3)*sin(beta))
          cn2(j,3)=(-cn1(j,2)*sin(beta)+cn1(j,3)*cos(beta))
          c(j,1)=(+cn2(j,3)*cos(gama)+cn2(j,1)*sin(gama))
          c(j,2)= +cn2(j,2)
          c(j,3)=(-cn2(j,3)*sin(gama)+cn2(j,1)*cos(gama))
        enddo
 
!        open (3,file='data3',status='replace')
!        do i=1,nat
!         write(3,'(A6,3F10.4,I6,3F10.4,I6,F10.4)') spec(i),c(i,1),c(i,2)&
!     &,c(i,3),coord(i),ene(i),distanza(i),n_v(i),sig(i),pres_atom(i)
!        enddo
!        close (3,status='keep')
        
        ! Starting Grouping criteria
        if (index_l.eq.0) then
!          write(99,*) 'CASE 0'
          ! 0: Spherical groups,cas1.eq.0
          do j=1,nat
            distanza(j)=dsqrt(c(j,1)**2+c(j,2)**2+c(j,3)**2)
          enddo
          goto 112
        elseif (index_l.eq.1) then
!          write(99,*) 'CASE 1'
          ! 1: Layered,cas1.eq.1
          do j=1,nat
            distanza(j)=c(j,3)
          enddo
          goto 112
        elseif (index_l.eq.2) then
!          write(99,*) 'CASE 2'
          ! 2: Cilinder,cas2.eq.2
          do j=1,nat
            distanza(j)=dsqrt(c(j,1)**2+c(j,2)**2)
          enddo
          goto 112
        elseif (index_l.eq.3) then
!          write(99,*) 'CASE 3'
          ! 3: Combination of 1 and 2
          do j=1,nat
            distanza(j)=anint(c(j,3))*100+dsqrt(c(j,1)**2+c(j,2)**2)
          enddo
          goto 112
        elseif (index_l.eq.4) then
!          write(99,*) 'CASE 4'
          ! 4: Janus y
          do j=1,nat
            distanza(j)=dsqrt(c(j,1)**2+c(j,2)**2+c(j,3)**2)
            if (c(j,2).ge.0.d0) then
              distanza(j)=distanza(j)
            else
              distanza(j)=-distanza(j)
            endif
          enddo
          goto 112
        elseif (index_l.eq.5) then
!          write(99,*) 'CASE 5'
          ! 5: Janus z
          do j=1,nat
            distanza(j)=dsqrt(c(j,1)**2+c(j,2)**2+c(j,3)**2)
            if (c(j,3).ge.0.d0) then
              distanza(j)=distanza(j)
            else
              distanza(j)=-distanza(j)
            endif
          enddo
          goto 112
        elseif (index_l.eq.6) then
!          write(99,*) 'CASE 6'
          ! 6: Janus x
          do j=1,nat
            distanza(j)=dsqrt(c(j,1)**2+c(j,2)**2+c(j,3)**2)
            if (c(j,1).ge.0.d0) then
              distanza(j)=distanza(j)
            else
              distanza(j)=-distanza(j)
            endif
          enddo
          goto 112
        elseif (index_l.eq.7) then
!          write(99,*) 'CASE 7'
          ! 7: Pressure
          ! Ordering on pressure for bulk 
          do i=1,nbulk
            pmin=100.d0
            do j=i,nbulk
              if (pres_atom(j).lt.pmin) then
                pmin=pres_atom(j)
                indmin=j
              endif
            enddo
            tmp1=c(i,1)
            tmp2=c(i,2)
            tmp3=c(i,3)
            tmp4=ene(i)
            tmp5=distanza(i)
            tmp6=coord(i)
            tmp7=spec(i)
            tmp8=n_v(i)
            tmp9=sig(i)
            tmp10=pres_atom(i)
            c(i,1)=c(indmin,1)
            c(i,2)=c(indmin,2)
            c(i,3)=c(indmin,3)
            coord(i)=coord(indmin)
            ene(i)=ene(indmin)
            distanza(i)=distanza(indmin)
            spec(i)=spec(indmin)
            n_v(i)=n_v(indmin)
            sig(i)=sig(indmin)
            pres_atom(i)=pres_atom(indmin)
            c(indmin,1)=tmp1
            c(indmin,2)=tmp2
            c(indmin,3)=tmp3
            ene(indmin)=tmp4
            distanza(indmin)=tmp5
            coord(indmin)=tmp6
            spec(indmin)=tmp7
            n_v(indmin)=tmp8
            sig(indmin)=tmp9
            pres_atom(indmin)=tmp10
          enddo
          goto 1121
        endif
 
 112    continue
 
        ! Ordering on distance for bulk 
        do i=1,nbulk
          dmin=100.d0
          do j=i,nbulk
            if (distanza(j).lt.dmin) then
              dmin=distanza(j)
              indmin=j
            endif
          enddo
          tmp1=c(i,1)
          tmp2=c(i,2)
          tmp3=c(i,3)
          tmp4=ene(i)
          tmp5=distanza(i)
          tmp6=coord(i)
          tmp7=spec(i)
          tmp8=n_v(i)
          tmp9=sig(i)
          c(i,1)=c(indmin,1)
          c(i,2)=c(indmin,2)
          c(i,3)=c(indmin,3)
          coord(i)=coord(indmin)
          ene(i)=ene(indmin)
          distanza(i)=distanza(indmin)
          spec(i)=spec(indmin)
          n_v(i)=n_v(indmin)
          sig(i)=sig(indmin)
          c(indmin,1)=tmp1
          c(indmin,2)=tmp2
          c(indmin,3)=tmp3
          ene(indmin)=tmp4
          distanza(indmin)=tmp5
          coord(indmin)=tmp6
          spec(indmin)=tmp7
          n_v(indmin)=tmp8
          sig(indmin)=tmp9
        enddo

 1121   continue
 
        ! Ordering on energy for surface 
        do i=nbulk+1,nat
          emin=100.d0
          do j=i,nat
            if (ene(j).lt.emin) then
              emin=ene(j)
              indmin=j
            endif
          enddo
          tmp1=c(i,1)
          tmp2=c(i,2)
          tmp3=c(i,3)
          tmp4=ene(i)
          tmp5=distanza(i)
          tmp6=coord(i)
          tmp7=spec(i)
          tmp8=n_v(i)
          tmp9=sig(i)
          c(i,1)=c(indmin,1)
          c(i,2)=c(indmin,2)
          c(i,3)=c(indmin,3)
          coord(i)=coord(indmin)
          ene(i)=ene(indmin)
          distanza(i)=distanza(indmin)
          spec(i)=spec(indmin)
          n_v(i)=n_v(indmin)
          sig(i)=sig(indmin)
          c(indmin,1)=tmp1
          c(indmin,2)=tmp2
          c(indmin,3)=tmp3
          ene(indmin)=tmp4
          distanza(indmin)=tmp5
          coord(indmin)=tmp6
          spec(indmin)=tmp7
          n_v(indmin)=tmp8
          sig(indmin)=tmp9
        enddo
 
        ! Ridefinition of bulk and surface quantities 
        do i=1,nat
          if (i.le.orb_fin(1)) then
            specb(i)=spec(i)
            cb(i,1)=c(i,1)
            cb(i,2)=c(i,2)
            cb(i,3)=c(i,3)
            coordb(i)=coord(i)
            eneb(i)=ene(i)
            distb(i)=distanza(i)
            n_vb(i)=n_v(i)
            sigb(i)=sig(i)
            pres_atomb(i)=pres_atom(i)
          else
            specs(i)=spec(i)
            cs(i,1)=c(i,1)
            cs(i,2)=c(i,2)
            cs(i,3)=c(i,3)
            coords(i)=coord(i)
            enes(i)=ene(i)
            dists(i)=distanza(i)
            n_vs(i)=n_v(i)
            sigs(i)=sig(i)
            pres_atoms(i)=pres_atom(i)
          endif
        enddo
 
!        open (3,file='data1',status='replace')
!        do i=1,nbulk
!          write(3,'(A6,3F10.4,I6,3F10.4,I6,F10.4)') specb(i),cb(i,1),   &
!     &cb(i,2),cb(i,3),coordb(i),eneb(i),distb(i),n_vb(i),sigb(i),       &
!     &pres_atomb(i)
!        enddo
!        close (3,status='keep')
!        open (3,file='data2',status='replace')
!        do i=nbulk+1,nat
!          write(3,'(A6,3F10.4,I6,3F10.4,I6,F10.4)') specs(i),cs(i,1),   &
!     &cs(i,2),cs(i,3),coords(i),enes(i),dists(i),n_vs(i),sigs(i),       &
!     &pres_atoms(i)
!        enddo
!        close (3,status='keep')
        
        norb=1
        orb_ini(norb)=1
        if (index_l.eq.0) then
          pch=p2
        elseif (index_l.eq.1) then
          pch=p3
        elseif (index_l.eq.2) then
          pch=p4
        elseif (index_l.eq.3) then
          pch=p4
        elseif (index_l.eq.4) then
          pch=p2
        elseif (index_l.eq.5) then
          pch=p2
        elseif (index_l.eq.6) then
          pch=p2
        elseif (index_l.eq.7) then
          pch=p5
        endif
        do i=2,nbulk
          if (index_l.eq.7) then
            if ((abs(pres_atom(i)-pres_atom(i-1))).gt.pch) then
              orb_fin(norb)=i-1
              norb=norb+1
              orb_ini(norb)=i
            endif
          else
            if ((abs(distanza(i)-distanza(i-1))).gt.pch) then
              orb_fin(norb)=i-1
              norb=norb+1
              orb_ini(norb)=i
            endif
          endif
        enddo
        orb_fin(norb)=nbulk
        norb=norb+1
        orb_ini(norb)=nbulk+1
        do i=nbulk+2,nat
          if ((ene(i)-ene(i-1)).gt.p1) then             
            orb_fin(norb)=i-1
            norb=norb+1
            orb_ini(norb)=i
          endif
        enddo
        orb_fin(norb)=nat
 
        do i=1,norb
          if (spec(orb_ini(i)).eq.nomemet(1)) then
            do j=orb_ini(i)+1,orb_fin(i)
              if (spec(j).eq.nomemet(2)) then
                do k=j,orb_fin(i)
                  if (spec(k).eq.nomemet(1)) then
                    tmp1=c(k,1)
                    tmp2=c(k,2)
                    tmp3=c(k,3)
                    tmp4=ene(k)
                    tmp5=distanza(k)
                    tmp6=coord(k)
                    tmp7=spec(k)
                    tmp8=n_v(k)
                    tmp9=sig(k)
                    c(k,1)=c(j,1)
                    c(k,2)=c(j,2)
                    c(k,3)=c(j,3)
                    coord(k)=coord(j)
                    ene(k)=ene(j)
                    distanza(k)=distanza(j)
                    spec(k)=spec(j)
                    n_v(k)=n_v(j)
                    sig(k)=sig(j)
                    c(j,1)=tmp1
                    c(j,2)=tmp2
                    c(j,3)=tmp3
                    ene(j)=tmp4
                    distanza(j)=tmp5
                    coord(j)=tmp6
                    spec(j)=tmp7
                    n_v(j)=tmp8
                    sig(j)=tmp9
                  endif
                enddo
              endif
            enddo
          else 
            do j=orb_ini(i)+1,orb_fin(i)
              if (spec(j).eq.nomemet(1)) then
                do k=j,orb_fin(i)
                  if (spec(k).eq.nomemet(2)) then
                    tmp1=c(k,1)
                    tmp2=c(k,2)
                    tmp3=c(k,3)
                    tmp4=ene(k)
                    tmp5=distanza(k)
                    tmp6=coord(k)
                    tmp7=spec(k)
                    tmp8=n_v(k)
                    tmp9=sig(k)
                    c(k,1)=c(j,1)
                    c(k,2)=c(j,2)
                    c(k,3)=c(j,3)
                    coord(k)=coord(j)
                    ene(k)=ene(j)
                    distanza(k)=distanza(j)
                    spec(k)=spec(j)
                    n_v(k)=n_v(j)
                    sig(k)=sig(j)
                    c(j,1)=tmp1
                    c(j,2)=tmp2
                    c(j,3)=tmp3
                    ene(j)=tmp4
                    distanza(j)=tmp5
                    coord(j)=tmp6
                    spec(j)=tmp7
                    n_v(j)=tmp8
                    sig(j)=tmp9
                  endif
                enddo
              endif
            enddo 
          endif
        enddo
 
        norb1=norb
 
        ii=1
        do i=1,norb
          orb_ini1(ii)=orb_ini(i)
          if ((orb_fin(i)-orb_ini(i)+1).eq.1) then
            orb_fin1(ii)=orb_fin(i)
            ii=ii+1
            goto 14
          endif
          do j=orb_ini(i)+1,orb_fin(i)
            if (spec(j).ne.(spec(j-1))) then
              orb_fin1(ii)=j-1
              norb1=norb1+1
              ii=ii+1
              orb_ini1(ii)=j
              orb_fin1(ii)=orb_fin(i)
              ii=ii+1
              goto 14
            endif
          enddo
          orb_fin1(ii)=orb_fin(i)
          ii=ii+1
  14      continue
        enddo
        
        norb=norb1
!        if (norb.ge.200) goto 113
        do i=1,norb1
          orb_ini(i)=orb_ini1(i)      
          orb_fin(i)=orb_fin1(i)      
        enddo
 
!        write(99,*) 'Final Orbits:',norb
! 
!        do i=1,norb
!          write(99,'(A6,3I6,A6,I6)') 'Orbit',i,orb_ini(i),orb_fin(i),   &
!          &'ndim=',orb_fin(i)-orb_ini(i)+1
!        enddo
! 
!        open (3,file='data',status='replace')
!        do i=1,nat
!           write(3,'(A6,3F10.4,I6,3F10.4,I6,F10.4)') spec(i),c(i,1),    &
!     &c(i,2),c(i,3),coord(i),ene(i),distanza(i),n_v(i),sig(i),          &
!     &pres_atom(i)
!        enddo
!        close (3,status='keep')
!        
!        open (1,file='orb.dat',status='replace')
!        write(1,'(I6)') norb
        do i=1,norb
          occ1(i)=0
          occ2(i)=0
          do j=orb_ini(i),orb_fin(i)
            if (spec(j).eq.nomemet(1)) then
              occ1(i)=occ1(i)+1
            else
              occ2(i)=occ2(i)+1
            endif
          enddo
!          write(1,'(3I6)') orb_ini(i),orb_fin(i),orb_fin(i)-orb_ini(i)+1
        enddo
!        close (1,status='keep')
        
        nmet1=0
        nmet2=0
        do i=1,norb
          if (occ1(i).ge.occ2(i)) then
            string(i)=0
            nmet1=nmet1+orb_fin(i)-orb_ini(i)+1
          else 
            string(i)=1
            nmet2=nmet2+orb_fin(i)-orb_ini(i)+1
          endif
        enddo
!        write(99,'(500I1)') (string(i),i=1,norb)
!        write(res_file,'(A3,I1)')'RES',index_l
!        open (2,file=res_file,status='unknown')

        call minimization_rgl(3*nat,0)

!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,0,2,'Energy of starting config',(string(j),j=1,norb)
 
        enerold=ener
        nmet1old=nmet1
 
!        open (5,file='paramMC.dat',status='old')
!        read(5,*) parameterMC
!        read(5,*) nsteps
!        read(5,*) p21,p22
!        read(5,*) p31,p32
!        read(5,*) p41,p42
!        read(5,*) p51,p52
!        close (5,status='keep')
 
        ind_good=0
         
        !write(99,*) 'BH steps=',nsteps

        do k=1,nsteps  ! Loop BH
 
          CALL RANDOM_NUMBER(random)
          prob=nint(random*100.d0)
 
          if ((prob.ge.p21).and.(prob.le.p22)) then  ! Two moves 
 
            do ii=1,10000  ! Until it finds a proper combination
            
              CALL RANDOM_NUMBER(rndm1)
              n1=nint(rndm1*norb)
              if (n1.eq.0) n1=1
              CALL RANDOM_NUMBER(rndm2)
              n2=nint(rndm2*norb)
              if (n2.eq.0) n2=1
            
              deltamet2=0
              deltamet1=0
            
              if (n1.ne.n2) then
                n1old=string(n1)
                n2old=string(n2)
                do i=1,nat
                  specsave(i)=spec(i)
                  csave(i,1)=c(i,1)
                  csave(i,2)=c(i,2)
                  csave(i,3)=c(i,3)
                enddo
                if (string(n1).eq.0) then
                  string(n1)=1
                  deltamet1=deltamet1+orb_fin(n1)-orb_ini(n1)+1
                else
                  string(n1)=0
                  deltamet2=deltamet2+orb_fin(n1)-orb_ini(n1)+1
                endif
                if (string(n2).eq.0) then
                  string(n2)=1
                  deltamet1=deltamet1+orb_fin(n2)-orb_ini(n2)+1
                else
                  string(n2)=0
                  deltamet2=deltamet2+orb_fin(n2)-orb_ini(n2)+1
                endif
                if (deltamet1.eq.deltamet2) then
                  nmet1=0
                  nmet2=0
                  do kk=1,norb
                    do kkk=orb_ini(kk),orb_fin(kk)
                      if (string(kk).eq.0) then
                        spec(kkk)=nomemet(1)
                        nmet1=nmet1+1
                      else
                        spec(kkk)=nomemet(2)
                        nmet2=nmet2+1
                      endif
                    enddo
                  enddo
                  do i=1,nat
                     xin(i)=c(i,1)
                     yin(i)=c(i,2)
                     zin(i)=c(i,3)
                     xat(i)=xin(i)/arete(1)
                     yat(i)=yin(i)/arete(1)
                     zat(i)=zin(i)/arete(1)
                     if (spec(i).eq.nomemet(1)) then !assigning chemical species
                        itype(i)=1
                     else
                        itype(i)=2
                     endif
                  enddo
                  call minimization_rgl(3*nat,0)
                  if (ener.le.enerold) then
                    enerold=ener
                    !eneratt=ener
                    do i=1,nat
                       xin(i)=xat(i)*arete(1)
                       yin(i)=yat(i)*arete(1)
                       zin(i)=zat(i)*arete(1)
                       c(i,1)=xin(i)
                       c(i,2)=yin(i)
                       c(i,3)=zin(i)
                    enddo
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,2,'Energy is  LOWER, accepted',(string(j),j=1,norb)
!                    if (ener.le.enerMin) then
!                      ind_good=ind_good+1
!                     !glo###.xyz file name:
!                     if (index_l.eq.0) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG0',ind_good,'.xyz'
!                     elseif (index_l.eq.1) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG1',ind_good,'.xyz'
!                     elseif (index_l.eq.2) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG2',ind_good,'.xyz'
!                     elseif (index_l.eq.3) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG3',ind_good,'.xyz'
!                     elseif (index_l.eq.4) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG4',ind_good,'.xyz'
!                     elseif (index_l.eq.5) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG5',ind_good,'.xyz'
!                     elseif (index_l.eq.6) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG6',ind_good,'.xyz'
!                     endif
!                     open(1000,file=glo_file,status='unknown')
!                     write(1000,*) nat
!              write(1000,'(a2,1x,a2,f16.8)') nomemet(1),nomemet(2),ener
!                     do i=1,nat
!                   write(1000,'(a2,3f16.8)')spec(i),xin(i),yin(i),zin(i)
!                     enddo
!                     call flush(1000)
!                     close(1000)
!                    endif
                    goto 22
                  else
                    deltae=ener-enerold
                    deltae=deltae/(abs(deltamet2))
                    deltae=deltae/parameterMC
                    expdeltae=dexp(-deltae)
                    CALL RANDOM_NUMBER(rndm)
                    if (expdeltae.gt.rndm) then
                      enerold=ener
                      !eneratt=ener
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,2,'Energy is HIGHER, accepted',(string(j),j=1,norb)
                      do i=1,nat
                         xin(i)=xat(i)*arete(1)
                         yin(i)=yat(i)*arete(1)
                         zin(i)=zat(i)*arete(1)
                         c(i,1)=xin(i)
                         c(i,2)=yin(i)
                         c(i,3)=zin(i)
                      enddo
                      goto 22
                    else
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,2,'Energy is HIGHER,  refused',(string(j),j=1,norb)
                      ener=enerold
                      string(n1)=n1old
                      string(n2)=n2old
                      do i=1,nat
                        spec(i)=specsave(i)
                        c(i,1)=csave(i,1)
                        c(i,2)=csave(i,2)
                        c(i,3)=csave(i,3)
                      enddo
                      goto 22
                    endif
                  endif
                else
                  string(n1)=n1old
                  string(n2)=n2old
                endif
              endif
            enddo ! Until it finds a proper combination
 22         continue
             
          elseif ((prob.ge.p31).and.(prob.le.p32)) then   ! Three moves

            do ii=1,10000  ! Until it finds a proper combination
 
              CALL RANDOM_NUMBER(rndm1)
              n1=nint(rndm1*norb)
              if (n1.eq.0) n1=1
              CALL RANDOM_NUMBER(rndm2)
              n2=nint(rndm2*norb)
              if (n2.eq.0) n2=1
              CALL RANDOM_NUMBER(rndm3)
              n3=nint(rndm3*norb)
              if (n3.eq.0) n3=1
              
              deltamet2=0
              deltamet1=0
              
              if (n1.ne.n2) then
              if (n1.ne.n3) then
              if (n2.ne.n3) then
                n1old=string(n1)
                n2old=string(n2)
                n3old=string(n3)
                do i=1,nat
                  specsave(i)=spec(i)
                  csave(i,1)=c(i,1)
                  csave(i,2)=c(i,2)
                  csave(i,3)=c(i,3)
                enddo
                if (string(n1).eq.0) then
                  string(n1)=1
                  deltamet1=deltamet1+orb_fin(n1)-orb_ini(n1)+1
                else
                  string(n1)=0
                  deltamet2=deltamet2+orb_fin(n1)-orb_ini(n1)+1
                endif
                if (string(n2).eq.0) then
                  string(n2)=1
                  deltamet1=deltamet1+orb_fin(n2)-orb_ini(n2)+1
                else
                  string(n2)=0
                  deltamet2=deltamet2+orb_fin(n2)-orb_ini(n2)+1
                endif
                if (string(n3).eq.0) then
                   string(n3)=1
                   deltamet1=deltamet1+orb_fin(n3)-orb_ini(n3)+1
                else
                   string(n3)=0
                   deltamet2=deltamet2+orb_fin(n3)-orb_ini(n3)+1
                endif
                if (deltamet1.eq.deltamet2) then
                  nmet1=0
                  nmet2=0
                  do kk=1,norb
                    do kkk=orb_ini(kk),orb_fin(kk)
                      if (string(kk).eq.0) then
                        spec(kkk)=nomemet(1)
                        nmet1=nmet1+1
                      else
                        spec(kkk)=nomemet(2)
                        nmet2=nmet2+1
                      endif
                    enddo
                  enddo
                  do i=1,nat
                     xin(i)=c(i,1)
                     yin(i)=c(i,2)
                     zin(i)=c(i,3)
                     xat(i)=xin(i)/arete(1)
                     yat(i)=yin(i)/arete(1)
                     zat(i)=zin(i)/arete(1)
                     if (spec(i).eq.nomemet(1)) then !assigning chemical species
                        itype(i)=1
                     else
                        itype(i)=2
                     endif
                  enddo
                  call minimization_rgl(3*nat,0)
                  if (ener.le.enerold) then
                    enerold=ener
                    !eneratt=ener
                    do i=1,nat
                       xin(i)=xat(i)*arete(1)
                       yin(i)=yat(i)*arete(1)
                       zin(i)=zat(i)*arete(1)
                       c(i,1)=xin(i)
                       c(i,2)=yin(i)
                       c(i,3)=zin(i)
                    enddo
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,3,'Energy is  LOWER, accepted',(string(j),j=1,norb)
!                    if (ener.le.enerMin) then
!                      ind_good=ind_good+1
!                     !glo###.xyz file name:
!                     if (index_l.eq.0) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG0',ind_good,'.xyz'
!                     elseif (index_l.eq.1) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG1',ind_good,'.xyz'
!                     elseif (index_l.eq.2) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG2',ind_good,'.xyz'
!                     elseif (index_l.eq.3) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG3',ind_good,'.xyz'
!                     elseif (index_l.eq.4) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG4',ind_good,'.xyz'
!                     elseif (index_l.eq.5) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG5',ind_good,'.xyz'
!                     elseif (index_l.eq.6) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG6',ind_good,'.xyz'
!                     endif
!                     open(1000,file=glo_file,status='unknown')
!                     write(1000,*) nat
!              write(1000,'(a2,1x,a2,f16.8)') nomemet(1),nomemet(2),ener
!                     do i=1,nat
!                   write(1000,'(a2,3f16.8)')spec(i),xin(i),yin(i),zin(i)
!                     enddo
!                     call flush(1000)
!                     close(1000)
!                    endif
                    goto 33 
                  else
                    deltae=ener-enerold
                    deltae=deltae/(abs(deltamet2))
                    deltae=deltae/parameterMC
                    expdeltae=dexp(-deltae)
                    CALL RANDOM_NUMBER(rndm)
                    if (expdeltae.gt.rndm) then
                      enerold=ener
                      !eneratt=ener
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,3,'Energy is HIGHER, accepted',(string(j),j=1,norb)
                      do i=1,nat
                         xin(i)=xat(i)*arete(1)
                         yin(i)=yat(i)*arete(1)
                         zin(i)=zat(i)*arete(1)
                         c(i,1)=xin(i)
                         c(i,2)=yin(i)
                         c(i,3)=zin(i)
                      enddo
                      goto 33 
                    else
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,3,'Energy is HIGHER,  refused',(string(j),j=1,norb)
                      ener=enerold
                      string(n1)=n1old
                      string(n2)=n2old
                      string(n3)=n3old
                      do i=1,nat
                        spec(i)=specsave(i)
                        c(i,1)=csave(i,1)
                        c(i,2)=csave(i,2)
                        c(i,3)=csave(i,3)
                      enddo
                      goto 33 
                    endif
                  endif
                else
                  string(n1)=n1old
                  string(n2)=n2old
                  string(n3)=n3old
                endif
              endif
              endif
              endif
            enddo ! Until it finds a proper combination
 33         continue
 
          elseif ((prob.ge.p41).and.(prob.le.p42)) then   ! Four moves
 
            do ii=1,10000  ! Until it finds a proper combination
 
              CALL RANDOM_NUMBER(rndm1)
              n1=nint(rndm1*norb)
              if (n1.eq.0) n1=1
              CALL RANDOM_NUMBER(rndm2)
              n2=nint(rndm2*norb)
              if (n2.eq.0) n2=1
              CALL RANDOM_NUMBER(rndm3)
              n3=nint(rndm3*norb)
              if (n3.eq.0) n3=1
              CALL RANDOM_NUMBER(rndm4)
              n4=nint(rndm4*norb)
              if (n4.eq.0) n4=1
 
              deltamet2=0
              deltamet1=0
              
              if (n1.ne.n2) then
              if (n1.ne.n3) then
              if (n1.ne.n4) then
              if (n2.ne.n3) then
              if (n2.ne.n4) then
              if (n3.ne.n4) then
                n1old=string(n1)
                n2old=string(n2)
                n3old=string(n3)
                n4old=string(n4)
                do i=1,nat
                  specsave(i)=spec(i)
                  csave(i,1)=c(i,1)
                  csave(i,2)=c(i,2)
                  csave(i,3)=c(i,3)
                enddo
                if (string(n1).eq.0) then
                   string(n1)=1
                   deltamet1=deltamet1+orb_fin(n1)-orb_ini(n1)+1
                else
                   string(n1)=0
                   deltamet2=deltamet2+orb_fin(n1)-orb_ini(n1)+1
                endif
                if (string(n2).eq.0) then
                   string(n2)=1
                   deltamet1=deltamet1+orb_fin(n2)-orb_ini(n2)+1
                else
                   string(n2)=0
                   deltamet2=deltamet2+orb_fin(n2)-orb_ini(n2)+1
                endif
                if (string(n3).eq.0) then
                   string(n3)=1
                   deltamet1=deltamet1+orb_fin(n3)-orb_ini(n3)+1
                else
                   string(n3)=0
                   deltamet2=deltamet2+orb_fin(n3)-orb_ini(n3)+1
                endif
                if (string(n4).eq.0) then
                   string(n4)=1
                   deltamet1=deltamet1+orb_fin(n4)-orb_ini(n4)+1
                else
                   string(n4)=0
                   deltamet2=deltamet2+orb_fin(n4)-orb_ini(n4)+1
                endif
                if (deltamet1.eq.deltamet2) then
                  nmet1=0
                  nmet2=0
                  do kk=1,norb
                    do kkk=orb_ini(kk),orb_fin(kk)
                      if (string(kk).eq.0) then
                        spec(kkk)=nomemet(1)
                        nmet1=nmet1+1
                      else
                        spec(kkk)=nomemet(2)
                        nmet2=nmet2+1
                      endif
                    enddo
                  enddo
                  do i=1,nat
                     xin(i)=c(i,1)
                     yin(i)=c(i,2)
                     zin(i)=c(i,3)
                     xat(i)=xin(i)/arete(1)
                     yat(i)=yin(i)/arete(1)
                     zat(i)=zin(i)/arete(1)
                     if (spec(i).eq.nomemet(1)) then !assigning chemical species
                        itype(i)=1
                     else
                        itype(i)=2
                     endif
                  enddo
                  call minimization_rgl(3*nat,0)
                  if (ener.le.enerold) then
                    enerold=ener
                    !eneratt=ener
                    do i=1,nat
                       xin(i)=xat(i)*arete(1)
                       yin(i)=yat(i)*arete(1)
                       zin(i)=zat(i)*arete(1)
                       c(i,1)=xin(i)
                       c(i,2)=yin(i)
                       c(i,3)=zin(i)
                    enddo
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,4,'Energy is  LOWER, accepted',(string(j),j=1,norb)
!                    if (ener.le.enerMin) then
!                      ind_good=ind_good+1
!                     !glo###.xyz file name:
!                     if (index_l.eq.0) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG0',ind_good,'.xyz'
!                     elseif (index_l.eq.1) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG1',ind_good,'.xyz'
!                     elseif (index_l.eq.2) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG2',ind_good,'.xyz'
!                     elseif (index_l.eq.3) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG3',ind_good,'.xyz'
!                     elseif (index_l.eq.4) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG4',ind_good,'.xyz'
!                     elseif (index_l.eq.5) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG5',ind_good,'.xyz'
!                     elseif (index_l.eq.6) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG6',ind_good,'.xyz'
!                     endif
!                     open(1000,file=glo_file,status='unknown')
!                     write(1000,*) nat
!              write(1000,'(a2,1x,a2,f16.8)') nomemet(1),nomemet(2),ener
!                     do i=1,nat
!                   write(1000,'(a2,3f16.8)')spec(i),xin(i),yin(i),zin(i)
!                     enddo
!                     call flush(1000)
!                     close(1000)
!                    endif
                    goto 44 
                  else
                    deltae=ener-enerold
                    deltae=deltae/(abs(deltamet2))
                    deltae=deltae/parameterMC
                    expdeltae=dexp(-deltae)
                    CALL RANDOM_NUMBER(rndm)
                    if (expdeltae.gt.rndm) then
                      enerold=ener
                      !eneratt=ener
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,          &
!     &nmet2,k,4,'Energy is HIGHER, accepted',(string(j),j=1,norb)
                      do i=1,nat
                         xin(i)=xat(i)*arete(1)
                         yin(i)=yat(i)*arete(1)
                         zin(i)=zat(i)*arete(1)
                         c(i,1)=xin(i)
                         c(i,2)=yin(i)
                         c(i,3)=zin(i)
                      enddo
                      goto 44 
                    else
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,4,'Energy is HIGHER,  refused',(string(j),j=1,norb)
                      ener=enerold
                      string(n1)=n1old
                      string(n2)=n2old
                      string(n3)=n3old
                      string(n4)=n4old
                      do i=1,nat
                        spec(i)=specsave(i)
                        c(i,1)=csave(i,1)
                        c(i,2)=csave(i,2)
                        c(i,3)=csave(i,3)
                      enddo
                      goto 44 
                    endif
                  endif
                else
                  string(n1)=n1old
                  string(n2)=n2old
                  string(n3)=n3old
                  string(n4)=n4old
                endif
              endif
              endif
              endif
              endif
              endif
              endif
            enddo ! Until it finds a proper combination
 44         continue
 
          elseif ((prob.ge.p51).and.(prob.le.p52)) then   ! Five moves
 
            do ii=1,10000  ! Until it finds a proper combination
 
              CALL RANDOM_NUMBER(rndm1)
              n1=nint(rndm1*norb)
              if (n1.eq.0) n1=1
              CALL RANDOM_NUMBER(rndm2)
              n2=nint(rndm2*norb)
              if (n2.eq.0) n2=1
              CALL RANDOM_NUMBER(rndm3)
              n3=nint(rndm3*norb)
              if (n3.eq.0) n3=1
              CALL RANDOM_NUMBER(rndm4)
              n4=nint(rndm4*norb)
              if (n4.eq.0) n4=1
              CALL RANDOM_NUMBER(rndm5)
              n5=nint(rndm5*norb)
              if (n5.eq.0) n5=1
              
              deltamet2=0
              deltamet1=0
              
              if (n1.ne.n2) then
              if (n1.ne.n3) then
              if (n1.ne.n4) then
              if (n1.ne.n5) then
              if (n2.ne.n3) then
              if (n2.ne.n4) then
              if (n2.ne.n5) then
              if (n3.ne.n4) then
              if (n3.ne.n5) then
              if (n4.ne.n5) then
                n1old=string(n1)
                n2old=string(n2)
                n3old=string(n3)
                n4old=string(n4)
                n5old=string(n5)
                do i=1,nat
                  specsave(i)=spec(i)
                  csave(i,1)=c(i,1)
                  csave(i,2)=c(i,2)
                  csave(i,3)=c(i,3)
                enddo
                if (string(n1).eq.0) then
                   string(n1)=1
                   deltamet1=deltamet1+orb_fin(n1)-orb_ini(n1)+1
                else
                   string(n1)=0
                   deltamet2=deltamet2+orb_fin(n1)-orb_ini(n1)+1
                endif
                if (string(n2).eq.0) then
                   string(n2)=1
                   deltamet1=deltamet1+orb_fin(n2)-orb_ini(n2)+1
                else
                   string(n2)=0
                   deltamet2=deltamet2+orb_fin(n2)-orb_ini(n2)+1
                endif
                if (string(n3).eq.0) then
                   string(n3)=1
                   deltamet1=deltamet1+orb_fin(n3)-orb_ini(n3)+1
                else
                   string(n3)=0
                   deltamet2=deltamet2+orb_fin(n3)-orb_ini(n3)+1
                endif
                if (string(n4).eq.0) then
                   string(n4)=1
                   deltamet1=deltamet1+orb_fin(n4)-orb_ini(n4)+1
                else
                   string(n4)=0
                   deltamet2=deltamet2+orb_fin(n4)-orb_ini(n4)+1
                endif
                if (string(n5).eq.0) then
                   string(n5)=1
                   deltamet1=deltamet1+orb_fin(n5)-orb_ini(n5)+1
                else
                   string(n5)=0
                   deltamet2=deltamet2+orb_fin(n5)-orb_ini(n5)+1
                endif
                if (deltamet1.eq.deltamet2) then
                  nmet1=0
                  nmet2=0
                  do kk=1,norb
                    do kkk=orb_ini(kk),orb_fin(kk)
                      if (string(kk).eq.0) then
                        spec(kkk)=nomemet(1)
                        nmet1=nmet1+1
                      else
                        spec(kkk)=nomemet(2)
                        nmet2=nmet2+1
                      endif
                    enddo
                  enddo
                  do i=1,nat
                     xin(i)=c(i,1)
                     yin(i)=c(i,2)
                     zin(i)=c(i,3)
                     xat(i)=xin(i)/arete(1)
                     yat(i)=yin(i)/arete(1)
                     zat(i)=zin(i)/arete(1)
                     if (spec(i).eq.nomemet(1)) then !assigning chemical species
                        itype(i)=1
                     else
                        itype(i)=2
                     endif
                  enddo
                  call minimization_rgl(3*nat,0)
                  if (ener.le.enerold) then
                    enerold=ener
                    !eneratt=ener
                    do i=1,nat
                       xin(i)=xat(i)*arete(1)
                       yin(i)=yat(i)*arete(1)
                       zin(i)=zat(i)*arete(1)
                       c(i,1)=xin(i)
                       c(i,2)=yin(i)
                       c(i,3)=zin(i)
                    enddo
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,5,'Energy is  LOWER, accepted',(string(j),j=1,norb)
!                    if (ener.le.enerMin) then
!                      ind_good=ind_good+1
!                     !glo###.xyz file name:
!                     if (index_l.eq.0) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG0',ind_good,'.xyz'
!                     elseif (index_l.eq.1) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG1',ind_good,'.xyz'
!                     elseif (index_l.eq.2) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG2',ind_good,'.xyz'
!                     elseif (index_l.eq.3) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG3',ind_good,'.xyz'
!                     elseif (index_l.eq.4) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG4',ind_good,'.xyz'
!                     elseif (index_l.eq.5) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG5',ind_good,'.xyz'
!                     elseif (index_l.eq.6) then
!                   write(glo_file,'(A4,I3.3,A4)')'glG6',ind_good,'.xyz'
!                     endif
!                     open(1000,file=glo_file,status='unknown')
!                     write(1000,*) nat
!              write(1000,'(a2,1x,a2,f16.8)') nomemet(1),nomemet(2),ener
!                     do i=1,nat
!                   write(1000,'(a2,3f16.8)')spec(i),xin(i),yin(i),zin(i)
!                     enddo
!                     call flush(1000)
!                     close(1000)
!                    endif
                    goto 55 
                  else
                    deltae=ener-enerold
                    deltae=deltae/(abs(deltamet2))
                    deltae=deltae/parameterMC
                    expdeltae=dexp(-deltae)
                    CALL RANDOM_NUMBER(rndm)
                    if (expdeltae.gt.rndm) then
                      enerold=ener
                      !eneratt=ener
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,5,'Energy is HIGHER, accepted',(string(j),j=1,norb)
                      do i=1,nat
                         xin(i)=xat(i)*arete(1)
                         yin(i)=yat(i)*arete(1)
                         zin(i)=zat(i)*arete(1)
                         c(i,1)=xin(i)
                         c(i,2)=yin(i)
                         c(i,3)=zin(i)
                      enddo
                      goto 55 
                    else
!        write(2,'(1F14.6,I20,3I6,2X,A30,2X,500I1)') ener,nmet1,         &
!     &nmet2,k,5,'Energy is HIGHER,  refused',(string(j),j=1,norb)
                      ener=enerold
                      string(n1)=n1old
                      string(n2)=n2old
                      string(n3)=n3old
                      string(n4)=n4old
                      string(n5)=n5old
                      do i=1,nat
                        spec(i)=specsave(i)
                        c(i,1)=csave(i,1)
                        c(i,2)=csave(i,2)
                        c(i,3)=csave(i,3)
                      enddo
                      goto 55 
                    endif
                  endif
                else
                  string(n1)=n1old
                  string(n2)=n2old
                  string(n3)=n3old
                  string(n4)=n4old
                  string(n5)=n5old
                endif
              endif
              endif
              endif
              endif
              endif
              endif
              endif
              endif
              endif
              endif
            enddo ! Until it finds a proper combination
 55         continue

          endif
 
        enddo 
 
 113    continue
 
        index_l=index_l+1
 
        if (index_l.eq.8) goto 111
 
        goto 110
 
 111    continue

        do i=1,nat
           xin(i)=c(i,1)
           yin(i)=c(i,2)
           zin(i)=c(i,3)
           xat(i)=xin(i)/arete(1)
           yat(i)=yin(i)/arete(1)
           zat(i)=zin(i)/arete(1)
           if (spec(i).eq.nomemet(1)) then !assigning chemical species
              itype(i)=1
           else
              itype(i)=2
           endif
        enddo

!        close (2,status='keep')
!        close (99,status='keep')

        ind_tot=ind_tot+8*nsteps
                
        end subroutine grouping 
!***************************************
       SUBROUTINE init_random_seed()
         INTEGER :: i, n, clock
         INTEGER, DIMENSION(:), ALLOCATABLE :: seed
       
         CALL RANDOM_SEED(size = n)
         ALLOCATE(seed(n))
       
         CALL SYSTEM_CLOCK(COUNT=clock)
       
         seed = clock + 37 * (/ (i - 1, i = 1, n) /)
         CALL RANDOM_SEED(PUT = seed)
       
         DEALLOCATE(seed)
       END SUBROUTINE
!***************************************
