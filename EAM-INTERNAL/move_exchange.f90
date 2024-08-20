SUBROUTINE MOVE_EXCHANGE(temperature_exchange)
  
  !****************************************************************************
  !By "exchange" we mean that the positions of two atoms of different species are swapped
  !Several types of exchange move are possibile:
  !
  !1. Exchange move ALL
  !Two atoms of different species are randomly picked and their positions are swapped
  !
  !2. Exchange move SUR
  !Two surface atoms are of different species are randomly picked and their positions swapped
  !
  !3. Exchange move CSH
  !This move favors core-shell ordering. It depends on two parameters: bignn_csh and smallnn_csh. 
  !A small surface atom, that have at least bignn_csh neighbors, is swapped with a big volume atom that 
  !has at least smallnn_csh small neighbors.
  !
  !4. Exchange move INT
  !This move swaps atoms that have more than bignn (for small atoms) or more than smallnn (for big atoms) 
  !neighbors. It thus acts on atoms that are isolated within or at the interface with the other species
  !
  !5. Exchange move SEP
  !This move works in two opposite ways:
  !small_inside=T: the move swaps a small atom that is in the volume and has many big neighbors, 
  !or is in the surface and has several big neighbors, with a big atom that is in the volume and has many small neighbors
  !small_inside=F: as before, but the roles of the small and big atoms are swapped. This is useful 
  !for binary systems where the small atoms are not the more cohesive, such as Cu-Rh

  !6. Exchange move MIX
  !This move swaps atoms of different species that have many atom of their same species
  !7.  Exchange move MUL
  !****************************************************************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  
  Real(8), intent(OUT) :: temperature_exchange
  !local variables
  Integer, Parameter :: npvmax=15
  Integer :: iatom1,iatom2,aa,cont,ind,indd,npv,nclass1,nclass2,aa1,aa2,ip1,ip2
  Integer :: pochivicini,pochivicini_counter,tantiAg_counter,pochiCu_counter
  Integer :: pochivicini_array(nat3d),tantivicini_Ag(nat3d),pochivicini_Cu(nat3d)
  Integer :: who_small(0:npvmax,nat3d),who_big(0:npvmax,nat3d)
  Integer :: cpv_small(0:npvmax),cpv_big(0:npvmax),ipos1(nat3d),ipos2(nat3d)
  Integer :: imult1(nat3d),imult2(nat3d)
  Integer :: iequal1,iequal2,nspecmin,nchosen,ntotchosen
  Real(8), Parameter :: epsicl=1.d-2
  Real(8) :: weight_small(-1:npvmax),weight_big(-1:npvmax)
  Real(8) :: wbkl_small(0:nat3d),wbkl_small_tot(0:nat3d),wbkl_big(0:nat3d),wbkl_big_tot(0:nat3d)
  Real(8) :: proba(0:nat3d)
  Real(8) :: pr1, pr2, pr3, pr4, pr5, pr6, pr7, pr8!probabilities of the different kinds of exchange moves
  Real(8) :: xtemp,ytemp,ztemp,bb1,bb2,normprob,rexch
  Logical :: couple

  !Exchange moves all, sur, csh, int will create 2 lists, one containing small atoms and one containing big atoms with the desired 
  !characteristics (ex: small atoms being on the surface, big atoms being inside the volume). Then, atoms to be swapped 
  !will be randomly extracted from these lists.
  Integer  :: exlist_small(nsmall), exlist_big(nbig) !these are maximum dimensions for the lists, that might be shorter
  Integer  :: exlist_small_dim, exlist_big_dim 
  Integer  :: lt1neigh,lt3neigh !needed by the call to comparison.f90
  Character(3) :: exchange_type

 
  exlist_small(:)=0
  exlist_big(:)=0

  !Choice of the type of exchange move to be executed
  pr1=prob_exch_all !           zzz < pr1   -->   exchange move ALL
  pr2=pr1+prob_exch_sur  ! pr1 < zzz < pr2   -->   exchange move SUR
  pr3=pr2+prob_exch_csh  ! pr2 < zzz < pr3   -->   exchange move CSH
  pr4=pr3+prob_exch_int  ! pr3 < zzz < pr4   -->   exchange move INT
  pr5=pr4+prob_exch_sep  ! pr4 < zzz < pr5   -->   exchange move SEP
  pr6=pr5+prob_exch_mix  ! pr5 < zzz < pr6   -->   exchange move MIX
  pr7=pr6+prob_exch_mult ! pr6 < zzz < pr7   -->   exchange move MULT
  pr8=pr7+prob_exch_gru  ! pr7 < zzz < pr8   -->   exchange move GRU 

  call ranmar(zzz)
        	if(zzz.le.pr1) then
        	   exchange_type='all'
        	elseif((zzz.gt.pr1).and.(zzz.le.pr2)) then
        	   exchange_type='sur'
		elseif((zzz.gt.pr2).and.(zzz.le.pr3)) then
        	   exchange_type='csh'
        	elseif((zzz.gt.pr3).and.(zzz.le.pr4)) then
        	   exchange_type='int'
		elseif((zzz.gt.pr4).and.(zzz.le.pr5)) then
        	   exchange_type='sep'
        	elseif((zzz.gt.pr5).and.(zzz.le.pr6)) then
        	   exchange_type='mix'
        	elseif((zzz.gt.pr6).and.(zzz.le.pr7)) then
        	   exchange_type='mul'
        	elseif((zzz.gt.pr7).and.(zzz.le.pr8)) then
        	   exchange_type='gru'
        	endif

  !Informations on which are the small and the big atoms are stored in the global variable lists list_big(1:nbig) and list_small(1:nsmall)
  SELECT CASE (exchange_type)

          CASE ('gru')

               call grouping

          CASE ('all')
                
                temperature_exchange=temp_all
                call ranmar(zzz)
                iatom1=list_small(min(int(zzz*nsmall+1),nsmall))
                call ranmar(zzz)
                iatom2=list_big(min(int(zzz*nbig+1),nbig))

          CASE ('sur')
                

                temperature_exchange=temp_sur
                call comparison(lt1neigh,lt3neigh)

                !create small atom list
                exlist_small_dim=0
                do ind=1,nsmall
                   if(quantivicini(list_small(ind)) .le. 9)then 
                     exlist_small_dim=exlist_small_dim+1
                     exlist_small(exlist_small_dim)=list_small(ind)
                   endif
                enddo
                if(exlist_small_dim .eq. 0)then
                  write(*,*)'WARNING: at step ',ind_tot,' it was impossible to perform the exchange move, type sur'
                  return
                endif

                !create big atom list
                exlist_big_dim=0
                do ind=1,nbig
                   if(quantivicini(list_big(ind)) .le. 9)then 
                     exlist_big_dim=exlist_big_dim+1
                     exlist_big(exlist_big_dim)=list_big(ind)
                   endif
                enddo
                if(exlist_big_dim .eq. 0)then
                  write(*,*)'WARNING: at step ',ind_tot,' it was impossible to perform the exchange move, type sur'
                  return
                endif

                !random extraction from lists
                call ranmar(zzz)
                iatom1=exlist_small(min(int(zzz*exlist_small_dim+1),exlist_small_dim))
                call ranmar(zzz)
                iatom2=exlist_big(min(int(zzz*exlist_big_dim+1),exlist_big_dim))


          CASE ('csh')

                temperature_exchange=temp_csh

                call comparison(lt1neigh,lt3neigh)
                !create small atom list
                exlist_small_dim=0
                do ind=1,nsmall
                   if((quantivicini(list_small(ind)) .le. 9) .and. (quantivi_neq(list_small(ind)) .gt. bignn_csh))then 
                     exlist_small_dim=exlist_small_dim+1
                     exlist_small(exlist_small_dim)=list_small(ind)
                   endif
                enddo

                if(exlist_small_dim .eq. 0)then
                  write(*,*)'WARNING: at step ',ind_tot,' it was impossible to perform the exchange move, type csh'
                  return
                endif

                !create big atom list
                exlist_big_dim=0
                do ind=1,nbig
                   if((quantivicini(list_big(ind)) .ge. 12) .and. (quantivi_neq(list_big(ind)) .gt. smallnn_csh))then 
                     exlist_big_dim=exlist_big_dim+1
                     exlist_big(exlist_big_dim)=list_big(ind)
                   endif
                enddo

                if(exlist_big_dim .eq. 0)then
                  write(*,*)'WARNING: at step ',ind_tot,' it was impossible to perform the exchange move, type csh'
                  return
                endif

                !random extraction from lists
                call ranmar(zzz)
                iatom1=exlist_small(min(int(zzz*exlist_small_dim+1),exlist_small_dim))
                call ranmar(zzz)
                iatom2=exlist_big(min(int(zzz*exlist_big_dim+1),exlist_big_dim))

          CASE ('int')

                temperature_exchange=temp_int

                call comparison(lt1neigh,lt3neigh)

                !create small atom list
                exlist_small_dim=0
                do ind=1,nsmall
                   if(quantivi_neq(list_small(ind)) .gt. bignn_int)then 
                     exlist_small_dim=exlist_small_dim+1
                     exlist_small(exlist_small_dim)=list_small(ind)
                   endif
                enddo

                if(exlist_small_dim .eq. 0)then
                  write(*,*)'WARNING: at step ',ind_tot,' it was impossible to perform the exchange move, type int'
                  return
                endif

                !create big atom list
                exlist_big_dim=0
                do ind=1,nbig
                   if(quantivi_neq(list_big(ind)) .gt. smallnn_int)then 
                     exlist_big_dim=exlist_big_dim+1
                     exlist_big(exlist_big_dim)=list_big(ind)
                   endif
                enddo

                if(exlist_big_dim .eq. 0)then
                  write(*,*)'WARNING: at step ',ind_tot,' it was impossible to perform the exchange move, type int'
                  return
                endif

                !random extraction from lists
                call ranmar(zzz)
                iatom1=exlist_small(min(int(zzz*exlist_small_dim+1),exlist_small_dim))
                call ranmar(zzz)
                iatom2=exlist_big(min(int(zzz*exlist_big_dim+1),exlist_big_dim))

          CASE ('sep')

                temperature_exchange=temp_sep

                call comparison(lt1neigh,lt3neigh)

                wbkl_small(:)=0.d0
                wbkl_small_tot(:)=0.d0
                wbkl_big(:)=0.d0
                wbkl_big_tot(:)=0.d0

                IF(small_inside)THEN

                  !Each atom is assigned a weight.
                  !small atoms will weight more the more they are isolated or on the surface:
                  do ind=1,nsmall
                     if (quantivicini(list_small(ind)).ge.12) then
                        wbkl_small(ind)=1.d0*(quantivi_neq(list_small(ind))+1)**exp_sep
                     else
                        wbkl_small(ind)=1.d0*(quantivi_neq(list_small(ind))+1+12-quantivicini(list_small(ind)))**exp_sep
                     endif
                     wbkl_small_tot(ind)= wbkl_small_tot(ind-1)+wbkl_small(ind)
                  enddo

                  !big atoms will weight more the more they are in contact with small atoms. 
                  do ind=1,nbig
                     wbkl_big(ind)=1.d0*(quantivi_neq(list_big(ind))+1)**exp_sep
                     wbkl_big_tot(ind)= wbkl_big_tot(ind-1)+wbkl_big(ind)
                  enddo

                ELSEIF(.NOT. small_inside)THEN

                  !Each atom is assigned a weight.
                  !big atoms will weight more the more they are isolated or on the surface:
                  do ind=1,nbig
                     if (quantivicini(list_big(ind)).ge.12) then
                        wbkl_big(ind)=1.d0*(quantivi_neq(list_big(ind))+1)**exp_sep
                     else
                        wbkl_big(ind)=1.d0*(quantivi_neq(list_big(ind))+1+12-quantivicini(list_big(ind)))**exp_sep
                     endif
                     wbkl_big_tot(ind)= wbkl_big_tot(ind-1)+wbkl_big(ind)
                  enddo

                  !small atoms will weight more the more they are in contact with big atoms. 
                  do ind=1,nsmall
                     wbkl_small(ind)=1.d0*(quantivi_neq(list_small(ind))+1)**exp_sep
                     wbkl_small_tot(ind)=wbkl_small_tot(ind-1)+wbkl_small(ind)
                  enddo


                ENDIF


               !Final choice of the atoms to be swapped:
               !The small one:
               call ranmar(zzz)
               bb1=wbkl_small_tot(nsmall)*zzz
               if (bb1.lt.epsicl) bb1=epsicl
               do ind=1,nsmall
                  if ((bb1.gt.wbkl_small_tot(ind-1)).and.(bb1.le.wbkl_small_tot(ind))) iatom1=list_small(ind)
               enddo

               !The big one:
               call ranmar(zzz)
               bb2=wbkl_big_tot(nbig)*zzz
               if (bb2.lt.epsicl) bb2=epsicl
               do ind=1,nbig
                  if ((bb2.gt.wbkl_big_tot(ind-1)).and.(bb2.le.wbkl_big_tot(ind))) iatom2=list_big(ind)
               enddo

          CASE ('mix')

                temperature_exchange=temp_mix

                call comparison(lt1neigh,lt3neigh)

                wbkl_small(:)=0.d0
                wbkl_small_tot(:)=0.d0
                wbkl_big(:)=0.d0
                wbkl_big_tot(:)=0.d0

                ip1=0 
                ip2=0

                !Each atom is assigned a weight.
                !small atoms will weight more the more they are in the volume and/or surrounded by small neighbors
                do ind=1,nsmall
                   if (quantivicini(list_small(ind)).ge.12) then
                      wbkl_small(ind)=1.d0*(quantivi_eq(list_small(ind))+1)**exp_mix
                   else
                      wbkl_small(ind)=1.d0*(quantivi_eq(list_small(ind))+1+12-quantivicini(list_small(ind)))**exp_mix
                   endif
                   wbkl_small_tot(ind)= wbkl_small_tot(ind-1)+wbkl_small(ind)


                enddo

                !big atoms will weight more the more they are in contact with big atoms. 
                do ind=1,nbig
                   wbkl_big(ind)=1.d0*(quantivi_eq(list_big(ind))+1)**exp_mix
                   wbkl_big_tot(ind)= wbkl_big_tot(ind-1)+wbkl_big(ind)
                enddo

               !Final choice of the atoms to be swapped:
               !The small one:
               call ranmar(zzz)
               bb1=wbkl_small_tot(nsmall)*zzz
               if (bb1.lt.epsicl) bb1=epsicl
               do ind=1,nsmall
                  if ((bb1.gt.wbkl_small_tot(ind-1)).and.(bb1.le.wbkl_small_tot(ind))) iatom1=list_small(ind)
               enddo

               !The big one:
               call ranmar(zzz)
               bb2=wbkl_big_tot(nbig)*zzz
               if (bb2.lt.epsicl) bb2=epsicl
               do ind=1,nbig
                  if ((bb2.gt.wbkl_big_tot(ind-1)).and.(bb2.le.wbkl_big_tot(ind))) iatom2=list_big(ind)
               enddo
     CASE('mul')
               temperature_exchange=temp_mult
               ! choose how may atoms are to be exchanged
               nspecmin=min0(nbig,nsmall)
               normprob=1.d0
               proba(0)=0.d0
               proba(1)=1.d0
               do ind=2,nspecmin
                  normprob=normprob+1.d0/dfloat(ind)**exp_mult
                  proba(ind)=proba(ind-1)+1.d0/dfloat(ind)**exp_mult
               enddo
               call ranmar(zzz)
               rexch=normprob*zzz
               do ind=1,nspecmin
                  if ((rexch.gt.proba(ind-1)).and.(rexch.le.proba(ind))) then
                     ntotchosen=ind
                  endif
               enddo
               nchosen=0
               do while (nchosen.lt.ntotchosen)              
                  call ranmar(zzz)
                  iatom1=list_small(min(int(zzz*nsmall+1),nsmall))
                  iequal1=0
                  do ind=1,nchosen
                     if(iatom1.eq.imult1(ind)) iequal1=iequal1+1                     
                  enddo
                  if (iequal1.eq.0) then
                     nchosen=nchosen+1
                     imult1(nchosen)=iatom1
                  endif
               enddo
               nchosen=0
               do while (nchosen.lt.ntotchosen)              
                  call ranmar(zzz)
                  iatom2=list_big(min(int(zzz*nbig+1),nbig))
                  iequal2=0
                  do ind=1,nchosen
                     if(iatom2.eq.imult2(ind)) iequal2=iequal2+1                     
                  enddo
                  if (iequal2.eq.0) then
                     nchosen=nchosen+1
                     imult2(nchosen)=iatom2
                  endif
               enddo
               write (88,*) ind_tot,ntotchosen,nchosen
  END SELECT

  if(exchange_type.eq.'gru') goto 27

  !exchanging coordinates
  if(exchange_type.ne.'mul') then
     xtemp=xin(iatom1)
     ytemp=yin(iatom1)
     ztemp=zin(iatom1)
     xin(iatom1)=xin(iatom2)
     yin(iatom1)=yin(iatom2)
     zin(iatom1)=zin(iatom2)
     xin(iatom2)=xtemp
     yin(iatom2)=ytemp
     zin(iatom2)=ztemp
  else
     do ind=1,ntotchosen
        xtemp=xin(imult1(ind))
        ytemp=yin(imult1(ind))
        ztemp=zin(imult1(ind))
        xin(imult1(ind))=xin(imult2(ind))
        yin(imult1(ind))=yin(imult2(ind))
        zin(imult1(ind))=zin(imult2(ind))
        xin(imult2(ind))=xtemp
        yin(imult2(ind))=ytemp
        zin(imult2(ind))=ztemp
     enddo
  endif

27 continue
  
END SUBROUTINE MOVE_EXCHANGE

 
