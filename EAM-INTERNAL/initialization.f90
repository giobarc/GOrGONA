SUBROUTINE INITIALIZATION
  
  USE PARAMETERS
  USE BH_MOD
  USE HISTO
  USE FUNCTION_PAR
  
  implicit none
  
  integer :: iwk,ind
  integer :: lt1neigh,lt3neigh
  integer :: number_of_parameters
!debug
  integer :: n_pass,n_ini=0
  number_of_parameters=3*nat3d

  !if there is a seed:
  IF (choice_seeding .eq. 'seeded') THEN
     !assigning the coordinates of the seed to the first walker:
     do ind=1,nat3d
        xat(ind)=xin(ind)/arete(1)
        yat(ind)=yin(ind)/arete(1)
        zat(ind)=zin(ind)/arete(1)
        if (spec(ind).eq.nomemet(1)) then !assigning chemical species
           itype(ind)=1
        else
           itype(ind)=2
        endif
     enddo

     n_pass=0
     is_nan=.true.
     do while (is_nan .and. n_pass<2000)
	n_pass=n_pass+1
        call minimization_rgl(number_of_parameters,n_ini)
     enddo
     if( n_pass==2000 )then
        write(*,*) 'subroutine initialization: I can not find minima of seed.in'
        ener=0
     endif
     
     do ind=1,nat3d
        xin(ind)=xat(ind)*arete(1)
        yin(ind)=yat(ind)*arete(1)
        zin(ind)=zat(ind)*arete(1)
     enddo


     call force_rgl
     call comparison(lt1neigh,lt3neigh)
     enerwalker(1)=ener
     esubwalker(1)=e_met_mgo
     call struttura(xin,yin,zin) !structural analysis of the initial configuration (order parameter is evaluated)
     call calcola_histo !initial configuration's order parameter is the first value to fill histogram
     walker_parameter(1,:)=parameter_value(:)
     !saving walkers coordinates in *pos matrix
     do ind=1,nat3d
        xpos(ind,1)=xin(ind)
        ypos(ind,1)=yin(ind)
        zpos(ind,1)=zin(ind)
        specpos(ind,1)=spec(ind)
     enddo

     !saving coordinates, energy and order parameter of the minimized seed in minimized_seed.out
     open(13,file='minimized_seed.out',status='unknown')
     write(13,*)nat3d
     write(13,'(a2,1x,a2,3f16.8)')nomemet(1),nomemet(2),ener, e_met_mgo !,parameter_value(1),parameter_value(2)
     do ind=1,nat3d
     write(13,'(a2,3f16.8)')spec(ind),xpos(ind,1),ypos(ind,1),zpos(ind,1)
     enddo
     call flush(13)
     close(13)

     !the other walkers are initialized randomly:
     do iwk=2,nwalkers
        lt3neigh=1
        do while (lt3neigh.gt.0)
           call initial_configuration !atoms randomly displaced
           do ind=1,nat3d
              xat(ind)=xin(ind)/arete(1)
              yat(ind)=yin(ind)/arete(1)
              zat(ind)=zin(ind)/arete(1)
              if (spec(ind).eq.nomemet(1)) then !assigning chemical species
                 itype(ind)=1
              else
                 itype(ind)=2
              endif
           enddo
 	   call minimization_rgl(number_of_parameters,n_ini)
           do ind=1,nat3d
              xin(ind)=xat(ind)*arete(1)
              yin(ind)=yat(ind)*arete(1)
              zin(ind)=zat(ind)*arete(1)
           enddo
           call comparison(lt1neigh,lt3neigh)!cicle again if 
           !there are atoms with less than 3 nearest neighbors 
        enddo
        
        !when a random and compact cluster has been built,
        !initial coordinates are called xin,yin,zin
        enerwalker(iwk)=ener
        esubwalker(iwk)=e_met_mgo
        call struttura(xin,yin,zin) !structural analysis of the initial configuration (order parameter is evaluated)
        call calcola_histo !initial configuration's order parameter is the first value to fill histogram
        walker_parameter(iwk,:)=parameter_value(:)
        
        !saving walkers coordinates in *pos matrix
        do ind=1,nat3d
           xpos(ind,iwk)=xin(ind)
           ypos(ind,iwk)=yin(ind)
           zpos(ind,iwk)=zin(ind)
           specpos(ind,iwk)=spec(ind)
        enddo
        
     enddo
     
     !if there is no seed:
  ELSEIF (choice_seeding .eq. 'unseeded') THEN
     !nwalkers clusters are initialized:
     do iwk=1,nwalkers
        lt3neigh=1
        
        do while (lt3neigh.gt.0)

           call initial_configuration !atoms randomly displaced
           do ind=1,nat3d
              xat(ind)=xin(ind)/arete(1)
              yat(ind)=yin(ind)/arete(1)
              zat(ind)=zin(ind)/arete(1)
              if (spec(ind).eq.nomemet(1)) then !assigning chemical species
                 itype(ind)=1
              else
                 itype(ind)=2
              endif
           enddo
  	   call minimization_rgl(number_of_parameters,n_ini)

           do ind=1,nat3d
              xin(ind)=xat(ind)*arete(1)
              yin(ind)=yat(ind)*arete(1)
              zin(ind)=zat(ind)*arete(1)
           enddo
	   call force_rgl
           call comparison(lt1neigh,lt3neigh)!cicle again if 

           !there are atoms with less than 3 nearest neighbors 
        enddo
        !when a random and compact cluster has been built,
        !initial coordinates are called xin,yin,zin
        enerwalker(iwk)=ener
        esubwalker(iwk)=e_met_mgo
        call struttura(xin,yin,zin) !structural analysis of the initial configuration (order parameter is evaluated)
        call calcola_histo !initial configuration's order parameter is the first value to fill histogram
        walker_parameter(iwk,:)=parameter_value(:)
        !saving walkers coordinates in *pos matrix
        do ind=1,nat3d
           xpos(ind,iwk)=xin(ind)
           ypos(ind,iwk)=yin(ind)
           zpos(ind,iwk)=zin(ind)
           specpos(ind,iwk)=spec(ind)
        enddo
        
     enddo

  ELSE
	write(*,*)'I accept either seeded or unseeded. You wrote: ',choice_seeding
	stop
  ENDIF

END SUBROUTINE INITIALIZATION
