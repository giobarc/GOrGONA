SUBROUTINE COMPARISON(lt1neigh,lt3neigh)

  !********************
  !Last change: June 6th, 2013 by Giulia Rossi
  !This routine checks how many atoms in the cluster
  !have less than 1 nearest neighbors (lt1neigh) or 3 neighbors (lt3neigh)
  !Moreover, for each atom it calculates and stores in global variables the total number of 
  !neighbors (quantivicini(1:nat3d)), the number of neighbors of the same species (quantivi_eq(1:nat3d))
  !and the number of neighbors of the other species (quantivi_neq(1:nat3d))
  !********************

  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR

  implicit none

  !local variables
  integer :: ind,ind1,ind2
  real(8) :: ccc,nnmax2
  integer :: lt1neigh,lt3neigh

  quantivicini(1:nat3d)=0
  quantivi_eq(1:nat3d)=0
  quantivi_neq(1:nat3d)=0

  nnmax2=dmax1(nn(1),nn(2))*dmax1(nn(1),nn(2))!square of the biggest nearest neighbor distance
  nnmax2=nnmax2+0.3d0*nnmax2 !14% tolerance on nnmax2, it's enough for systems with a big size mismatch
  !write(*,*)'nnmax2=',nnmax2

  !If the cluster is homogeneous, the routine assigns to each atom a number of neighbors, stored in "quantivicini"
  IF((nspec1 .eq. 0) .or. (nspec2 .eq. 0))THEN
  !loop on all atoms
  do ind1=1,nat3d
     do ind2=ind1+1,nat3d
        ccc= (xin(ind1)-xin(ind2))**2+(yin(ind1)-yin(ind2))**2+(zin(ind1)-zin(ind2))**2
        if(ccc.le.nnmax2) then !the number of neighbors of both atoms needs to be updated:
           quantivicini(ind1)=quantivicini(ind1)+1
           quantivicini(ind2)=quantivicini(ind2)+1
        endif
     enddo
  enddo

  !If the cluster is NOT homogeneous, the routine calculates the total and partial (same species, different species) 
  !number of neighbors of each atom
  ELSE

  !loop on the small atoms
  do ind1=1,nsmall

     do ind2=1,nsmall !looking for neighbors of the same species
        ccc= (xin(list_small(ind1))-xin(list_small(ind2)))**2+(yin(list_small(ind1))-yin(list_small(ind2)))**2+&
        & (zin(list_small(ind1))-zin(list_small(ind2)))**2
        if(ccc.le.nnmax2) then !the number of neighbors of both atoms needs to be updated:
           quantivicini(list_small(ind1))=quantivicini(list_small(ind1))+1
           quantivi_eq(list_small(ind1))=quantivi_eq(list_small(ind1))+1
        endif
     enddo

     do ind2=1,nbig !looking for neighbors of the other species
        ccc= (xin(list_small(ind1))-xin(list_big(ind2)))**2+(yin(list_small(ind1))-yin(list_big(ind2)))**2+&
        & (zin(list_small(ind1))-zin(list_big(ind2)))**2
        if(ccc.le.nnmax2) then !the number of neighbors of both atoms needs to be updated:
           quantivicini(list_small(ind1))=quantivicini(list_small(ind1))+1
           quantivicini(list_big(ind2))=quantivicini(list_big(ind2))+1
           quantivi_neq(list_small(ind1))=quantivi_neq(list_small(ind1))+1
           quantivi_neq(list_big(ind2))=quantivi_neq(list_big(ind2))+1
        endif
     enddo
  enddo

  !loop on the big atoms
  do ind1=1,nbig

     do ind2=1,nbig !looking for neighbors of the same species
        ccc= (xin(list_big(ind1))-xin(list_big(ind2)))**2+(yin(list_big(ind1))-yin(list_big(ind2)))**2+&
        & (zin(list_big(ind1))-zin(list_big(ind2)))**2
        if(ccc.le.nnmax2) then !the number of neighbors of both atoms needs to be updated:
           !write(*,*)list_big(ind1),'and',list_big(ind2),'are neighbors'
           quantivicini(list_big(ind1))=quantivicini(list_big(ind1))+1
           quantivi_eq(list_big(ind1))=quantivi_eq(list_big(ind1))+1
        endif
     enddo

  enddo

  !Taking into account that no atom is neighbor to itself:
  quantivicini(1:nat3d)=quantivicini(1:nat3d)-1
  quantivi_eq(1:nat3d)=quantivi_eq(1:nat3d)-1

  ENDIF

!To check neighbors' count
!  do ind1=1,nat3d
!     write(*,*)'atomo ',ind1,'ha',quantivicini(ind1),'vicini'
!     write(*,*)'dei quali',quantivi_eq(ind1),'della sua stessa specie, e'
!     write(*,*)quantivi_neq(ind1),'di specie diversa.'
!  enddo
!stop

  lt1neigh=0
  lt3neigh=0
  do ind1=1,nat3d
     if (quantivicini(ind1).le.1) then
        lt1neigh=lt1neigh+1
        !write(*,*)'atomo',ind1,'ha le 1 neigh',quantivicini(ind1)
     elseif(quantivicini(ind1).le.3) then
        lt3neigh=lt3neigh+1
     endif
  enddo

END SUBROUTINE COMPARISON
