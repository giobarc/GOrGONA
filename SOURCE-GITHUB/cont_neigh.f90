SUBROUTINE CONT_NEIGH(iwalk,number_of_neighbours)
  
  !*************************************************
  !This subroutine evaluates the number of neighbors
  !of the iwalk walker, before attempting a move towards
  !another minimum. Cont_neigh acts without any information about 
  !the choice of the order parameter.
  !*************************************************
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  Integer :: iwalk,icc
  Integer :: number_of_neighbours
  Real(8) :: lowlimit(2),uplimit(2)
  
  number_of_neighbours=0
  
  maxenw=1

  do icc=1,nwalkers
     if (icc.ne.iwalk) then
        if (enerwalker(icc).gt.enerwalker(iwalk)) maxenw=0
        uplimit(1)=walker_parameter(icc,1)+widthwalkers(1)
        lowlimit(1)=walker_parameter(icc,1)-widthwalkers(1)
	uplimit(2)=walker_parameter(icc,2)+widthwalkers(2)
        lowlimit(2)=walker_parameter(icc,2)-widthwalkers(2)
        if (((walker_parameter(iwalk,1).ge.lowlimit(1)).and.(walker_parameter(iwalk,1).le.uplimit(1))) &
	& .and. &
	& (walker_parameter(iwalk,2).ge.lowlimit(2)).and.(walker_parameter(iwalk,2).le.uplimit(2))) then
           number_of_neighbours=number_of_neighbours+1
        endif
     endif
  enddo
  
END SUBROUTINE CONT_NEIGH
