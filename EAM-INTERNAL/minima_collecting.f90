Subroutine minima_collecting(l_ener,l_parameter,ind_mc)
  
USE BH_MOD, only: param_min,param_max,parameter_tolerance,parameter_interval,parameter_partitioning

  Implicit None

  Integer, Intent(INOUT) :: ind_mc
  Real(8), Intent(IN)  :: l_ener,l_parameter


  Integer :: nl,nnl
  Real(8) :: nl_max
  Real(8) :: l_best_ener(1000),l_best_param(1000)
  Logical :: label_collecting_minima,out_call

!assegnazione del minimo a una casella dello spazio del parametro d'ordine
!valori del p.o. variabili tra 0 e 1

if(ind_mc .eq. 1)then
  l_best_ener(:)=0.
  l_best_param(:)=1000.
  out_call=.true.
  write(*,*)'MINIMA_COLLECTING',out_call
else
  out_call=.false.
endif

open(22,file='errors_bh.out',status='unknown',position='append')
if(parameter_partitioning .ge. 1000)then
  write(*,*)'WARNING: problema di dimensionamento di l_best_ener in minima_collecting.f90'
  write(22,*)'WARNING: problema di dimensionamento di l_best_ener in minima_collecting.f90'
  close(22)
  stop
endif

parameter_interval=(param_max-param_min)/parameter_partitioning !param_max, param_min and parameter_partitioning are read in input_bh.in

do nl=1,parameter_partitioning

   if ((l_parameter .ge. param_min+(nl-1)*parameter_interval) .and. (l_parameter .lt. param_min+nl*parameter_interval))then

      !controllo se l'energia del nuovo minimo e' inferiore a quello trovato in precedenza in questo intervallo del parametro d'ordine
      !al primo ingresso nell'intervallo "nl" la condizione e' sempre verificata:
        if(l_ener .lt. l_best_ener(nl))then !sto migliorando l'energia di questo param ordine

	  l_best_ener(nl)=l_ener !aggiornamento minima energia per nl
          l_best_param(nl)=l_parameter !aggiornamento p ordine del miglior cluster per nl

	  call output_lista(out_call,nl,l_ener,l_best_param(nl))
	  open(15,file='param_best.out',status='unknown')
	  do nnl=1,parameter_partitioning
	     !scrivo nel file solo i valori del p.o. per i quali ho trovato un cluster:
             if(l_best_ener(nnl) .lt. 0.)write(15,'(i6,3f21.13)')nnl,nnl*parameter_interval,l_best_param(nnl),l_best_ener(nnl)
	  enddo
          call flush(15)
          close(15)

        endif
    endif
enddo

End Subroutine minima_collecting
