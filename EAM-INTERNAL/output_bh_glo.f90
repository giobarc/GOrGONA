SUBROUTINE OUTPUT_GLO(iwalk,output_instruction)
  
  USE PARAMETERS
  USE BH_MOD
  USE FUNCTION_PAR
  
  implicit none
  !local variables
  Integer :: ind,ncall_output,iwalk
  Integer :: output_instruction
             ! = 1 first call to the output subroutine, before MC loop starts - file headings
             ! = 2 current lowest energy minimum has been found - its coordinates are saved
             ! = 3 energy and order parameter at the current MC step are recorded
             ! = 4 every 1000 steps the percentage of accepted/rejected moves are recorded
  Character(LEN=100) :: file_name,glo_file_name
  Logical :: label_output_coordinates  

  IF(output_instruction .eq. 1)THEN !program writes file headers for ener_oparameter.out, moves.out, ener_best.out
     open(100,file='ener_oparameter.out',status='unknown')
     write(100,*)'#MC step, ener_att-ener_old, ener,order par(1), order par(2),accflag'
     if(output_walkers .eq. 'si')then
 	if(nwalkers .gt. 9)then
            write(*,*)'####################################################'
            write(*,*)'ERROR'
            write(*,*)'Too many walkers to produce separate output files'
            write(*,*)'####################################################'
            open(22,file='errors_bh.out',status='unknown',position='append')
            write(22,*)'####################################################'
            write(22,*)'ERROR'
            write(22,*)'Too many walkers to produce separate output files'
            write(22,*)'####################################################'
            close(22)
            stop
 	endif
	do ind=1,nwalkers
           write(file_name,'(A15,I2.2,A4)')'ener_oparameter',ind,'.out'
           open(100+ind,file=file_name,status='unknown')
           write(100+ind,*)'# walker n.',ind 
           write(100+ind,*)'# MC step, ener_att-ener_old,ener,order par(1), order par(2), accflag, label'
	enddo
     endif
     open(200,file='moves.out',status='unknown')
     write(200,*)'# MC step, % accepted moves, % true moves, % rejected moves'
     close(200)
     open(400,file='ener_best.out',status='unknown')
     write(400,*)'#MC step, glo#.xyz, Etot, order par(1), order par(2)'
     close(400)
     return
  ELSEIF(output_instruction .eq. 2)THEN
     !updating file ener_best.out
     open(400,file='ener_best.out',status='unknown',position='append')
     write(400,'(2i10,3f16.6)')ind_tot,nglo-1,ener,&
&walker_parameter(iwalk,1),walker_parameter(iwalk,2)
     call flush(400)
     close(400)
     !glo###.xyz file name:
     write(glo_file_name,'(A3,I4.4,A4)')'glo',nglo-1,'.xyz'
     open(1000,file=glo_file_name,status='unknown')
     write(1000,*)nat3d
     write(1000,'(a2,1x,a2,1f16.6,2f10.5)')nomemet(1),nomemet(2),ener,&
&walker_parameter(iwalk,1),walker_parameter(iwalk,2)
     do ind=1,nat3d
     write(1000,'(a2,3f16.8)')spec(ind),xin(ind),yin(ind),zin(ind)
     enddo
     call flush(1000)
     close(1000)
  ELSEIF(output_instruction .eq. 3)THEN
     !write(100,'(i10,2f16.6)')ind_tot,ener_att,ener_old
     write(100,'(i10,2f16.6,2f12.5,i6,i4,2x,a6)')ind_tot,ener_att-ener_old,ener,&
&walker_parameter(iwalk,1),walker_parameter(iwalk,2),accflag,isucc,string_move
     call flush(100)
     if(output_walkers .eq. 'si')then !different files for different 
        !walkers have been requested in the input file
        write(100+iwalk,'(i10,2f16.6,2f12.5,i6,i4,2x,a6)')ind_tot,ener_att-ener_old,ener,&
&walker_parameter(iwalk,1),walker_parameter(iwalk,2),accflag,isucc,string_move
        call flush(100+iwalk)
     endif
  ELSEIF(output_instruction .eq. 4)THEN
     open(200,file='moves.out',status='unknown',position='append')
     write(200,'(i8,3f8.2)')ind_tot,100*n_accepted_moves/dfloat(ind_tot), &
                            100*n_true_moves/dfloat(ind_tot),100*n_rejected_moves/dfloat(ind_tot)
     call flush(200)
     close(200)
  ENDIF

END SUBROUTINE OUTPUT_GLO
