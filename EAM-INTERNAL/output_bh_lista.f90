SUBROUTINE OUTPUT_LISTA(out_call,nlista,ener_or_para,para)

  USE BH_MOD, only : ind_tot,nat3d,nomemet,spec,xin,yin,zin
  USE FUNCTION_PAR, only: ener, e_met_mgo
  
  implicit none
  !local variables
  Integer, Intent(IN) :: nlista
  Real(8), Intent(IN) :: ener_or_para,para
  Integer :: i
  Character*24 :: l_file_name
  Logical :: out_call

    !if this is the first call to output_lista, the file where the minima files are listed has to be initialized:
    if(out_call)then
      open(600,file='lista_minimi.out',status='unknown')
      write(600,*)'# MC step,   nlista, Etot, E met-ox, E met-met, order par(1)'
      call flush(600)
      close(600)
    endif
    open(600,file='lista_minimi.out',status='unknown',position='append')
    write(600,'(2i6,3f16.3,f10.5)')ind_tot,nlista,ener,e_met_mgo, ener-e_met_mgo,para
    call flush(600)
    close(600)

    !the coordinates of the new minimum are going to be stored in min###.xyz:
    write(l_file_name,'(A3,I3.3,A4)')'min',nlista,'.xyz'
    open(1000,file=l_file_name,status='unknown')
    write(1000,*)nat3d
    write(1000,'(a2,1x,a2,3f16.6,f10.5)')nomemet(1),nomemet(2),ener,e_met_mgo,ener-e_met_mgo,para
    do i=1,nat3d
       write(1000,'(a2,3f16.8)')spec(i),xin(i),yin(i),zin(i)
    enddo
    call flush(1000)
    close(1000)


END SUBROUTINE OUTPUT_LISTA
