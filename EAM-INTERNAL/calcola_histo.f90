SUBROUTINE CALCOLA_HISTO
  
  !***********************
  !This subroutine builds the histogram of the order parameter
  !
  !***********************
  
  USE HISTO
  
  Implicit none
  
  Integer :: i,j
  
  !every renew_histo MC steps, histogram is renewed:
  if (mod(nisto,renew_histo).eq.0) then
     do i=0,nbars(1)
     do j=0,nbars(2)
        histogram(i+1,j+1)=0
     enddo
     enddo
     write (*,*) 'Histogram has been renewed at step',nisto
     nisto=0.
  endif
  
  !Histogram is filled with data coming from the current cluster:
  !histogram is NOT a simmetric matrix
  do i=0,nbars(1)-1
  do j=0,nbars(2)-1
     !first order parameter range:
     hi_l(1)=hist_step(1)*i+start_histo(1)
     hi_s(1)=hi_l(1)+hist_step(1)
     !second order parameter range:
     hi_l(2)=hist_step(2)*j+start_histo(2)
     hi_s(2)=hi_l(2)+hist_step(2)
     if (((parameter_value(1).ge.hi_l(1)).and.(parameter_value(1).lt.hi_s(1))) .and.  &
         &(parameter_value(2).ge.hi_l(2)).and.(parameter_value(2).lt.hi_s(2))) then
        histogram(i+1,j+1)=histogram(i+1,j+1)+1
     endif
     histogram(j+1,i+1)=histogram(i+1,j+1)
  enddo
  enddo
  
  nisto=nisto+1
  
  do i=1,nbars(1)+1
  do j=1,nbars(2)+1
     perc_histo(i,j)=histogram(i,j)/float(nisto)
  enddo
  enddo
  
  
  !if (mod(nisto,10) .eq. 0)then
  !open(100,file='isto.out',status='unknown',position='append')
  !do i=1,nbars(1)+1
  !write(100,*)i,perc_histo(i,1)
  !enddo
  !call flush(100)
  !close(100)
  !endif
  
  
  
END SUBROUTINE CALCOLA_HISTO
