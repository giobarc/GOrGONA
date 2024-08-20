        program orbits 
        implicit none 

        integer::i
        character(len=10),dimension(1000)::a1,a2,a3,a4,a5 

        open (1,file='res-pures.dat',status='old')
        do i=1,2
          read(1,'(A4,A1,A2,A1,A8)') a1(i),a2(i),a3(i),a4(i),a5(i)
        enddo

        close(1)
         
        open (1,file='023_copy.sh',status='replace')
        write(1,'(A11)') '#!/bin/bash'
        do i=1,2
          write(1,'(A8,A4,A1,A2,A1,A8,A2)') 'cp -r ./',a1(i),'/',a3(i),'&
     &/',a5(i),' .'
        enddo
        close(1)
        
        call system ('chmod a+x 023_copy.sh')
        end program orbits
