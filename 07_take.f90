        program orbits 
        implicit none 

        integer::ns,i
        integer,dimension(1000)::n1,n2,n3,n4
        real(kind=8),dimension(1000)::f1,f2 
        character(len=10),dimension(1000)::a1 

        ns=0
        open (1,file='result.out',status='old')
        do i=1,1000
          read(1,10,end=2) n1(i),f1(i),f2(i),n2(i),n3(i),n4(i),a1(i) 
          ns=ns+1
        enddo

  10    format(I6,2F10.4,3I6,4X,A10)
 
   2    continue
        close(1)
         
        call system ('rm -r DBFN')
        call system ('mkdir DBFN')

        open (1,file='08_extract.sh',status='replace')
        write(1,'(A11)') '#!/bin/bash'
        do i=1,ns
          if (n2(i).ne.0) then
            if (n1(i).eq.0) then
              write(1,'(A6,I3,A2,I1,A8,I1,A1,A10,A14,I1,A4)') 'cp ./0',n&
     &2(i),'/0',n3(i),'/run_000',n1(i),'/',a1(i),' ./DBFN/MIN000',n1(i),&
     &'.xyz'
            elseif ((n1(i).gt.0).and.(n1(i).lt.100)) then
              write(1,'(A6,I3,A2,I1,A7,I2,A1,A10,A13,I2,A4)') 'cp ./0',n&
     &2(i),'/0',n3(i),'/run_00',n1(i),'/',a1(i),' ./DBFN/MIN00',n1(i),'.&
     &xyz'
            elseif ((n1(i).ge.100).and.(n1(i).lt.1000)) then
              write(1,'(A6,I3,A2,I1,A6,I3,A1,A10,A12,I3,A4)') 'cp ./0',n&
     &2(i),'/0',n3(i),'/run_0',n1(i),'/',a1(i),' ./DBFN/MIN0',n1(i),'.xy&
     &z'
            elseif ((n1(i).ge.1000).and.(n1(i).lt.10000)) then
              write(1,'(A6,I3,A2,I1,A5,I4,A1,A10,A11,I4,A4)') 'cp ./0',n&
     &2(i),'/0',n3(i),'/run_',n1(i),'/',a1(i),' ./DBFN/MIN',n1(i),'.xyz'
            endif
          else
            if (n4(i).lt.10) then
              if (n1(i).eq.0) then
                write(1,'(A15,I1,A8,I1,A1,A10,A14,I1,A4)') 'cp ./motifs/&
     &000',n4(i),'/run_000',n1(i),'/',a1(i),' ./DBFN/MIN000',n1(i),'.xyz&
     &'
              elseif ((n1(i).gt.0).and.(n1(i).lt.100)) then
                write(1,'(A15,I1,A7,I2,A1,A10,A13,I2,A4)') 'cp ./motifs/&
     &000',n4(i),'/run_00',n1(i),'/',a1(i),' ./DBFN/MIN00',n1(i),'.xyz'
              elseif ((n1(i).ge.100).and.(n1(i).lt.1000)) then
                write(1,'(A11,I1,A6,I3,A1,A10,A12,I3,A4)') 'cp ./motifs/&
     &000',n4(i),'/run_0',n1(i),'/',a1(i),' ./DBFN/MIN0',n1(i),'.xyz'
              elseif ((n1(i).ge.1000).and.(n1(i).lt.10000)) then
                write(1,'(A15,I1,A5,I4,A1,A10,A11,I4,A4)') 'cp ./motifs/&
     &000',n4(i),'/run_',n1(i),'/',a1(i),' ./DBFN/MIN',n1(i),'.xyz'
              endif
            elseif ((n4(i).ge.10).and.(n4(i).lt.100)) then
              if (n1(i).eq.0) then
                write(1,'(A14,I2,A8,I1,A1,A10,A14,I1,A4)') 'cp ./motifs/&
     &00',n4(i),'/run_000',n1(i),'/',a1(i),' ./DBFN/MIN000',n1(i),'.xyz'
              elseif ((n1(i).gt.0).and.(n1(i).lt.100)) then
                write(1,'(A14,I2,A7,I2,A1,A10,A13,I2,A4)') 'cp ./motifs/&
     &00',n4(i),'/run_00',n1(i),'/',a1(i),' ./DBFN/MIN00',n1(i),'.xyz'
              elseif ((n1(i).ge.100).and.(n1(i).lt.1000)) then
                write(1,'(A14,I2,A6,I3,A1,A10,A12,I3,A4)') 'cp ./motifs/&
     &00',n4(i),'/run_0',n1(i),'/',a1(i),' ./DBFN/MIN0',n1(i),'.xyz'
              elseif ((n1(i).ge.1000).and.(n1(i).lt.10000)) then
                write(1,'(A14,I2,A5,I4,A1,A10,A11,I4,A4)') 'cp ./motifs/&
     &00',n4(i),'/run_',n1(i),'/',a1(i),' ./DBFN/MIN',n1(i),'.xyz'
              endif
            elseif ((n4(i).ge.100).and.(n4(i).lt.1000)) then
              if (n1(i).eq.0) then
                write(1,'(A13,I3,A8,I1,A1,A10,A14,I1,A4)') 'cp ./motifs/&
     &0',n4(i),'/run_000',n1(i),'/',a1(i),' ./DBFN/MIN000',n1(i),'.xyz'
              elseif ((n1(i).gt.0).and.(n1(i).lt.100)) then
                write(1,'(A13,I3,A7,I2,A1,A10,A13,I2,A4)') 'cp ./motifs/&
     &0',n4(i),'/run_00',n1(i),'/',a1(i),' ./DBFN/MIN00',n1(i),'.xyz'
              elseif ((n1(i).ge.100).and.(n1(i).lt.1000)) then
                write(1,'(A13,I3,A6,I3,A1,A10,A12,I3,A4)') 'cp ./motifs/&
     &0',n4(i),'/run_0',n1(i),'/',a1(i),' ./DBFN/MIN0',n1(i),'.xyz'
              elseif ((n1(i).ge.1000).and.(n1(i).lt.10000)) then
                write(1,'(A13,I3,A5,I4,A1,A10,A11,I4,A4)') 'cp ./motifs/&
     &0',n4(i),'/run_',n1(i),'/',a1(i),' ./DBFN/MIN',n1(i),'.xyz'
              endif
            endif
          endif
        enddo
        close(1)
        
        call system ('chmod a+x 08_extract.sh')
        end program orbits
