         program parallel_orbits 
         implicit none
         integer:: i,ni,nf,ns,j,nn,nmet1,diff,nd
         real(kind=8),dimension(1000):: x,y,z
         character(len=2),dimension(1000):: lab

         open(3,file='range.dat',status='old')
         read(3,*) ni,nf,ns
         do i=ni,nf,ns
           open(1,file='move',status='replace')
           open(4,file='seed.in.temp',status='old')
           read(4,*) nn 
           read(4,*)  
           nmet1=0
           do j=1,nn
             read(4,*) lab(j),x(j),y(j),z(j)
             if (lab(j).eq.'Xx') then
               nmet1=nmet1+1
             endif
           enddo
           diff=i-nmet1
           close(4)
           open(5,file='seed.in',status='replace')
           write(5,*) nn
           write(5,*) 'Xx Yy' 
           if (diff.eq.0) then
             do j=1,nn
               write(5,'(A2,3F10.4)') lab(j),x(j),y(j),z(j)
             enddo
           elseif (diff.gt.0) then
             nd=diff
             do j=1,nn
               if (lab(j).eq.'Xx') then
                 write(5,'(A2,3F10.4)') lab(j),x(j),y(j),z(j)
               else
                 if (nd.ne.0) then
                   write(5,'(A2,3F10.4)') 'Xx',x(j),y(j),z(j)
                   nd=nd-1
                 else
                   write(5,'(A2,3F10.4)') 'Yy',x(j),y(j),z(j)
                 endif
               endif
             enddo
           elseif (diff.lt.0) then
             nd=diff
             do j=1,nn
               if (lab(j).eq.'Yy') then
                 write(5,'(A2,3F10.4)') lab(j),x(j),y(j),z(j)
               else
                 if (nd.ne.0) then
                   write(5,'(A2,3F10.4)') 'Yy',x(j),y(j),z(j)
                   nd=nd+1
                 else
                   write(5,'(A2,3F10.4)') 'Xx',x(j),y(j),z(j)
                 endif
               endif
             enddo
           endif
           close(5)

           if (i.eq.0) then
             write(1,'(A13,I1)') 'mkdir run_000',i
             write(1,'(A24,I1)') 'cp XxYy.in       run_000',i
             write(1,'(A24,I1)') 'cp bh_v07        run_000',i
             write(1,'(A24,I1)') 'cp input_bh.in   run_000',i
             write(1,'(A24,I1)') 'cp input_gru.dat run_000',i
             write(1,'(A24,I1)') 'cp seed.in       run_000',i
             write(1,'(A10,I1)')  'cd run_000',i
             write(1,'(A20)')    './bh_v07 > logfile &'
             close(1)
             call system ('chmod a+x move') 
             call system ('./move') 
           elseif ((i.gt.0).and.(i.lt.100)) then
             write(1,'(A12,I2)') 'mkdir run_00',i
             write(1,'(A23,I2)') 'cp XxYy.in       run_00',i
             write(1,'(A23,I2)') 'cp bh_v07        run_00',i
             write(1,'(A23,I2)') 'cp input_bh.in   run_00',i
             write(1,'(A23,I2)') 'cp input_gru.dat run_00',i
             write(1,'(A23,I2)') 'cp seed.in       run_00',i
             write(1,'(A9,I2)')  'cd run_00',i
             write(1,'(A20)')    './bh_v07 > logfile &'
             close(1)
             call system ('chmod a+x move') 
             call system ('./move') 
           elseif (i.ge.100) then
             write(1,'(A11,I3)') 'mkdir run_0',i
             write(1,'(A22,I3)') 'cp XxYy.in       run_0',i
             write(1,'(A22,I3)') 'cp bh_v07        run_0',i
             write(1,'(A22,I3)') 'cp input_bh.in   run_0',i
             write(1,'(A22,I3)') 'cp input_gru.dat run_0',i
             write(1,'(A22,I3)') 'cp seed.in       run_0',i
             write(1,'(A8,I3)')  'cd run_0',i
             write(1,'(A20)')    './bh_v07 > logfile &'
             close(1)
             call system ('chmod a+x move') 
             call system ('./move') 
           endif
         enddo 

         close(3)

         end program parallel_orbits
