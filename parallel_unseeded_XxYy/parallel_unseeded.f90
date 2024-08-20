         program parallel_orbits 
         implicit none
         integer:: i,ni,nf,ns,nat

         open(1,file='size.dat',status='old')
         read(1,*) nat
         close(1)

         open(3,file='range.dat',status='old')
         read(3,*) ni,nf,ns
         do i=ni,nf,ns
           open(1,file='noseed_rgl.in',status='replace')
           write(1,'(I6)') nat 
           write(1,'(A5)') 'Xx Yy' 
           write(1,'(2I6)') i,nat-i 
           write(1,'(A20)') '!cartesian ref_frame' 
           write(1,'(A22)') '-6.0, 6.0 !x_min,x_max' 
           write(1,'(A22)') '-6.0, 6.0 !y_min,y_max' 
           write(1,'(A22)') '-6.0, 6.0 !z_min,z_max' 
           close(1)
           open(1,file='move',status='replace')
           if (i.eq.0) then
             write(1,'(A14)') 'mkdir run_0000'
             write(1,'(A25)') 'cp XxYy.in       run_0000'
             write(1,'(A25)') 'cp bh_v07        run_0000'
             write(1,'(A37)') 'cp input_bh.in_P run_0000/input_bh.in'
             write(1,'(A25)') 'cp input_gru.dat run_0000'
             write(1,'(A25)') 'cp noseed_rgl.in run_0000'
             write(1,'(A11)') 'cd run_0000'
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
             write(1,'(A23,I2)') 'cp noseed_rgl.in run_00',i
             write(1,'(A9,I2)')  'cd run_00',i
             write(1,'(A20)')    './bh_v07 > logfile &'
             close(1)
             call system ('chmod a+x move') 
             call system ('./move') 
           elseif ((i.ge.100).and.(i.lt.1000).and.(i.ne.nat)) then 
             write(1,'(A11,I3)') 'mkdir run_0',i
             write(1,'(A22,I3)') 'cp XxYy.in       run_0',i
             write(1,'(A22,I3)') 'cp bh_v07        run_0',i
             write(1,'(A22,I3)') 'cp input_bh.in   run_0',i
             write(1,'(A22,I3)') 'cp input_gru.dat run_0',i
             write(1,'(A22,I3)') 'cp noseed_rgl.in run_0',i
             write(1,'(A8,I3)')  'cd run_0',i
             write(1,'(A20)')    './bh_v07 > logfile &'
             close(1)
             call system ('chmod a+x move') 
             call system ('./move') 
           elseif (i.eq.nat) then
             write(1,'(A11,I3)') 'mkdir run_0',i
             write(1,'(A22,I3)') 'cp XxYy.in       run_0',i
             write(1,'(A22,I3)') 'cp bh_v07        run_0',i
       write(1,'(A22,I3,A12)') 'cp input_bh.in_P run_0',i,'/input_bh.in'
             write(1,'(A22,I3)') 'cp input_gru.dat run_0',i
             write(1,'(A22,I3)') 'cp noseed_rgl.in run_0',i
             write(1,'(A8,I3)')  'cd run_0',i
             write(1,'(A20)')    './bh_v07 > logfile &'
             close(1)
             call system ('chmod a+x move') 
             call system ('./move') 
           endif
         enddo 

         close(3)

         end program parallel_orbits
