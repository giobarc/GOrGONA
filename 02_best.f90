        program orbits 
        implicit none 

        integer::ncurv,i,j,n1,m1,m2,m3,nmax,p10,p20,p30
        integer,dimension(1000,1000)::p1,p2,p3
        integer,dimension(1000)::p1min,p2min,p3min
        real(kind=8):: e1,d1,ene0,des10 
        real(kind=8),dimension(1000,1000):: ene,des1 
        real(kind=8),dimension(1000):: enemin,des1min 
        character(len=200),dimension(1000):: nfile
        character(len=10),dimension(1000,1000):: nf
        character(len=10),dimension(1000):: nfmin
        character(len=200):: nome
        character(len=10):: nomef,nf0

        open (2,file='result.out',status='replace')
        do i=1,1000
          do j=1,100
            ene(i,j)=100.d0
          enddo
        enddo

        call system ('rm lista')
        call system ('ls curve-0* >> lista')

        ncurv=0
        open (1,file='lista',status='old')
        rewind (unit=1)
        do i=1,10000
          read (1,*,end=3) nfile(i)
          ncurv=ncurv+1
        enddo
  3     continue
        close (1)

        write(*,*) ncurv,'curve'
        do i=1,ncurv
          write(*,*) nfile(i)
        enddo

        do i=1,ncurv
          nome=nfile(i)
          open (1,file=nome,status='old')
          do j=1,10000
            read(1,*,end=4) n1,e1,d1,m1,m2,m3,nomef
            if (n1.eq.0) then
              ene0=e1
              des10=d1
              p10=m1
              p20=m2
              p30=m3
              nf0=nomef
            else
              ene(n1,i)=e1
              des1(n1,i)=d1
              p1(n1,i)=m1
              p2(n1,i)=m2
              p3(n1,i)=m3
              nf(n1,i)=nomef
              nmax=n1
            endif
          enddo
   4      continue
          close (1)
        enddo
        
        write(*,*) nmax
         
        write(2,10) 0,ene0,des10,p10,p20,p30,nf0

        do i=1,nmax
          enemin(i)=80.d0
          do j=1,ncurv
            if (ene(i,j).le.enemin(i)) then
              enemin(i)=ene(i,j)
              des1min(i)=des1(i,j)
              p1min(i)=p1(i,j)
              p2min(i)=p2(i,j)
              p3min(i)=p3(i,j)
              nfmin(i)=nf(i,j)
            endif
          enddo
          if (enemin(i).le.70.d0) then
           write(2,10) i,enemin(i),des1min(i),p1min(i),p2min(i),p3min(i)&
     &,nfmin(i)
          endif
        enddo

        close(2)
        
  10    format(I6,2F10.4,3I6,4X,A10)
 
        end program orbits
