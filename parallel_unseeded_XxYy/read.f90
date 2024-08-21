        program orbits 
        implicit none 

        integer::i,ns,ntot,m1,m2,m3,nat
        integer,dimension(1000)::n1
        real(kind=8),dimension(1000):: exc,ene,parm
        real(kind=8):: emet1,emet2
        character(len=20), dimension(1000)::a1,a2,a3,a4,a5,a6,a7,a21,a31
        character(len=10), dimension(1000):: nome 

        open (1,file='res',status='old')
        open (2,file='res1',status='replace')
        open (4,file='param.dat',status='old')
        open (3,file='curve.dat',status='replace')

        read(4,*) m1,m2,m3
        ns=0
        do i=1,10000
          read (1,10,end=3) a1(i),a2(i),a21(i),a3(i),a31(i),a4(i),a5(i),&
     &a6(i),a7(i) 
          ns=ns+1
        enddo 
   3    continue

        do i=1,ns
          write(2,30) a2(i),a4(i),a7(i),a3(i)
        enddo
        rewind(2)

        do i=1,ns
          read(2,*) n1(i),ene(i),parm(i),nome(i)
        enddo

        ntot=n1(ns)
        write(*,*) ntot

        do i=1,ns
          if (n1(i).eq.0) then
            emet2=ene(i)/ntot
          endif
          if (n1(i).eq.ntot) then
            emet1=ene(i)/ntot
          endif
        enddo
        write(*,*) emet2,emet1 

        do i=1,ns
          exc(i)=ene(i)-n1(i)*emet1-(ntot-n1(i))*emet2
          write(3,20) n1(i),exc(i),parm(i),m1,m2,m3,nome(i)
        enddo

 10     format(A6,A4,A1,A10,A6,3A16,A10)
 30     format(A4,A16,A10,4X,A10)
 20     format(I6,2F10.4,3I6,4X,A10)

        close (1)
        close (2)
        close (3)
        close (4)

        end program orbits
