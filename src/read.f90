module readfile

        use parameter

contains

subroutine read_config
! ====== Reading Main File ====

   !open(unit=90,file="./input/ptb7cord4.inp")
   !open(unit=91,file="./input/unit4.inp")


         !read(91,*) n
         read(91,*) nbonds1,nbonds2,nbonds3,nbonds4!,nbonds5 
         nbnd=nbonds1+nbonds2+nbonds3+nbonds4!+nbonds5
         print*,"total no of bonds",nbnd
         !nbonds1 = C-S; nbonds2 = C-C; nbonds3 = C=C; nbonds4 = benzene bonds; nbonds5
         != interchain bonds
         read(91,*) t1,t2,t3,tb,kappa,eta
        !write(*,*) 't1,t2,t3',t1,t2,t3
        !Define this array variable: ilbnd(2,100)
        read(91,*) (ilbond(1,k),ilbond(2,k),k=1,nbnd)
        !write(*,*) 'links',(ilbond(1,k),ilbond(2,k),k=1,nbnd)
        !Define siten, Hub-U, and chem-pot arrays
         read(91,*) (eps(i),i=1,n)
         read(91,*) (U_c(i),i=1,n)
         read(91,*) (z(i),i=1,n)
         !write(*,*) 'z',(z(i),i=1,n)


         print*,"HOMO",n/2+sulp/2,"LUMO",n/2+sulp/2+1
        ! do i=1,n
        ! print*,x_pos(i),z(i),U(i)
        ! end do


 print*,"ended reading coordinates"

end subroutine read_config


end module readfile
