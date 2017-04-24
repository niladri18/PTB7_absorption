module make_huckel

        use parameter
contains

subroutine make_huckel_level

! ========= Constructing Huckel Matrix =======

        arr = 0.0

        do i=1,n
                arr(i,i)=eps(i)
        enddo

        do i=1,nbonds1
                i1=ilbond(1,i)
                i2=ilbond(2,i)
                arr(i1,i2)=t1
                arr(i2,i1)=t1
                !print*,i1,i2,arr(i1,i2),arr(i2,i1)
        enddo

        do i=nbonds1+1,nbonds1+nbonds2
                i1=ilbond(1,i)
                i2=ilbond(2,i)
                arr(i1,i2)=t2
                arr(i2,i1)=t2
                !print*,i1,i2,arr(i1,i2),arr(i2,i1)
        enddo


        do i=nbonds1+nbonds2+1,nbonds1+nbonds2+nbonds3
                i1=ilbond(1,i)
                i2=ilbond(2,i)
                arr(i1,i2)=t3
                arr(i2,i1)=t3
        !print*,i1,i2,arr(i1,i2),arr(i2,i1)
        enddo

        do i=nbonds1+nbonds2+nbonds3+1,nbonds1+nbonds2+nbonds3+nbonds4
                i1=ilbond(1,i)
                i2=ilbond(2,i)
                arr(i1,i2)=tb
                arr(i2,i1)=tb
                !print*,i1,i2,arr(i1,i2),arr(i2,i1)
        enddo


        print*,"ended constructing huckel matrix"


        f = 0.0
        v_ij = 0.0

        f0=arr

        call DSYEV('V','U', N, arr, N, tmp, WORK, lwork, info )

        print*,"HUCKEL ENERGY LEVELS"

        do i=1,n
                print*,i,tmp(i)
        end do

        print*,"Huckel calculation ends"

end subroutine make_huckel_level


end module make_huckel
