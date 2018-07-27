module Dumbbell_util

    implicit none

    contains

    function shift_xor(val,shift)
        integer*8 :: shift_xor
        integer*8, intent(in) :: val, shift
        shift_xor = ieor(val,ishft(val,shift))
    end function

    function rand_floats(seed, N)
        implicit none
        integer*8, intent(in) :: N
        integer*8, intent(inout) :: seed
        real*8, dimension(N) :: rand_floats
        integer :: i

        do i=1,N
            !Generates a random number between 0 and 1
            !Using xorshift and one round of 64-bit MCG
            seed = shift_xor(shift_xor(shift_xor(seed, 13_8),-17_8),43_8)
            rand_floats(i) = seed * 2685821657736338717_8 * 5.42101086242752217D-20 + 0.5D0
        end do

    end function rand_floats

    function Wiener_step(seed, dt)
        implicit none
        integer*8, intent(inout) :: seed
        real*8, intent(in) :: dt
        real*8, dimension(3) :: Wiener_step
        real*8, dimension(3) :: dW
        integer :: i

        do i=1,3
            !Generates a random number between -0.5 and 0.5
            !Using xorshift and one round of 64-bit MCG
            seed = shift_xor(shift_xor(shift_xor(seed, 13_8),-17_8),43_8)
            dW(i) = seed * 2685821657736338717_8 * 5.42101086242752217D-20
        end do

        !Generates an approximately gaussian-distributed number dW
        Wiener_step = dW*sqrt(dt)*(14.14855378D0*dW*dW + 1.21569221D0)
        return

    end function Wiener_step

    pure function beta(x,y)
        implicit none
        real*8, intent(in) :: x, y
        real*8 :: beta

        beta = log_gamma(x) + log_gamma(y) - log_gamma(x+y)
        beta = exp(beta)
    end function beta

    function betai(a,b,x)
        !From numerical recipes, uses betacf
        !Outputs the incomplete (and regularised) beta function
        implicit none
        real*8, intent(in) :: a, b, x
        real*8 :: betai
        real*8 :: bt

        if(x.eq.0.D0.or.x.eq.1.D0) then
            bt = 0.D0
        else
            bt = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) &
                     + a*log(x) + b*log(1.D0 - x))
        end if

        if(x.lt.(a+1.D0)/(a+b+2.D0)) then
            betai = bt*betacf(a,b,x)/a
            return
        else
            betai = 1.D0 - bt*betacf(b, a, 1.D0-x)/b
            return
        end if

    end function

    function betacf(a,b,x)
        real*8, intent(in) :: a, b, x
        real*8 :: betacf
        integer*8, parameter :: maxit = 100
        real*8, parameter :: eps = 3.D-7
        real*8, parameter :: fpmin = 1.D-30
        integer*8 :: m, m2
        real*8 :: aa, c, d, del, h, qab, qam, qap

        qab = a+b
        qap = a+1.D0
        qam = a-1.D0
        c = 1.D0
        d = 1.D0 - qab*x/qap
        if (abs(d).lt.fpmin) d=fpmin
        d = 1.D0/d
        h = d
        do m = 1,maxit
            m2 = 2.D0*m
            ! step one of the recurrence
            aa = m*(b-m)*x/((qam+m2)*(a+m2))
            d = 1.D0+aa*d
            if (abs(d).lt.fpmin) d=fpmin
            c = 1.D0+aa/c
            if (abs(c).lt.fpmin) c=fpmin
            d = 1.D0/d
            h = h*d*c
            ! step two (odd step) of recurrence
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
            d = 1.D0+aa*d
            if (abs(d).lt.fpmin) d=fpmin
            c = 1.D0+aa/c
            if (abs(c).lt.fpmin) c=fpmin
            d = 1.D0/d

            del = d*c
            h = h*del
            if (abs(del-1.D0).lt.eps) EXIT
            if (m.eq.maxit) print *, "reached max iterations betacf"
        end do
        betacf = h

    end function

end module
