module Dumbbell_util
    implicit none
    integer, parameter :: quad = SELECTED_REAL_KIND(32)
    real*8, parameter :: PI = 4.D0*atan(1.0D0)
    real*8, parameter, dimension(3,3) :: delT = reshape((/1, 0, 0, &
                                                          0, 1, 0, &
                                                          0, 0, 1/), &
                                                        shape(delT))

    type measured_variables
        real*8 :: Qavg, Vqavg, S, Serr, Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2
    end type

    contains

    subroutine measure_shear_no_VR(meas, Q, Q0, alpha, sr, Ntraj)
        implicit none
        type(measured_variables), intent(out) :: meas
        type(measured_variables) :: meas_tmp
        integer*8, intent(in) :: Ntraj
        integer :: i
        real*8, intent(in) :: Q0, alpha, Q(3,Ntraj), sr
        real*8 :: tau(3,3), F(3), Bs, Ql, Ql2, B_eta, Bpsi, Bpsi2

        ! These variables are all global and shared between threads
        meas%Aeta = 0.D0; meas%Apsi = 0.D0; meas%Apsi2 = 0.D0
        meas%Veta = 0.D0; meas%Vpsi = 0.D0; meas%Vpsi2 = 0.D0
        meas%Qavg = 0.D0; meas%Vqavg = 0.D0
        meas%S = 0.D0; meas%Serr = 0.D0

        !These variables SHOULD(!) be private
        tau = 0.D0
        meas_tmp%Aeta = 0.D0; meas_tmp%Apsi = 0.D0; meas_tmp%Apsi2 = 0.D0
        meas_tmp%Veta = 0.D0; meas_tmp%Vpsi = 0.D0; meas_tmp%Vpsi2 = 0.D0
        meas_tmp%Qavg = 0.D0; meas_tmp%Vqavg = 0.D0
        meas_tmp%S = 0.D0; meas_tmp%Serr = 0.D0

        !$OMP DO
        do i=1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product( (/1,0,0/), Q(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            meas_tmp%Qavg = meas_tmp%Qavg + Ql2
            meas_tmp%Vqavg = meas_tmp%Vqavg + Ql
            meas_tmp%S = meas_tmp%S + 0.5*(3*Bs - 1)
            meas_tmp%Serr = meas_tmp%Serr + 0.25*(9*Bs**2 - 6*Bs + 1)

            B_eta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            meas_tmp%Aeta = meas_tmp%Aeta + B_eta
            meas_tmp%Apsi = meas_tmp%Apsi + Bpsi
            meas_tmp%Apsi2 = meas_tmp%Apsi2 + Bpsi2
            meas_tmp%Veta = meas_tmp%Veta + B_eta**2
            meas_tmp%Vpsi = meas_tmp%Vpsi + Bpsi**2
            meas_tmp%Vpsi2 = meas_tmp%Vpsi2 + Bpsi2**2
        end do
        !$OMP END DO

        !$OMP atomic
        meas%Aeta = meas%Aeta + meas_tmp%Aeta
        !$OMP atomic
        meas%Veta = meas%Veta + meas_tmp%Veta
        !$OMP atomic
        meas%Apsi = meas%Apsi + meas_tmp%Apsi
        !$OMP atomic
        meas%Vpsi = meas%Vpsi + meas_tmp%Vpsi
        !$OMP atomic
        meas%Apsi2 = meas%Apsi2 + meas_tmp%Apsi2
        !$OMP atomic
        meas%Vpsi2 = meas%Vpsi2 + meas_tmp%Vpsi2
        !$OMP atomic
        meas%Qavg = meas%Qavg + meas_tmp%Qavg
        !$OMP atomic
        meas%Vqavg = meas%Vqavg + meas_tmp%Vqavg
        !$OMP atomic
        meas%S = meas%S + meas_tmp%S
        !$OMP atomic
        meas%Serr = meas%Serr + meas_tmp%Serr

        !$OMP barrier

        !$OMP single
        meas%Aeta = meas%Aeta/(Ntraj*sr)
        meas%Veta = meas%Veta/(Ntraj*sr**2)
        meas%Veta = sqrt((meas%Veta - meas%Aeta**2)/(Ntraj-1))

        meas%Apsi = meas%Apsi/(Ntraj*sr**2)
        meas%Vpsi = meas%Vpsi/(Ntraj*sr**4)
        meas%Vpsi = sqrt((meas%Vpsi - meas%Apsi**2)/(Ntraj-1))

        meas%Apsi2 = meas%Apsi2/(Ntraj*sr**2)
        meas%Vpsi2 = meas%Vpsi2/(Ntraj*sr**4)
        meas%Vpsi2 = sqrt((meas%Vpsi2 - meas%Apsi2**2)/(Ntraj-1))

        meas%Qavg = sqrt(meas%Qavg/Ntraj)
        meas%Vqavg = meas%Vqavg/Ntraj
        meas%Vqavg = sqrt((meas%Qavg**2 - meas%Vqavg**2)/(Ntraj-1))

        meas%S = meas%S/Ntraj
        meas%Serr = meas%Serr/Ntraj
        meas%Serr = sqrt((meas%Serr - meas%S**2)/(Ntraj-1))
        !$OMP end single

    end subroutine

    subroutine measure_shear_with_VR(meas, Q, Q_eq_VR, Q0, alpha, sr, Ntraj)
        implicit none
        type(measured_variables), intent(out) :: meas
        type(measured_variables) :: meas_tmp
        integer*8, intent(in) :: Ntraj
        integer :: i
        real*8, intent(in) :: Q0, alpha, Q(3,Ntraj), Q_eq_VR(3,Ntraj), sr
        real*8 :: tau(3,3), F(3), Bs, Ql, Ql2, B_eta, Bpsi, Bpsi2

        ! These variables are all global and shared between threads
        meas%Aeta = 0.D0; meas%Apsi = 0.D0; meas%Apsi2 = 0.D0
        meas%Veta = 0.D0; meas%Vpsi = 0.D0; meas%Vpsi2 = 0.D0
        meas%Qavg = 0.D0; meas%Vqavg = 0.D0
        meas%S = 0.D0; meas%Serr = 0.D0

        !These variables SHOULD(!) be private
        tau = 0.D0
        meas_tmp%Aeta = 0.D0; meas_tmp%Apsi = 0.D0; meas_tmp%Apsi2 = 0.D0
        meas_tmp%Veta = 0.D0; meas_tmp%Vpsi = 0.D0; meas_tmp%Vpsi2 = 0.D0
        meas_tmp%Qavg = 0.D0; meas_tmp%Vqavg = 0.D0
        meas_tmp%S = 0.D0; meas_tmp%Serr = 0.D0

        !$OMP DO
        do i=1,Ntraj
            !Calculated shear-flow values
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product( (/1,0,0/), Q(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            B_eta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))

            meas_tmp%Qavg = meas_tmp%Qavg + Ql2
            meas_tmp%Vqavg = meas_tmp%Vqavg + Ql
            meas_tmp%S = meas_tmp%S + 0.5*(3*Bs - 1)
            meas_tmp%Serr = meas_tmp%Serr + 0.25*(9*Bs**2 - 6*Bs + 1)

            !subtract equilibrium values from shear-flow values
            Ql2 = Q_eq_VR(1,i)**2 + Q_eq_VR(2,i)**2 + Q_eq_VR(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product( (/1,0,0/), Q_eq_VR(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q_eq_VR(:,i)/Ql
            tau(:,:) = dyadic_prod(Q_eq_VR(:,i), F)

            B_eta = B_eta - tau(1,2)
            Bpsi = Bpsi - (tau(1,1) - tau(2,2))
            Bpsi2 = Bpsi2 - (tau(2,2) - tau(3,3))

            meas_tmp%S = meas_tmp%S - 0.5*(3*Bs - 1)
            meas_tmp%Serr = meas_tmp%Serr - 0.25*(9*Bs**2 - 6*Bs + 1)
            meas_tmp%Aeta = meas_tmp%Aeta + B_eta
            meas_tmp%Apsi = meas_tmp%Apsi + Bpsi
            meas_tmp%Apsi2 = meas_tmp%Apsi2 + Bpsi2
            meas_tmp%Veta = meas_tmp%Veta + B_eta**2
            meas_tmp%Vpsi = meas_tmp%Vpsi + Bpsi**2
            meas_tmp%Vpsi2 = meas_tmp%Vpsi2 + Bpsi2**2
        end do
        !$OMP END DO

        !$OMP atomic
        meas%Aeta = meas%Aeta + meas_tmp%Aeta
        !$OMP atomic
        meas%Veta = meas%Veta + meas_tmp%Veta
        !$OMP atomic
        meas%Apsi = meas%Apsi + meas_tmp%Apsi
        !$OMP atomic
        meas%Vpsi = meas%Vpsi + meas_tmp%Vpsi
        !$OMP atomic
        meas%Apsi2 = meas%Apsi2 + meas_tmp%Apsi2
        !$OMP atomic
        meas%Vpsi2 = meas%Vpsi2 + meas_tmp%Vpsi2
        !$OMP atomic
        meas%Qavg = meas%Qavg + meas_tmp%Qavg
        !$OMP atomic
        meas%Vqavg = meas%Vqavg + meas_tmp%Vqavg
        !$OMP atomic
        meas%S = meas%S + meas_tmp%S
        !$OMP atomic
        meas%Serr = meas%Serr + meas_tmp%Serr

        !$OMP barrier

        !$OMP single
        meas%Aeta = meas%Aeta/(Ntraj*sr)
        meas%Veta = meas%Veta/(Ntraj*sr**2)
        meas%Veta = sqrt((meas%Veta - meas%Aeta**2)/(Ntraj-1))

        meas%Apsi = meas%Apsi/(Ntraj*sr**2)
        meas%Vpsi = meas%Vpsi/(Ntraj*sr**4)
        meas%Vpsi = sqrt((meas%Vpsi - meas%Apsi**2)/(Ntraj-1))

        meas%Apsi2 = meas%Apsi2/(Ntraj*sr**2)
        meas%Vpsi2 = meas%Vpsi2/(Ntraj*sr**4)
        meas%Vpsi2 = sqrt((meas%Vpsi2 - meas%Apsi2**2)/(Ntraj-1))

        meas%Qavg = sqrt(meas%Qavg/Ntraj)
        meas%Vqavg = meas%Vqavg/Ntraj
        meas%Vqavg = sqrt((meas%Qavg**2 - meas%Vqavg**2)/(Ntraj-1))

        meas%S = meas%S/Ntraj
        meas%Serr = meas%Serr/Ntraj
        meas%Serr = sqrt((meas%Serr - meas%S**2)/(Ntraj-1))
        !$OMP end single

    end subroutine

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

    pure function find_roots(a, b, c, lower_bound, upper_bound)
        implicit none
        real*8, intent(in) :: a, b, c, lower_bound, upper_bound
        real*8 :: find_roots
        real*8 :: Q, R, theta, x, Au, Bu
        integer :: i, iter

        Q = (a**2 - 3.D0*b)/9.D0
        R = (2.D0*a**3 - 9.D0*a*b + 27.D0*c)/54.D0

        !If roots are all real, get root in range
        if (R**2.lt.Q**3) then
            iter = 0
            theta = acos(R/sqrt(Q**3))
            do i=-1,1
                iter = iter+1
                x = -2.D0*sqrt(Q)*cos((theta + dble(i)*PI*2.D0)/3.D0)-a/3.D0
                if ((x.ge.lower_bound).and.(x.le.upper_bound)) then
                    find_roots = x
                    return
                end if
            end do
        !Otherwise, two imaginary roots and one real root, return real root
        else
            Au = -sign(1.D0, R)*(abs(R)+sqrt(R**2-Q**3))**(1.D0/3.D0)
            if (Au.eq.0.D0) then
                Bu = 0.D0
            else
                Bu = Q/Au
            end if
            find_roots = (Au+Bu)-a/3.D0
            return
        end if

    end function find_roots

    pure function find_roots_quad(a, b, c, lower_bound, upper_bound)
        implicit none
        real(quad), intent(in) :: a, b, c, lower_bound, upper_bound
        real*8 :: find_roots_quad
        real(quad) :: Q, R, theta, x
        integer :: i

        Q = (a**2 - 3.Q0*b)/9.Q0
        R = (2.Q0*a**3 - 9.Q0*a*b + 27.Q0*c)/54.Q0

        theta = acos(R/sqrt(Q**3))
        do i=-1,1
            x = -2.Q0*sqrt(Q)*cos((theta + real(i, kind=quad)*PI*2.Q0)/3.Q0)-a/3.Q0
            if ((x.ge.lower_bound).and.(x.le.upper_bound)) then
                find_roots_quad = dble(x)
                return
            end if
        end do

    end function

    function generate_lookup_table(beta, Q0, alpha, Nsteps)
        implicit none
        real(quad), intent(in) :: beta, Q0, alpha
        integer, intent(in) :: Nsteps
        integer :: i
        real*8, dimension(:, :), allocatable :: generate_lookup_table
        real(quad) :: step
        real(quad), dimension(:), allocatable :: Y

        allocate(Y(1:2*Nsteps), generate_lookup_table(1:2*Nsteps, 2))

        step = 20.Q0*sqrt(alpha)/(real(Nsteps,kind=quad)-1.Q0)

        Y(1:Nsteps) = (/((Q0-10.Q0*sqrt(alpha)+(real(i,kind=quad)-1.Q0)*step),i=1,Nsteps)/)

        step = (log10(10000.Q0*sqrt(alpha)+10.Q0*Q0)-log10(10.Q0*sqrt(alpha)+Q0+step))/(real(Nsteps,kind=quad)-1.Q0)

        Y(Nsteps+1:2*Nsteps) = (/((log10(10.Q0*sqrt(alpha)+Q0+step)+(real(i,kind=quad)-1.Q0)*step),i=1,Nsteps)/)
        Y(Nsteps+1:2*Nsteps) = 10**(Y(Nsteps+1:2*Nsteps))

        do i=1,size(Y)
            generate_lookup_table(i, 2) = find_roots_quad(-(2.Q0*Q0+Y(i)), &
                               -(alpha-Q0**2-2.Q0*Y(i)*Q0+beta), &
                               (beta*Q0+Y(i)*alpha-Y(i)*Q0**2), &
                               Q0-sqrt(alpha), Q0+sqrt(alpha))
        end do

        generate_lookup_table(:, 1) = dble(Y)

    end function

    pure function locate(xvals, x, m)
        ! Given a value x, return a value j such that x is (insofar as possible) centered in the subrange xvals(j..j+m-1)
        ! The values in xvals must be monotonic ascending. The returned value is not less than 0, nor greater than n-1.
        implicit none
        integer :: locate
        real*8, intent(in) :: xvals(:), x
        integer, intent(in) :: m
        integer :: ju, jm, jl, n


        n = size(xvals)
        jl = 1
        ju = n

        do while(ju-jl.gt.1)
            jm = (ju+jl)/2
            if (x.ge.xvals(jm)) then
                jl = jm
            else
                ju=jm
            end if
        end do

        locate = max(1, min(n-m, (jl-((m-2)/2))))

    end function

    function hunt(xvals, x, m, jguess)
        implicit none
        integer :: hunt
        real*8, intent(in) :: xvals(:), x
        integer, intent(in) :: m, jguess
        integer :: ju, jm, jl, n, inc
        ! Given a value x, return a value j such that x is (insofar as possible) centered in the subrange xvals(j..j+m-1)
        ! The values in xvals must be monotonic ascending. The returned value is not less than 0, nor greater than n-1.

        n = size(xvals)
        inc = 1
        jl = jguess
        !Input guess not useful, go to bisection
        if ((jl.lt.0).or.(jl.gt.n-1)) then
            jl=0
            ju=n-1
        else
            !hunt up
            if (x.ge.xvals(jl)) then
                do while(.true.)
                    ju = jl + inc
                    if (ju.ge.n-1) then
                        ju = n-1
                        EXIT
                    elseif (x.lt.xvals(ju)) then
                        exit
                    else
                        jl = ju
                        inc = inc + 1
                    end if
                end do
            !hunt down
            else
                do while(.true.)
                    ju = jl - inc
                    if (ju.le.0) then
                        jl = 0
                        EXIT
                    elseif (x.ge.xvals(jl)) then
                        exit
                    else
                        ju = jl
                        inc = inc + 1
                    end if
                end do
            end if
        end if

        !do bisection for final answer
        do while(ju-jl > 1)
            jm = (ju+jl)/2
            if (x.ge.xvals(jm)) then
                jl = jm
            else
                ju=jm
            end if
        end do

        hunt = max(0, min(n-m, (jl-((m-2)/2))))

    end function

    pure function lin_interp_bs(xvals, yvals, x)
        implicit none
        real*8 :: lin_interp_bs
        real*8, intent(in) :: xvals(:), yvals(:), x
        integer :: j

        j = locate(xvals, x, 2)
        lin_interp_bs = yvals(j) + ((x-xvals(j))/(xvals(j+1)-xvals(j)))*(yvals(j+1)-yvals(j))
    end function

    function poly_interp_bs(xvals, yvals, x, mm)
        implicit none
        real*8 :: poly_interp_bs
        real*8, intent(in) :: xvals(:), yvals(:), x
        integer, intent(in) :: mm
        integer :: j, i, m, ns
        real*8 :: xa(mm), ya(mm), c(mm), d(mm), ho, hp, den, w, dy, dif, dift

        j = locate(xvals, x, mm)

        do i=1,mm
            xa(i) = xvals(j+i-1)
            ya(i) = yvals(j+i-1)
        end do

        !find index ns of closest table entry, initialise values
        dif = abs(x-xa(1))
        do i=1,mm
            dift = abs(x-xa(i))
            if (dift.le.dif) then
                ns = i
                dif = dift
            end if
            c(i) = ya(i)
            d(i) = ya(i)
        end do

        !initial guess
        poly_interp_bs = ya(ns)
        ns = ns - 1

        do m=1,mm-1
            do i=1,mm-m
                ho = xa(i) - x
                hp = xa(i+m) - x
                w = c(i+1) - d(i)
                den = ho-hp
                den = w/den
                d(i) = hp*den
                c(i) = ho*den
            end do
            dy = 2.D0*dble(ns-1)+1.D0
            if (dy.lt.(mm-m)) then
                poly_interp_bs = poly_interp_bs + c(ns+1)
            else
                poly_interp_bs = poly_interp_bs + d(ns)
                ns = ns - 1
            end if
        end do
    end function

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

    pure function dyadic_prod(vec1,vec2)
        implicit none
        real*8, dimension(3), intent(in) :: vec1, vec2
        real*8, dimension(3,3) :: dyadic_prod
        integer :: i, j

        do i=1,3
            do j=1,3
                dyadic_prod(i,j) = vec1(i)*vec2(j)
            end do
        end do
    end function dyadic_prod

    pure function ten_vec_dot(tensor, vector)
        implicit none
        real*8, dimension(3,3), intent(in) :: tensor
        real*8, dimension(3), intent(in) :: vector
        real*8, dimension(3) :: ten_vec_dot
        integer :: i

        do i=1,3
            ten_vec_dot(i) = dot_product(tensor(i,:), vector)
        end do
    end function ten_vec_dot

end module
