module Generate_Initial_Conditions
    use Dumbbell_util
    implicit none

    contains

    function psiQ_FF(Q, alpha, Q0, Jeq)
        implicit none
        real*8, intent(in) :: Q, alpha, Q0, Jeq
        real*8 :: psiQ_FF

    !    Jeq = (1.D0/(alpha+3.D0)+Q0**2/alpha)*beta(0.5D0,(alpha+2.D0)/2.D0)*alpha**(1.5D0)

        psiQ_FF = Q**2*(1.D0-(Q-Q0)**2.D0/alpha)**(alpha/2.D0)/Jeq
    end function psiQ_FF

    function integral_psiQ_FF(Q, alpha, Q0, Jeq)
        !Simple cumulative trapezoidal integral of psiQ_FF at
        !points specified in Q
        implicit none
        real*8, dimension(:), intent(in) :: Q
        real*8, intent(in) :: alpha, Q0, Jeq
        real*8, dimension(size(Q)) :: integral_psiQ_FF
        integer :: k

        integral_psiQ_FF(1) = 0.D0
        do k=2,size(Q)
            integral_psiQ_FF(k) = integral_psiQ_FF(k-1) + &
                                  (psiQ_FF(Q(k-1),alpha,Q0, Jeq) + psiQ_FF(Q(k),alpha,Q0, Jeq))*(Q(k)-Q(k-1))/2.D0
        end do
        !Trapezoidal rule is far from perfect, but we must have int from 0 to 1
    !    integral_psiQ_FF = integral_psiQ_FF/integral_psiQ_FF(size(Q))

    end function integral_psiQ_FF

    function generate_Ql_eq_FF(N, alpha, Q0, seed, Nsteps)
        implicit none
        real*8, intent(in) :: alpha, Q0
        integer*8, intent(in) :: N, Nsteps
        integer*8, intent(inout) :: seed
        integer :: k
        real*8, dimension(N) :: generate_Ql_eq_FF, rands
        real*8, dimension(Nsteps) :: Q, intpsiQ
        real*8 :: width, Jeq

        !Generate a Q vector between Q0-sqrt(alpha) and Q0+sqrt(alpha) with Nsteps steps
        !Round off a little at the end to prevent singularities
        Q = 0.D0
        if (Q0-sqrt(alpha).le.0.D0) then
            width = (Q0+sqrt(alpha))/(Nsteps-1)
            Q(1) = 0.D0
        else
            width = 2.D0*sqrt(alpha)/(Nsteps-1)
            Q(1) = Q0-sqrt(alpha)
        end if
        do k=2,Nsteps
            Q(k) = Q(k-1) + width
        end do
        Q(1) = Q(1) + 0.0000001D0
        Q(Nsteps) = Q(Nsteps) - 0.0000001D0

    !   determine normalisation constant
    !   Is this needed? It's normalised anyway... TODO
        intpsiQ = integral_psiQ_FF(Q,alpha,Q0, 1.D0)
        Jeq = intpsiQ(size(Q))
    !   perform actual integration and normalise to 1
        intpsiQ = integral_psiQ_FF(Q,alpha,Q0,Jeq)
        Jeq = intpsiQ(size(Q))
        intpsiQ = intpsiQ/Jeq

        rands(:) = rand_floats(seed, N)

        !Inverse of integral returns distribution of psiQ_FF
        do k=1,N
            generate_Ql_eq_FF(k) = inverse_lin_interp(Q,intpsiQ,rands(k))
        end do

    end function generate_Ql_eq_FF

    function inverse_lin_interp(x, fx, fxval)
        implicit none
        real*8, dimension(:), intent(in) :: x, fx
        real*8, intent(in) :: fxval
        real*8 :: inverse_lin_interp
        integer :: i

        do i=1,size(x)
            if (fx(i) > fxval) then
                inverse_lin_interp = (fxval-fx(i-1))/(fx(i)-fx(i-1))*(x(i)-x(i-1)) + x(i-1)
                EXIT
            end if
       end do

    end function inverse_lin_interp

    function generate_Q_FF(Q0,alpha, N, seed, Nsteps)
        implicit none
        real*8, intent(in) :: Q0, alpha
        integer*8, intent(inout) :: seed
        integer*8, intent(in) :: N, Nsteps
        real*8, dimension(3,N) :: generate_Q_FF
        real*8 :: Ql(N)

        generate_Q_FF(:,:) = spherical_unit_vectors(N, seed)

        Ql(:) = generate_Ql_eq_FF(N, alpha, Q0, seed, Nsteps)

        generate_Q_FF(1,:) = generate_Q_FF(1,:)*Ql
        generate_Q_FF(2,:) = generate_Q_FF(2,:)*Ql
        generate_Q_FF(3,:) = generate_Q_FF(3,:)*Ql

    end function generate_Q_FF

    function spherical_unit_vectors(N, seed)

        implicit none
        integer*8, intent(inout) :: seed
        integer*8, intent(in) :: N
        integer*8 :: i
        real*8 :: x1, x2, r1(1), r2(1)
        real*8, dimension(3,N) :: spherical_unit_vectors

        do i=1,N
            do while (.True.)
                r1 = rand_floats(seed,1)
                r2 = rand_floats(seed,1)
                x1 = r1(1)*2.D0 - 1.D0
                x2 = r2(1)*2.D0 - 1.D0
                if ((x1**2 + x2**2).lt.1.D0) EXIT
            end do
            spherical_unit_vectors(1,i) = 2.D0*x1*sqrt(1.D0-x1**2-x2**2)
            spherical_unit_vectors(2,i) = 2.D0*x2*sqrt(1.D0-x1**2-x2**2)
            spherical_unit_vectors(3,i) = 1.D0-2.D0*(x1**2 + x2**2)

        end do

    end function spherical_unit_vectors
end module
