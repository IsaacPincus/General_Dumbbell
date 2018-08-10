module Integrate
    use Dumbbell_util
    implicit none
    real*8, parameter :: PI = 4.D0*atan(1.0D0)
    integer, parameter :: quad = SELECTED_REAL_KIND(32)
    real*8, parameter, dimension(3,3) :: delT = reshape((/1, 0, 0, &
                                                          0, 1, 0, &
                                                          0, 0, 1/), &
                                                        shape(delT))

    interface step
        module procedure step_semimp_FF
        module procedure step_semimp_FF_old
        module procedure step_Euler_Hookean
    end interface

    contains

    function construct_B_ROB(Q, a)
        implicit none
        real*8, intent(in) :: Q(3), a
        real*8, dimension(3,3) :: construct_B_ROB
        real*8 :: Ql, aux, g, g_til, Ql2, Ql4, Ql6, a2, a4
        real*8, save :: C43 = 4.D0/3.D0
        real*8, save :: C83 = 8.D0/3.D0
        real*8, save :: C143 = 14.D0/3.D0

        a2 = a**2
        a4 = a2**2

        Ql2 = Q(1)**2 + Q(2)**2 + Q(3)**2
        Ql = sqrt(Ql2)
        Ql4 = Ql2**2
        Ql6 = Ql2**3

        aux = a/(C43*Ql*(Ql2 + C43*a2)**3)
        g = 1.D0-aux*(Ql6 + C143*a2*Ql4 + 8.D0*a4*Ql2)
        g_til = -aux*(Ql6 + 2.D0*a2*Ql4 - C83*a4*Ql2)
        construct_B_ROB = sqrt(g)*delT + (sqrt(g+g_til)-sqrt(g))*dyadic_prod(Q,Q)/Ql2
        return

    end function construct_B_ROB

    function construct_B_RPY(Q, a)
        implicit none
        real*8, intent(in) :: Q(3), a
        real*8, dimension(3,3) :: construct_B_RPY
        real*8 :: Ql, aux, g, g_til, Ql2, Ql4, Ql6, a2, Atmp, Btmp, temp
        real*8, save :: C43 = 4.D0/3.D0
        real*8, save :: C38 = 3.D0/8.D0
        real*8, save :: C18 = 1.D0/8.D0
        real*8, save :: C23 = 2.D0/3.D0

        a2 = a**2

        Ql2 = Q(1)**2 + Q(2)**2 + Q(3)**2
        Ql = sqrt(Ql2)

        if (Ql.ge.(2*a)) then
            temp = a2/Ql2
            Atmp = 1 + C23*temp
            Btmp = 1 - 2*temp
        elseif (Ql.lt.(2*a)) then
            temp = Ql2/a2
            Atmp = C43*Ql/a-C38*temp
            Btmp = C18*temp
        end if

        aux = a/(C43*Ql)
        g = 1.D0-aux*Atmp
        g_til = -aux*Btmp
        construct_B_RPY = sqrt(g)*delT + (sqrt(g+g_til)-sqrt(g))*dyadic_prod(Q,Q)/Ql2
    end function construct_B_RPY

    function step_Euler_Hookean(Q, k, dt, a, dW)
        implicit none
        real*8, intent(in) :: Q(3), dt, a, dW(3)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql
        real*8, dimension(3) :: F
        real*8, dimension(3,3) :: B, BdotB
        real*8, dimension(3) :: step_Euler_Hookean

        B = construct_B_RPY(Q, a)

        BdotB = matmul(B, B)

        step_Euler_Hookean = Q + (ten_vec_dot(k, Q) - 0.5*ten_vec_dot(BdotB,Q))*dt &
                               + ten_vec_dot(B, dW)

    end function

    function step_semimp_FF(Q, k, dt, Q0, alpha, a, dW)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, alpha, a, dW(3)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, L, RdotB2_L, Qlength, temp_R1, temp_R2, Q_test
        real*8, dimension(3) :: F, u, RHS, Qpred, RdotB2
        real*8, dimension(3,3) :: B1, B2, B1B1, B2B2
        real*8, dimension(3) :: step_semimp_FF
        integer*8, save :: failure_count = 0

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                  + ten_vec_dot(B1,dW)

        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - ten_vec_dot(B1B1,F) + F)*dt + ten_vec_dot(B1,dW)

        L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/L

        if (L.gt.(10000.D0*sqrt(alpha)+Q0)) then
            Qlength = find_roots_cubic_newton( -(2.D0*Q0+L), -alpha+Q0**2-dt*alpha/2.D0+2.d0*L*Q0, &
                                dt*alpha*Q0/2.D0+L*(alpha-Q0**2), &
                                Q0-sqrt(alpha), Q0+sqrt(alpha), Q0+sqrt(alpha))
        else
            Qlength = find_roots( -(2.D0*Q0+L), -alpha+Q0**2-dt*alpha/2.D0+2.d0*L*Q0, &
                            dt*alpha*Q0/2.D0+L*(alpha-Q0**2), Q0-sqrt(alpha), Q0+sqrt(alpha) )
        end if

        step_semimp_FF = RHS*Qlength/L

        Q_test = 1.D0+dt/2.D0*(((Qlength-Q0)/(1.D0+(Qlength-Q0)**2/alpha))/Qlength)
        if (Q_test.le.0.D0) then
            print *, "Q prefactor is negative or 0, =", Q_test
        end if

    end function step_semimp_FF

    function step_semimp_FF_old(Q, k, dt, Q0, alpha, a, dW, dummy)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, alpha, a, dW(3), dummy
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, L, RdotB2_L, Qlength, temp_R1, temp_R2, Q_test
        real*8, dimension(3) :: F, u, RHS, Qpred, RdotB2
        real*8, dimension(3,3) :: B1, B2, B1B1, B2B2
        real*8, dimension(3) :: step_semimp_FF_old

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                  + ten_vec_dot(B1,dW)

        B2 = construct_B_RPY(Qpred, a)
        B2B2 = matmul(B2,B2)
        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                + ten_vec_dot(B1,dW)

        L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/L
        RdotB2 = ten_vec_dot(alpha*B2B2*dt/4.D0,u)
        RdotB2_L = sqrt(RdotB2(1)**2 + RdotB2(2)**2 + RdotB2(3)**2)

        if (L.gt.(10000.D0*sqrt(alpha)+Q0)) then
            Qlength = find_roots_cubic_newton( -(2.D0*Q0+L), -alpha+Q0**2-RdotB2_L+2.d0*L*Q0, &
                                RdotB2_L*Q0+L*(alpha-Q0**2), &
                                Q0-sqrt(alpha), Q0+sqrt(alpha), Q0+sqrt(alpha))
        else
            Qlength = find_roots( -(2.D0*Q0+L), -alpha+Q0**2-RdotB2_L+2.d0*L*Q0, &
                                RdotB2_L*Q0+L*(alpha-Q0**2), Q0-sqrt(alpha), Q0+sqrt(alpha) )
        end if

        !Testing prefactor
        Q_test = 1+RdotB2_L/alpha*(((Qlength-Q0)/(1.D0+(Qlength-Q0)**2/alpha))/Qlength)
        if (Q_test.le.0.D0) then
            print *, "Q prefactor is negative", Q_test
        end if

        step_semimp_FF_old = RHS*Qlength/L
    end function step_semimp_FF_old

    function find_roots(a, b, c, lower_bound, upper_bound)
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

    function find_roots_cubic_newton(a, b, c, lower_bound, upper_bound, initial_guess)
        real*8, intent(in) :: a, b, c, lower_bound, upper_bound, initial_guess
        real*8 :: find_roots_cubic_newton
        real*8 :: x0, x1
        real*8, parameter :: exit_err = 1.D-10
        integer :: i

        x0 = initial_guess
        do i=1,1000
            x1 = x0-(x0**3+a*x0**2+b*x0+c)/(3.D0*x0**2+2.D0*a*x0+b)
            if (((abs(x1-x0)).lt.exit_err).and.((x1.le.upper_bound).and.(x1.ge.lower_bound))) then
                find_roots_cubic_newton = x1
                return
            else
                x0 = x1
            end if
        end do
        print *, "didn't find root, final guess is ", x1
        find_roots_cubic_newton = upper_bound
    end function

end module
