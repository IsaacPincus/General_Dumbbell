module Integrate
    use Dumbbell_util
    implicit none


    interface step
        module procedure step_semimp_FF
        module procedure step_Euler_Hookean
        module procedure step_semimp_FF_lookup
    end interface

    contains

    pure function construct_B_ROB(Q, a)
        implicit none
        real*8, intent(in) :: Q(3), a
        real*8, dimension(3,3) :: construct_B_ROB
        real*8 :: Ql, aux, g, g_til, Ql2, Ql4, Ql6, a2, a4
        real*8, parameter :: C43 = 4.D0/3.D0
        real*8, parameter :: C83 = 8.D0/3.D0
        real*8, parameter :: C143 = 14.D0/3.D0

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

    pure function construct_B_RPY(Q, a)
        implicit none
        real*8, intent(in) :: Q(3), a
        real*8, dimension(3,3) :: construct_B_RPY
        real*8 :: Ql, aux, g, g_til, Ql2, a2, Atmp, Btmp, temp
        real*8, parameter :: C43 = 4.D0/3.D0
        real*8, parameter :: C38 = 3.D0/8.D0
        real*8, parameter :: C18 = 1.D0/8.D0
        real*8, parameter :: C23 = 2.D0/3.D0

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
        real*8, dimension(3,3) :: B, BdotB
        real*8, dimension(3) :: step_Euler_Hookean

        B = construct_B_RPY(Q, a)

        BdotB = matmul(B, B)

        step_Euler_Hookean = Q + (ten_vec_dot(k, Q) - 0.5*ten_vec_dot(BdotB,Q))*dt &
                               + ten_vec_dot(B, dW)

    end function

    function step_semimp_FF_lookup(Q, k, dt, sigma, alpha, a, dW, Yvals, Qvals)
        implicit none
        real*8, intent(in) :: Q(3), dt, sigma, alpha, a, dW(3), Yvals(:), Qvals(:)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, Y, Qlength
        real*8, dimension(3) :: F, RHS, Qpred
        real*8, dimension(3,3) :: B1, B1B1
        real*8, dimension(3) :: step_semimp_FF_lookup
        integer :: n

        n = size(Yvals)

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - sigma)/(1.0D0-(Ql-sigma)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                  + ten_vec_dot(B1,dW)

        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - ten_vec_dot(B1B1,F) + 0.5D0*F)*dt + ten_vec_dot(B1,dW)
        Y =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)

        if (Y.lt.Yvals(1)) then
            Qlength = Qvals(1)
        elseif (Y.gt.Yvals(n)) then
            Qlength = Qvals(n)
        else
            Qlength = poly_interp_bs(Yvals, Qvals, Y, 4)
        end if

        step_semimp_FF_lookup = RHS*Qlength/Y

    end function

    pure function step_semimp_FF(Q, k, dt, sigma, alpha, a, dW)
        implicit none
        real*8, intent(in) :: Q(3), dt, sigma, alpha, a, dW(3)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, Y, Qlength
        real*8, dimension(3) :: F, u, RHS, Qpred
        real*8, dimension(3,3) :: B1, B1B1
        real*8, dimension(3) :: step_semimp_FF

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - sigma)/(1.0D0-(Ql-sigma)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                  + ten_vec_dot(B1,dW)

        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - ten_vec_dot(B1B1,F) + 0.5D0*F)*dt &
                + ten_vec_dot(B1,dW)

        Y =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/Y

        Qlength = find_roots( -(2.D0*sigma+Y), -alpha+sigma**2-alpha*dt/4.D0+2.d0*Y*sigma, &
                                alpha*dt/4.D0*sigma+Y*(alpha-sigma**2), sigma-sqrt(alpha), sigma+sqrt(alpha) )

        step_semimp_FF = RHS*Qlength/Y

    end function step_semimp_FF

end module
