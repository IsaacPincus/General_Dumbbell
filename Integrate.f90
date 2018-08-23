module Integrate
    use Dumbbell_util
    implicit none


    interface step
        module procedure step_semimp_FF
        module procedure step_semimp_FF_old
        module procedure step_Euler_Hookean
        module procedure step_semimp_FF_Y_test
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

    function step_semimp_FF_lookup(Q, k, dt, Q0, alpha, a, dW, Yvals, Qvals)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, alpha, a, dW(3), Yvals(:), Qvals(:)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, Y, Qlength
        real*8, dimension(3) :: F, RHS
        real*8, dimension(3,3) :: B1, B1B1
        real*8, dimension(3) :: step_semimp_FF_lookup
        integer :: n

        n = size(Yvals)

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)

        RHS = Q + 0.5D0*(ten_vec_dot(k,2.D0*Q) - ten_vec_dot(B1B1,F) + 0.5D0*F)*dt + ten_vec_dot(B1,dW)
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

    pure function step_semimp_FF(Q, k, dt, Q0, alpha, a, dW)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, alpha, a, dW(3)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, L, Qlength
        real*8, dimension(3) :: F, u, RHS, Qpred
        real*8, dimension(3,3) :: B1, B1B1
        real*8, dimension(3) :: step_semimp_FF
!        integer*8, save :: failure_count = 0
!        integer*8, save :: iter_count = 0
!        integer*8, save :: no_calls = 0

!        no_calls = no_calls + 1

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
!        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
!                  + ten_vec_dot(B1,dW)

        Qpred = Q
!        Ql = sqrt(Qpred(1)**2 + Qpred(2)**2 + Qpred(3)**2)
!
!        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Qpred/Ql

        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - ten_vec_dot(B1B1,F) + 0.5D0*F)*dt + ten_vec_dot(B1,dW)

        L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/L

        if (L.gt.(10000.D0*sqrt(alpha)+Q0)) then
            Qlength = find_roots_cubic_newton( -(2.D0*Q0+L), -alpha+Q0**2-dt*alpha/4.D0+2.d0*L*Q0, &
                                dt*alpha*Q0/4.D0+L*(alpha-Q0**2), &
                                Q0-sqrt(alpha), Q0+sqrt(alpha), Q0+sqrt(alpha))
        else
            Qlength = find_roots( -(2.D0*Q0+L), -alpha+Q0**2-dt*alpha/4.D0+2.d0*L*Q0, &
                            dt*alpha*Q0/4.D0+L*(alpha-Q0**2), Q0-sqrt(alpha), Q0+sqrt(alpha) )
        end if

!        do while(.true.)
!            iter_count = iter_count + 1
!            Qpred = RHS*Qlength/L
!            Ql = sqrt(Qpred(1)**2 + Qpred(2)**2 + Qpred(3)**2)
!            F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Qpred/Ql
!
!            RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - ten_vec_dot(B1B1,F) + F)*dt + ten_vec_dot(B1,dW)
!
!            L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
!            u = RHS/L
!
!            if (L.gt.(10000.D0*sqrt(alpha)+Q0)) then
!                Qlength = find_roots_cubic_newton( -(2.D0*Q0+L), -alpha+Q0**2-dt*alpha/4.D0+2.d0*L*Q0, &
!                                    dt*alpha*Q0/4.D0+L*(alpha-Q0**2), &
!                                    Q0-sqrt(alpha), Q0+sqrt(alpha), Q0+sqrt(alpha))
!            else
!                Qlength = find_roots( -(2.D0*Q0+L), -alpha+Q0**2-dt*alpha/4.D0+2.d0*L*Q0, &
!                                dt*alpha*Q0/4.D0+L*(alpha-Q0**2), Q0-sqrt(alpha), Q0+sqrt(alpha) )
!            end if
!
!            if (abs(Ql-Qlength).lt.1.D-10) then
!                if (mod(no_calls, 20000000).eq.0) then
!                    print *, "avg. iterations is: ", dble(iter_count)/dble(no_calls)
!                end if
!                EXIT
!            end if
!
!        end do

        step_semimp_FF = RHS*Qlength/L

!        Q_test = 1.D0+dt/2.D0*(((Qlength-Q0)/(1.D0+(Qlength-Q0)**2/alpha))/Qlength)
!        if (Q_test.le.0.D0) then
!            print *, "Q prefactor is negative or 0, =", Q_test
!        end if

    end function step_semimp_FF

    function step_semimp_FF_Y_test(Q, k, dt, Q0, alpha, a, dW, Ymin, Ymax)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, alpha, a, dW(3), Ymin, Ymax
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, RdotB2_L, Qlength, Q_test, Y
        real*8, dimension(3) :: F, u, RHS, Qpred, RdotB2
        real*8, dimension(3,3) :: B1, B1B1
        real*8, dimension(3) :: step_semimp_FF_Y_test


        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
!        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
!                  + ten_vec_dot(B1,dW)
        Qpred = Q

        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                + ten_vec_dot(B1,dW)

        Y =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/Y
        RdotB2 = ten_vec_dot(alpha*B1B1*dt/4.D0,u)
        RdotB2_L = sqrt(RdotB2(1)**2 + RdotB2(2)**2 + RdotB2(3)**2)

        !print *, "beta =", RdotB2_L

        Qlength = find_roots_Yminmax(Y, alpha, RdotB2_L, Q0, Ymin, Ymax)

        step_semimp_FF_Y_test = RHS*Qlength/Y

        Q_test = 1.D0+dt/2.D0*(((Qlength-Q0)/(1.D0+(Qlength-Q0)**2/alpha))/Qlength)
        if (Q_test.le.0.D0) then
            print *, "Q prefactor is negative or 0, =", Q_test
            print *, Y
        end if

    end function step_semimp_FF_Y_test

    function find_roots_Yminmax(Y, alpha, beta, Q0, Ymin, Ymax)
        implicit none
        real*8, intent(in) :: Y, alpha, beta, Q0, Ymin, Ymax
        real*8 :: find_roots_Yminmax
        integer*8, save :: failure_count = 0

        if (Y.gt.Ymax) then
            find_roots_Yminmax = find_roots_cubic_newton( -(2.D0*Q0+Y), &
                                -alpha+Q0**2-beta+2.d0*Y*Q0, &
                                beta*Q0+Y*(alpha-Q0**2), &
                                Q0-sqrt(alpha), Q0+sqrt(alpha), Q0+sqrt(alpha))
            failure_count = failure_count + 1
            if (mod(failure_count, 10).eq.0) then
                print *, "failure count is ", failure_count
            end if
        elseif (Y.lt.Ymin) then
            find_roots_Yminmax = find_roots_cubic_newton( -(2.D0*Q0+Y), &
                                -alpha+Q0**2-beta+2.d0*Y*Q0, &
                                beta*Q0+Y*(alpha-Q0**2), &
                                Q0-sqrt(alpha), Q0+sqrt(alpha), Q0-sqrt(alpha))
            failure_count = failure_count + 1
            if (mod(failure_count, 1000).eq.0) then
                print *, "failure count is ", failure_count
            end if
        else
            find_roots_Yminmax = find_roots( -(2.D0*Q0+Y), -alpha+Q0**2-beta+2.d0*Y*Q0, &
                            beta*Q0+Y*(alpha-Q0**2), Q0-sqrt(alpha), Q0+sqrt(alpha) )
        end if

    end function

    pure function step_semimp_FF_old(Q, k, dt, Q0, alpha, a, dW, dummy)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, alpha, a, dW(3), dummy
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, L, RdotB2_L, Qlength
        real*8, dimension(3) :: F, u, RHS, Qpred, RdotB2
        real*8, dimension(3,3) :: B1, B1B1
        real*8, dimension(3) :: step_semimp_FF_old

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q/Ql

        B1 = construct_B_RPY(Q, a)

        B1B1 = matmul(B1,B1)
!        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
!                  + ten_vec_dot(B1,dW)
        Qpred = Q
        !Ql = sqrt(Qpred(1)**2 + Qpred(2)**2 + Qpred(3)**2)

        !F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Qpred/Ql

!        B2 = construct_B_RPY(Qpred, a)
!        B2B2 = matmul(B2,B2)
        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                + ten_vec_dot(B1,dW)

        L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/L
        RdotB2 = ten_vec_dot(alpha*B1B1*dt/4.D0,u)
        RdotB2_L = sqrt(RdotB2(1)**2 + RdotB2(2)**2 + RdotB2(3)**2)
        !print *, "beta = ", RdotB2_L
!        print *, u
!        print *, B2B2
!        print *, B1
!        print *, B2

        if (L.gt.(10000.D0*sqrt(alpha)+Q0)) then
            Qlength = find_roots_cubic_newton( -(2.D0*Q0+L), -alpha+Q0**2-RdotB2_L+2.d0*L*Q0, &
                                RdotB2_L*Q0+L*(alpha-Q0**2), &
                                Q0-sqrt(alpha), Q0+sqrt(alpha), Q0+sqrt(alpha))
        else
            Qlength = find_roots( -(2.D0*Q0+L), -alpha+Q0**2-RdotB2_L+2.d0*L*Q0, &
                                RdotB2_L*Q0+L*(alpha-Q0**2), Q0-sqrt(alpha), Q0+sqrt(alpha) )
        end if

        !Testing prefactor
!        Q_test = 1+RdotB2_L/alpha*(((Qlength-Q0)/(1.D0+(Qlength-Q0)**2/alpha))/Qlength)
!        if (Q_test.le.0.D0) then
!            print *, "Q prefactor is negative", Q_test
!        end if

        step_semimp_FF_old = RHS*Qlength/L
    end function step_semimp_FF_old

end module
