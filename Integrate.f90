module Integrate
    use Random_Numbers
    implicit none
    real*8, parameter :: PI = 4*atan(1.0D0)
    real*8, parameter, dimension(3,3) :: delT = reshape((/1, 0, 0, &
                                                          0, 1, 0, &
                                                          0, 0, 1/), &
                                                        shape(delT))

    interface step
        module procedure step_semimp_FF
        module procedure step_Euler_Hookean
    end interface

    contains

    function construct_B(Q, a)
        implicit none
        real*8, intent(in) :: Q(3), a
        real*8, dimension(3,3) :: construct_B
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
        construct_B = sqrt(g)*delT + (sqrt(g+g_til)-sqrt(g))*dyadic_prod(Q,Q)/Ql2
        return

    end function construct_B

    function step_Euler_Hookean(Q, k, dt, a, dW)
        implicit none
        real*8, intent(in) :: Q(3), dt, a, dW(3)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql
        real*8, dimension(3) :: F
        real*8, dimension(3,3) :: B, BdotB
        real*8, dimension(3) :: step_Euler_Hookean

        B = construct_B(Q, a)

        BdotB = matmul(B, B)

        step_Euler_Hookean = Q + (ten_vec_dot(k, Q) - 0.5*ten_vec_dot(BdotB,Q))*dt &
                               + ten_vec_dot(B, dW)

    end function

    function step_semimp_FF(Q, k, dt, Q0, b, a, dW)
        implicit none
        real*8, intent(in) :: Q(3), dt, Q0, b, a, dW(3)
        real*8, intent(in) :: k(3,3)
        real*8 :: Ql, L, RdotB2_L, Qlength, temp_R1, temp_R2
        real*8, dimension(3) :: F, u, RHS, Qpred, RdotB2
        real*8, dimension(3,3) :: B1, B2, B1B1, B2B2
        real*8, dimension(3) :: step_semimp_FF

        Ql = sqrt(Q(1)**2 + Q(2)**2 + Q(3)**2)
        F = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q/Ql

        B1 = construct_B(Q, a)

        B1B1 = matmul(B1,B1)
        Qpred = Q(:) + (ten_vec_dot(k,Q) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                  + ten_vec_dot(B1,dW)

        B2 = construct_B(Qpred, a)
        B2B2 = matmul(B2,B2)
        RHS = Q + 0.5D0*(ten_vec_dot(k,Q+Qpred) - 0.5D0*ten_vec_dot(B1B1,F))*dt &
                + ten_vec_dot(B1,dW)

        L =  sqrt(RHS(1)**2 + RHS(2)**2 + RHS(3)**2)
        u = RHS/L
        RdotB2 = ten_vec_dot(b*B2B2*dt/4.D0,u)
        RdotB2_L = sqrt(RdotB2(1)**2 + RdotB2(2)**2 + RdotB2(3)**2)

        Qlength = find_roots( -(2.D0*Q0+L), -b+Q0**2-RdotB2_L+2.d0*L*Q0, &
                                RdotB2_L*Q0+L*(b-Q0**2), Q0-sqrt(b), Q0+sqrt(b) )

        step_semimp_FF = RHS*Qlength/L

    end function step_semimp_FF

    pure function find_roots(a, b, c, lower_bound, upper_bound)
        implicit none
        real*8, intent(in) :: a, b, c, lower_bound, upper_bound
        real*8 :: find_roots
        real*8 :: Q, R, theta, x
        integer :: i

        Q = (a**2 - 3.D0*b)/9.D0
        R = (2.D0*a**3 - 9.D0*a*b + 27.D0*c)/54.D0

        theta = acos(R/sqrt(Q**3))

        do i=-1,1
            x = -2.D0*sqrt(Q)*cos((theta + real(i)*PI*2.D0)/3.D0)-a/3.D0
            if ((x.ge.lower_bound).and.(x.le.upper_bound)) then
                find_roots = x
                EXIT
            end if
       end do

    end function find_roots

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
