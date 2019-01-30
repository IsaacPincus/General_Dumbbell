module Dumbbell_util
    implicit none
    integer, parameter :: quad = SELECTED_REAL_KIND(32)
    real*8, parameter :: PI = 4.D0*atan(1.0D0)
    real*8, parameter, dimension(3,3) :: delT = reshape((/1, 0, 0, &
                                                          0, 1, 0, &
                                                          0, 0, 1/), &
                                                        shape(delT))

    type measured_variables
        real*8 :: AvgQ, ErrQ, S, ErrS, AvgEta, ErrEta
        real*8 :: AvgPsi, ErrPsi, AvgPsi2, ErrPsi2, AvgF, ErrF
        real*8 :: AvgChiTau, ErrChiTau, AvgChiG, ErrChiG
    end type

    contains

    subroutine measure(N, Q, avgA, avgB, varA, varB, covAB)
        implicit none
        !This details a single-pass algorithm for calculation of statistical moments,
        !which is used for a variety of material properties in the other functions
        !from https://prod-ng.sandia.gov/techlib-noauth/access-control.cgi/2008/086212.pdf
        integer*8, intent(in) :: N
        integer :: i
        real*8, intent(in), dimension(:,:) :: Q
        real*8, intent(out) :: avgA, avgB, varA, varB, covAB
        real*8 :: avgAtmp, avgBtmp, covABtmp, varAtmp, varBtmp, A, B, delA, delB, li
!        real*8 :: muA1, muB1, muA2, muB2
        real*8, save :: si


        !shared between threads in external openMP construct
        avgA = 0.D0; avgB = 0.D0; covAB = 0.D0; varA = 0.D0; varB = 0.D0

        !private to subroutine and hence private to threads
        avgAtmp = 0.D0; avgBtmp = 0.D0; covABtmp = 0.D0; varAtmp = 0.D0; varBtmp = 0.D0

        !li (local increment) has a value for each thread, keeps track of number of samples
        !per thread
        li = 0.D0
        !si (shared increment) is declared using 'save' and so is shared between threads, used to
        !keep track of the current number of samples in moment during reduction
        si = 0.D0

        !$OMP Do
        do i=1,N
            li = li + 1

            A = Q(1,i)*Q(2,i)
            B = Q(1,i)*Q(1,i) - Q(2,i)*Q(2,i)
            delA = A - avgAtmp
            delB = B - avgBtmp
            avgAtmp = avgAtmp + delA/li
            avgBtmp = avgBtmp + delB/li

            varAtmp = varAtmp + delA*(A-avgAtmp)
            varBtmp = varBtmp + delB*(B-avgBtmp)

!            covABtmp = (covABtmp*(li-1) + (1.D0/li)*(A*(li-1) - avgAtmp)*(B*(li-1) - avgBtmp))/li
            covABtmp = covABtmp + (li-1)/(li)*delA*delB

        end do
        !$OMP END DO

        !$OMP CRITICAL
        delA = avgAtmp-avgA
        delB = avgBtmp-avgB
        avgA = avgA + li*(delA)/(li+si)
        avgB = avgB + li*(delB)/(li+si)
        varA = varA + varAtmp + (li*si)*delA**2/(li+si)
        varB = varB + varBtmp + (li*si)*delB**2/(li+si)

        covAB = covAB + covABtmp + (delA)*(delB)*(li*si)/(li+si)

        si = si + li
        !$OMP END CRITICAL

        !$OMP barrier

        !$OMP SINGLE
        varA = varA/(N-1)
        varB = varB/(N-1)
        covAB = covAB/(N-1)
        !$OMP END SINGLE



    end subroutine

    subroutine measure_shear_no_VR(meas, Q, Q0, alpha, sr, Ntraj)
        implicit none
        type(measured_variables), intent(out) :: meas
        integer*8, intent(in) :: Ntraj
        integer :: i
        real*8, intent(in) :: Q0, alpha, Q(3,Ntraj), sr
        real*8 :: tau(3,3), F(3), Bs, Ql, Ql2, Fl2, Fl, S, inTan, li
        !temporary variables, local to threads
        real*8 :: AvgEtaTmp, VarEtaTmp, AvgPsiTmp, VarPsiTmp, AvgPsi2Tmp, VarPsi2Tmp
        real*8 :: AvgQl2Tmp, VarQl2Tmp, AvgFTmp, VarFTmp, AvgSTmp, VarSTmp
        real*8 :: AvgQxQyTmp, VarQxQyTmp
        real*8 :: AvgTauxyTmp, VarTauxyTmp
        real*8 :: CovQxyAndQxxminyyTmp, CovTauxyAndTauxxminyyTmp
        real*8 :: delEta, delPsi, delPsi2, delS, delF, delQl2, delQxQy, delTauxy
        real*8 :: AvgTauxxminyyTmp, VarTauxxminyyTmp, delTauxxminyy
        real*8 :: AvgTauyyminzzTmp, VarTauyyminzzTmp, delTauyyminzz
        real*8 :: AvgQxxminyyTmp, VarQxxminyyTmp, delQxxminyy
        !shared variables created using save
        real*8, save :: AvgEta, VarEta, AvgPsi, VarPsi, AvgPsi2, VarPsi2
        real*8, save :: AvgQl2, VarQl2, AvgF, VarF, AvgS, VarS
        real*8, save:: AvgQxQy, VarQxQy
        real*8, save :: AvgTauxy, VarTauxy
        real*8, save :: AvgTauxxminyy, VarTauxxminyy
        real*8, save :: AvgTauyyminzz, VarTauyyminzz
        real*8, save :: AvgQxxminyy, VarQxxminyy
        real*8, save :: CovQxyAndQxxminyy, CovTauxyAndTauxxminyy
        real*8, save :: si

        ! These variables are all global and shared between threads
        AvgEta = 0.D0; VarEta = 0.D0; AvgPsi = 0.D0; VarPsi = 0.D0;
        AvgPsi2 = 0.D0; VarPsi2 = 0.D0
        AvgQl2 = 0.D0; VarQl2 = 0.D0; AvgF = 0.D0; VarF = 0.D0;
        AvgS = 0.D0; VarS = 0.D0
        AvgQxQy = 0.D0; VarQxQy = 0.D0
        AvgTauxy = 0.D0; VarTauxy = 0.D0
        AvgTauxxminyy = 0.D0; VarTauxxminyy = 0.D0;
        AvgTauyyminzz = 0.D0; VarTauyyminzz = 0.D0;
        AvgQxxminyy = 0.D0; VarQxxminyy = 0.D0;
        si = 0.D0

        !These variables SHOULD(!) be private
        tau = 0.D0
        li = 0.D0 ! stands for local i, aka private increment variable for each thread
        AvgEtaTmp = 0.D0; VarEtaTmp = 0.D0; AvgPsiTmp = 0.D0; VarPsiTmp = 0.D0;
        AvgPsi2Tmp = 0.D0; VarPsi2Tmp = 0.D0
        AvgQl2Tmp = 0.D0; VarQl2Tmp = 0.D0; AvgFTmp = 0.D0; VarFTmp = 0.D0;
        AvgSTmp = 0.D0; VarSTmp = 0.D0
        AvgQxQyTmp = 0.D0; VarQxQyTmp = 0.D0
        AvgTauxyTmp = 0.D0; VarTauxyTmp = 0.D0
        AvgTauxxminyyTmp = 0.D0; VarTauxxminyyTmp = 0.D0;
        AvgTauyyminzzTmp = 0.D0; VarTauyyminzzTmp = 0.D0;
        AvgQxxminyyTmp = 0.D0; VarQxxminyyTmp = 0.D0;


        !$OMP DO
        do i=1,Ntraj
            li = li + 1.D0

            !Calculate required quantities
            Ql2 = Q(1,i)**2.D0 + Q(2,i)**2.D0 + Q(3,i)**2.D0
            Ql = sqrt(Ql2)
            S = 0.5D0*(3.D0*dot_product( (/1,0,0/), Q(:,i))**2/Ql2-1.D0)
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q(:,i)/Ql
            Fl2 = F(1)**2 + F(2)**2 + F(3)**2
            Fl = sqrt(Fl2)
            tau(:,:) = dyadic_prod(Q(:,i), F)

            !Calculate on-line average, variance and covariance
            delTauxy = tau(1,2) - AvgTauxyTmp
            delTauxxminyy = tau(1,1) - tau(2,2) - AvgTauxxminyyTmp
            delTauyyminzz = tau(2,2) - tau(3,3) - AvgTauyyminzzTmp
            delS = S - AvgSTmp
            delQl2 = Ql2 - AvgQl2Tmp
            delF = Fl - AvgFTmp
            delQxQy = Q(1,i)*Q(2,i) - AvgQxQyTmp
            delQxxminyy = Q(1,i)*Q(1,i) - Q(2,i)*Q(2,i) - AvgQxxminyyTmp

            AvgTauxyTmp = AvgTauxyTmp           + delTauxy/li
            AvgTauxxminyyTmp = AvgTauxxminyyTmp + delTauxxminyy/li
            AvgTauyyminzzTmp = AvgTauyyminzzTmp + delTauyyminzz/li
            AvgSTmp = AvgSTmp                   + delS/li
            AvgQl2Tmp = AvgQl2Tmp               + delQl2/li
            AvgFTmp = AvgFTmp                   + delF/li
            AvgQxQyTmp = AvgQxQyTmp             + DelQxQy/li
            AvgQxxminyyTmp = AvgQxxminyyTmp     + delQxxminyy/li

            VarTauxyTmp = VarTauxyTmp           + delTauxy*(tau(1,2) - AvgTauxyTmp)
            VarTauxxminyyTmp = VarTauxxminyyTmp + delTauxxminyy*&
                                                &(tau(1,1) - tau(2,2) - AvgTauxxminyyTmp)
            VarTauyyminzzTmp = VarTauyyminzzTmp + delTauyyminzz*&
                                                &(tau(2,2) - tau(3,3) - AvgTauyyminzzTmp)
            VarSTmp = VarSTmp                   + delS*(S - AvgSTmp)
            VarQl2Tmp = VarQl2Tmp               + delQl2*(Ql2 - AvgQl2Tmp)
            VarFTmp = VarFTmp                   + delF*(Fl - AvgFTmp)
            VarQxQyTmp = VarQxQyTmp             + DelQxQy*(Q(1,i)*Q(2,i) - AvgQxQyTmp)
            VarQxxminyyTmp = VarQxxminyyTmp     + delQxxminyy*&
                                                &(Q(1,i)*Q(1,i) - Q(2,i)*Q(2,i) - AvgQxxminyyTmp)

            CovTauxyAndTauxxminyyTmp = CovTauxyAndTauxxminyyTmp + &
                        &(li-1.D0)/li*delTauxy*delTauxxminyy
            CovQxyAndQxxminyyTmp = CovQxyAndQxxminyyTmp + &
                        &(li-1.D0)/li*delQxQy*delQxxminyy

        end do
        !$OMP END DO

        !Combine results from threads into shared variables
        !$OMP critical
        delTauxy = AvgTauxyTmp - AvgTauxy
        delTauxxminyy = AvgTauxxminyyTmp - AvgTauxxminyy
        delTauyyminzz = AvgTauyyminzzTmp - AvgTauyyminzz
        delS = AvgSTmp - AvgS
        delQl2 = AvgQl2Tmp - AvgQl2
        delF = AvgFTmp - AvgF
        delQxQy = AvgQxQyTmp - AvgQxQy
        delQxxminyy = AvgQxxminyyTmp - AvgQxxminyy

        AvgTauxy = AvgTauxy + li*delTauxy/(li+si)
        AvgTauxxminyy = AvgTauxxminyy + li*delTauxxminyy/(li+si)
        AvgTauyyminzz = AvgTauyyminzz + li*delTauyyminzz/(li+si)
        AvgS = AvgS + li*delS/(li+si)
        AvgQl2 = AvgQl2 + li*delQl2/(li+si)
        AvgF = AvgF + li*delF/(li+si)
        AvgQxQy = AvgQxQy + li*delQxQy/(li+si)
        AvgQxxminyy = AvgQxxminyy + li*delQxxminyy/(li+si)

        VarTauxy = VarTauxy + VarTauxyTmp + (li*si)*delTauxy**2/(li+si)
        VarTauxxminyy = VarTauxxminyy + VarTauxxminyyTmp +&
                            & (li*si)*delTauxxminyy**2/(li+si)
        VarTauyyminzz = VarTauyyminzz + VarTauyyminzzTmp +&
                            & (li*si)*delTauyyminzz**2/(li+si)
        VarS = VarS + VarSTmp + (li*si)*delS**2/(li+si)
        VarQl2 = VarQl2 + VarQl2Tmp + (li*si)*delQl2**2/(li+si)
        VarF = VarF + VarFTmp + (li*si)*delF**2/(li+si)
        VarQxQy = VarQxQy + VarQxQyTmp + (li*si)*delQxQy**2/(li+si)
        VarQxxminyy = VarQxxminyy + VarQxxminyyTmp + &
                            & (li*si)*delQxxminyy**2/(li+si)

        CovQxyAndQxxminyy = CovQxyAndQxxminyy + CovQxyAndQxxminyyTmp + &
                    & delQxQy*delQxxminyy*(li*si)/(li+si)
        CovTauxyAndTauxxminyy = CovTauxyAndTauxxminyy + CovTauxyAndTauxxminyyTmp + &
                    & delTauxy*delTauxxminyy*(li*si)/(li+si)

        si = si + li
        !$OMP END critical

        !$OMP barrier

        !$OMP single
        !Quantities above are technically moments, not variances
        !We divide by N-1 to give an unbiased estimator of the variance
        VarTauxy = VarTauxy/(Ntraj-1)
        VarTauxxminyy = VarTauxxminyy/(Ntraj-1)
        VarTauyyminzz = VarTauyyminzz/(Ntraj-1)
        VarS = VarS/(Ntraj-1)
        VarQl2 = VarQl2/(Ntraj-1)
        VarF = VarF/(Ntraj-1)
        VarQxQy = VarQxQy/(Ntraj-1)
        VarQxxminyy = VarQxxminyy/(Ntraj-1)
        CovQxyAndQxxminyy = CovQxyAndQxxminyy/(Ntraj-1)
        CovTauxyAndTauxxminyy = CovTauxyAndTauxxminyy/(Ntraj-1)

        !Convert from variances to standard errors
        meas%AvgQ = sqrt(AvgQl2)
        meas%ErrQ = 0.5D0*(1.D0/meas%AvgQ)*sqrt(VarQl2/Ntraj)

        meas%S = AvgS
        meas%ErrS = sqrt(VarS/Ntraj)

        meas%AvgF = AvgF
        meas%ErrF = sqrt(VarF/Ntraj)

        meas%AvgEta = AvgTauxy/sr
        meas%ErrEta = sqrt(VarTauxy/(Ntraj*sr**2))

        meas%AvgPsi = AvgTauxxminyy/sr**2
        meas%ErrPsi = sqrt(VarTauxxminyy/(Ntraj*sr**4))

        meas%AvgPsi2 = AvgTauyyminzz/sr**2
        meas%ErrPsi2 = sqrt(VarTauyyminzz/(Ntraj*sr**4))

        inTan = 2.D0*AvgTauxy/(AvgTauxxminyy)
        meas%AvgChiTau = 0.5D0*atan(inTan)
        meas%ErrChiTau = 0.5D0*abs(inTan)/(1.D0+inTan**2)*&
                        &sqrt(VarTauxxminyy/AvgTauxxminyy**2 +&
                            & VarTauxy/AvgTauxy**2 - &
                            & 2.D0*CovTauxyAndTauxxminyy/(AvgTauxy*AvgTauxxminyy))/&
                            & sqrt(dble(Ntraj))

        inTan = 2.D0*AvgQxQy/AvgQxxminyy
        meas%AvgChiG = 0.5D0*atan(inTan)
        meas%ErrChiG = 0.5D0*abs(inTan)/(1.D0+inTan**2)*&
                        sqrt(VarQxxminyy/AvgQxxminyy**2 +&
                            & VarTauxy/AvgQxQy**2 - &
                            & 2.D0*CovQxyAndQxxminyy/(AvgQxQy*AvgQxxminyy))/&
                            & sqrt(dble(Ntraj))
        !$OMP end single

    end subroutine

    subroutine measure_shear_with_VR(meas, Q, Q_eq_VR, Q0, alpha, sr, Ntraj)
        implicit none
        type(measured_variables), intent(out) :: meas
        integer*8, intent(in) :: Ntraj
        integer :: i
        real*8, intent(in) :: Q0, alpha, Q(3,Ntraj), Q_eq_VR(3, Ntraj), sr
        real*8 :: tau(3,3), F(3), Bs, Ql, Ql2, Fl2, Fl, S, inTan, li
        real*8 :: Btauxy, Btauxxminyy, Btauyyminzz, Bqxy, Bqxxminyy
        !temporary variables, local to threads
        real*8 :: AvgEtaTmp, VarEtaTmp, AvgPsiTmp, VarPsiTmp, AvgPsi2Tmp, VarPsi2Tmp
        real*8 :: AvgQl2Tmp, VarQl2Tmp, AvgFTmp, VarFTmp, AvgSTmp, VarSTmp
        real*8 :: AvgQxQyTmp, VarQxQyTmp
        real*8 :: AvgTauxyTmp, VarTauxyTmp
        real*8 :: CovQxyAndQxxminyyTmp, CovTauxyAndTauxxminyyTmp
        real*8 :: delEta, delPsi, delPsi2, delS, delF, delQl2, delQxQy, delTauxy
        real*8 :: AvgTauxxminyyTmp, VarTauxxminyyTmp, delTauxxminyy
        real*8 :: AvgTauyyminzzTmp, VarTauyyminzzTmp, delTauyyminzz
        real*8 :: AvgQxxminyyTmp, VarQxxminyyTmp, delQxxminyy
        !shared variables created using save
        real*8, save :: AvgEta, VarEta, AvgPsi, VarPsi, AvgPsi2, VarPsi2
        real*8, save :: AvgQl2, VarQl2, AvgF, VarF, AvgS, VarS
        real*8, save:: AvgQxQy, VarQxQy
        real*8, save :: AvgTauxy, VarTauxy
        real*8, save :: AvgTauxxminyy, VarTauxxminyy
        real*8, save :: AvgTauyyminzz, VarTauyyminzz
        real*8, save :: AvgQxxminyy, VarQxxminyy
        real*8, save :: CovQxyAndQxxminyy, CovTauxyAndTauxxminyy
        real*8, save :: si

        ! These variables are all global and shared between threads
        AvgEta = 0.D0; VarEta = 0.D0; AvgPsi = 0.D0; VarPsi = 0.D0;
        AvgPsi2 = 0.D0; VarPsi2 = 0.D0
        AvgQl2 = 0.D0; VarQl2 = 0.D0; AvgF = 0.D0; VarF = 0.D0;
        AvgS = 0.D0; VarS = 0.D0
        AvgQxQy = 0.D0; VarQxQy = 0.D0
        AvgTauxy = 0.D0; VarTauxy = 0.D0
        AvgTauxxminyy = 0.D0; VarTauxxminyy = 0.D0;
        AvgTauyyminzz = 0.D0; VarTauyyminzz = 0.D0;
        AvgQxxminyy = 0.D0; VarQxxminyy = 0.D0;
        si = 0.D0

        !These variables SHOULD(!) be private
        tau = 0.D0
        li = 0.D0 ! stands for local i, aka private increment variable for each thread
        AvgEtaTmp = 0.D0; VarEtaTmp = 0.D0; AvgPsiTmp = 0.D0; VarPsiTmp = 0.D0;
        AvgPsi2Tmp = 0.D0; VarPsi2Tmp = 0.D0
        AvgQl2Tmp = 0.D0; VarQl2Tmp = 0.D0; AvgFTmp = 0.D0; VarFTmp = 0.D0;
        AvgSTmp = 0.D0; VarSTmp = 0.D0
        AvgQxQyTmp = 0.D0; VarQxQyTmp = 0.D0
        AvgTauxyTmp = 0.D0; VarTauxyTmp = 0.D0
        AvgTauxxminyyTmp = 0.D0; VarTauxxminyyTmp = 0.D0;
        AvgTauyyminzzTmp = 0.D0; VarTauyyminzzTmp = 0.D0;
        AvgQxxminyyTmp = 0.D0; VarQxxminyyTmp = 0.D0;


        !$OMP DO
        do i=1,Ntraj
            li = li + 1.D0

            !Calculate required quantities
            Ql2 = Q(1,i)**2.D0 + Q(2,i)**2.D0 + Q(3,i)**2.D0
            Ql = sqrt(Ql2)
            S = 0.5D0*(3.D0*dot_product( (/1,0,0/), Q(:,i))**2/Ql2-1.D0)
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q(:,i)/Ql
            Fl2 = F(1)**2 + F(2)**2 + F(3)**2
            Fl = sqrt(Fl2)
            tau(:,:) = dyadic_prod(Q(:,i), F)

            Btauxy = tau(1,2)
            Btauxxminyy = tau(1,1) - tau(2,2)
            Btauyyminzz = tau(2,2) - tau(3,3)
            Bqxy = Q(1,i)*Q(2,i)
            Bqxxminyy = Q(1,i)*Q(1,i) - Q(2,i)*Q(2,i)

            !subtract equilibrium values from shear-flow values
            Ql2 = Q_eq_VR(1,i)**2 + Q_eq_VR(2,i)**2 + Q_eq_VR(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product( (/1,0,0/), Q_eq_VR(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/alpha)*Q_eq_VR(:,i)/Ql
            tau(:,:) = dyadic_prod(Q_eq_VR(:,i), F)

            Btauxy = Btauxy - tau(1,2)
            Btauxxminyy = Btauxxminyy - (tau(1,1) - tau(2,2))
            Btauyyminzz = Btauyyminzz - (tau(2,2) - tau(3,3))
            Bqxy = Bqxy - (Q_eq_VR(1,i)*Q_eq_VR(2,i))
            Bqxxminyy = Bqxxminyy - (Q_eq_VR(1,i)*Q_eq_VR(1,i) - Q_eq_VR(2,i)*Q_eq_VR(2,i))

            !Calculate on-line average, variance and covariance
            delTauxy = Btauxy - AvgTauxyTmp
            delTauxxminyy = Btauxxminyy - AvgTauxxminyyTmp
            delTauyyminzz = Btauyyminzz - AvgTauyyminzzTmp
            delS = S - AvgSTmp
            delQl2 = Ql2 - AvgQl2Tmp
            delF = Fl - AvgFTmp
            delQxQy = Bqxy - AvgQxQyTmp
            delQxxminyy = Bqxxminyy - AvgQxxminyyTmp

            AvgTauxyTmp = AvgTauxyTmp           + delTauxy/li
            AvgTauxxminyyTmp = AvgTauxxminyyTmp + delTauxxminyy/li
            AvgTauyyminzzTmp = AvgTauyyminzzTmp + delTauyyminzz/li
            AvgSTmp = AvgSTmp                   + delS/li
            AvgQl2Tmp = AvgQl2Tmp               + delQl2/li
            AvgFTmp = AvgFTmp                   + delF/li
            AvgQxQyTmp = AvgQxQyTmp             + DelQxQy/li
            AvgQxxminyyTmp = AvgQxxminyyTmp     + delQxxminyy/li

            VarTauxyTmp = VarTauxyTmp           + delTauxy*(Btauxy - AvgTauxyTmp)
            VarTauxxminyyTmp = VarTauxxminyyTmp + delTauxxminyy*&
                                                &(Btauxxminyy - AvgTauxxminyyTmp)
            VarTauyyminzzTmp = VarTauyyminzzTmp + delTauyyminzz*&
                                                &(Btauyyminzz - AvgTauyyminzzTmp)
            VarSTmp = VarSTmp                   + delS*(S - AvgSTmp)
            VarQl2Tmp = VarQl2Tmp               + delQl2*(Ql2 - AvgQl2Tmp)
            VarFTmp = VarFTmp                   + delF*(Fl - AvgFTmp)
            VarQxQyTmp = VarQxQyTmp             + DelQxQy*(Bqxy - AvgQxQyTmp)
            VarQxxminyyTmp = VarQxxminyyTmp     + delQxxminyy*&
                                                &(Bqxxminyy - AvgQxxminyyTmp)

            CovTauxyAndTauxxminyyTmp = CovTauxyAndTauxxminyyTmp + &
                        &(li-1.D0)/li*delTauxy*delTauxxminyy
            CovQxyAndQxxminyyTmp = CovQxyAndQxxminyyTmp + &
                        &(li-1.D0)/li*delQxQy*delQxxminyy

        end do
        !$OMP END DO

        !Combine results from threads into shared variables
        !$OMP critical
        delTauxy = AvgTauxyTmp - AvgTauxy
        delTauxxminyy = AvgTauxxminyyTmp - AvgTauxxminyy
        delTauyyminzz = AvgTauyyminzzTmp - AvgTauyyminzz
        delS = AvgSTmp - AvgS
        delQl2 = AvgQl2Tmp - AvgQl2
        delF = AvgFTmp - AvgF
        delQxQy = AvgQxQyTmp - AvgQxQy
        delQxxminyy = AvgQxxminyyTmp - AvgQxxminyy

        AvgTauxy = AvgTauxy + li*delTauxy/(li+si)
        AvgTauxxminyy = AvgTauxxminyy + li*delTauxxminyy/(li+si)
        AvgTauyyminzz = AvgTauyyminzz + li*delTauyyminzz/(li+si)
        AvgS = AvgS + li*delS/(li+si)
        AvgQl2 = AvgQl2 + li*delQl2/(li+si)
        AvgF = AvgF + li*delF/(li+si)
        AvgQxQy = AvgQxQy + li*delQxQy/(li+si)
        AvgQxxminyy = AvgQxxminyy + li*delQxxminyy/(li+si)

        VarTauxy = VarTauxy + VarTauxyTmp + (li*si)*delTauxy**2/(li+si)
        VarTauxxminyy = VarTauxxminyy + VarTauxxminyyTmp +&
                            & (li*si)*delTauxxminyy**2/(li+si)
        VarTauyyminzz = VarTauyyminzz + VarTauyyminzzTmp +&
                            & (li*si)*delTauyyminzz**2/(li+si)
        VarS = VarS + VarSTmp + (li*si)*delS**2/(li+si)
        VarQl2 = VarQl2 + VarQl2Tmp + (li*si)*delQl2**2/(li+si)
        VarF = VarF + VarFTmp + (li*si)*delF**2/(li+si)
        VarQxQy = VarQxQy + VarQxQyTmp + (li*si)*delQxQy**2/(li+si)
        VarQxxminyy = VarQxxminyy + VarQxxminyyTmp + &
                            & (li*si)*delQxxminyy**2/(li+si)

        CovQxyAndQxxminyy = CovQxyAndQxxminyy + CovQxyAndQxxminyyTmp + &
                    & delQxQy*delQxxminyy*(li*si)/(li+si)
        CovTauxyAndTauxxminyy = CovTauxyAndTauxxminyy + CovTauxyAndTauxxminyyTmp + &
                    & delTauxy*delTauxxminyy*(li*si)/(li+si)

        si = si + li
        !$OMP END critical

        !$OMP barrier

        !$OMP single
        !Quantities above are technically moments, not variances
        !We divide by N-1 to give an unbiased estimator of the variance
        VarTauxy = VarTauxy/(Ntraj-1)
        VarTauxxminyy = VarTauxxminyy/(Ntraj-1)
        VarTauyyminzz = VarTauyyminzz/(Ntraj-1)
        VarS = VarS/(Ntraj-1)
        VarQl2 = VarQl2/(Ntraj-1)
        VarF = VarF/(Ntraj-1)
        VarQxQy = VarQxQy/(Ntraj-1)
        VarQxxminyy = VarQxxminyy/(Ntraj-1)
        CovQxyAndQxxminyy = CovQxyAndQxxminyy/(Ntraj-1)
        CovTauxyAndTauxxminyy = CovTauxyAndTauxxminyy/(Ntraj-1)

        !Convert from variances to standard errors
        meas%AvgQ = sqrt(AvgQl2)
        meas%ErrQ = 0.5D0*(1.D0/meas%AvgQ)*sqrt(VarQl2/Ntraj)

        meas%S = AvgS
        meas%ErrS = sqrt(VarS/Ntraj)

        meas%AvgF = AvgF
        meas%ErrF = sqrt(VarF/Ntraj)

        meas%AvgEta = AvgTauxy/sr
        meas%ErrEta = sqrt(VarTauxy/(Ntraj*sr**2))

        meas%AvgPsi = AvgTauxxminyy/sr**2
        meas%ErrPsi = sqrt(VarTauxxminyy/(Ntraj*sr**4))

        meas%AvgPsi2 = AvgTauyyminzz/sr**2
        meas%ErrPsi2 = sqrt(VarTauyyminzz/(Ntraj*sr**4))

        inTan = 2.D0*AvgTauxy/(AvgTauxxminyy)
        meas%AvgChiTau = 0.5D0*atan(inTan)
        meas%ErrChiTau = 0.5D0*abs(inTan)/(1.D0+inTan**2)*&
                        &sqrt(VarTauxxminyy/AvgTauxxminyy**2 +&
                            & VarTauxy/AvgTauxy**2 - &
                            & 2.D0*CovTauxyAndTauxxminyy/(AvgTauxy*AvgTauxxminyy))/&
                            & sqrt(dble(Ntraj))

        inTan = 2.D0*AvgQxQy/AvgQxxminyy
        meas%AvgChiG = 0.5D0*atan(inTan)
        meas%ErrChiG = 0.5D0*abs(inTan)/(1.D0+inTan**2)*&
                        sqrt(VarQxxminyy/AvgQxxminyy**2 +&
                            & VarTauxy/AvgQxQy**2 - &
                            & 2.D0*CovQxyAndQxxminyy/(AvgQxQy*AvgQxxminyy))/&
                            & sqrt(dble(Ntraj))
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

!    function rand_norm(seed)
!        implicit none
!
!
!    end function

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
