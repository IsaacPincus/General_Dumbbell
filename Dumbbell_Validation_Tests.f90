Program Dumbbell_Validation_tests
    use Dumbbell_util
    use Generate_Initial_Conditions
    use Integrate
    use fruit
    use omp_lib
    implicit none

    call init_fruit
    call Unit_tests
    call Validation_tests
    call fruit_summary
    call fruit_finalize

    contains

    subroutine Unit_tests()
        implicit none

    end subroutine

    subroutine Validation_tests()
        implicit none
        integer :: k, n
        real*8 :: k_r, n_r
        real*8 :: prob
        ! Validation tests pass if results are within two
        ! standard deviations of target result,
        ! i.e. a 95% confidence level
        k = 0
        n = 0
        print *, "We expect about 5% of tests to fail by chance, or ~1/20"
        print *, ""
        call test_eq_hookean_semimp(0.01D0, 1000, 10000)
        call test_Hookean_viscosity_semimp(0.01D0, 100, 100000)
        call test_Hookean_psi2_with_HI_semimp(0.05D0, 1000, 100000)
        call test_Hookean_psi2_with_HI_semimp_2nd_method(0.05D0, 1000, 100000)
        call test_eq_FENE_semimp(0.01D0, 1000, 10000)
        call test_semimp_euler_equal(0.01D0, 100, 10000)
        call test_FENE_HI_shear_semimp_vs_Kailash_code(0.01D0, 1000, 10000)

        !p-test on likelihood of k or more 'failures' in n trials
        call get_failed_count(k)
        call get_total_count(n)
        k_r = k
        n_r = n
        !Incomplete beta function gives cumulative binomial probability
        !https://en.wikipedia.org/wiki/Binomial_distribution#Cumulative_distribution_function
        prob = 1.D0 - betai(n_r-k_r+1.D0, k_r, 0.95D0)

        print "(A, F5.1, A, I2, A, I3, A)", &
         "There is a ", prob*100.D0, "% chance of getting", k, " or more failures in", &
          n, " tests, assuming the underlying failure probability is 5% (2 sigma)"
        print *, ""

    end subroutine

    subroutine test_eq_hookean_semimp(dt, Nsteps, Ntraj)
        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, Bq, Bs, Qavg, Vqavg, S, Serr, start_time, stop_time
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running zero-shear hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        k(:,:) = 0.D0

        delX = (/1.D0, 0.D0, 0.D0/)

        Qavg = 0.D0
        Vqavg = 0.D0
        S = 0.D0
        Serr = 0.D0

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, 0.D0, dW)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product(delX, Q(:,i))**2/Ql2

            Qavg = Qavg + Ql2
            Vqavg = Vqavg + Ql
            S = S + 0.5D0*(3.D0*Bs - 1)
            Serr = Serr + 0.25D0*(9.D0*Bs**2 - 6.D0*Bs + 1.D0)
        end do

        Qavg = sqrt(Qavg/Ntraj)
        Vqavg = Vqavg/Ntraj
        Vqavg = sqrt((Qavg**2 - Vqavg**2)/(Ntraj-1))

        S = S/Ntraj
        Serr = Serr/Ntraj
        Serr = sqrt((Serr - S**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Q-sqrt(3) = ", (Qavg-sqrt(3.D0)), " +- ", Vqavg
        print *, "S = ", S, " +- ", Serr

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(3.D0), Qavg, Vqavg*2.D0, &
                          "Qavg != sqrt(3) for sr=0 Hookean Dumbbell (semimp)")
        call assertEquals(0.D0, S, Serr*2.D0, &
                          "S != 0 for sr=0 Hookean Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_Hookean_viscosity_semimp(dt, Nsteps, Ntraj)
        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr
        real*8 :: Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2, Beta, Bpsi, Bpsi2
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running sr=1 hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        sr = 1.D0;

        k(:,:) = 0.D0
        k(1,2) = 1.D0
        k = sr*k

        delX = (/1.D0, 0.D0, 0.D0/)

        Aeta = 0.D0
        Apsi = 0.D0
        Apsi2 = 0.D0
        Veta = 0.D0
        Vpsi = 0.D0
        Vpsi2 = 0.D0

        tau = 0.D0

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, 0.D0, dW)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            F(:) = (Ql - 0.D0)/(1.0D0-(Ql-0.D0)**2/10000.D0)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            Beta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            Aeta = Aeta + Beta
            Apsi = Apsi + Bpsi
            Apsi2 = Apsi2 + Bpsi2
            Veta = Veta + Beta**2
            Vpsi = Vpsi + Bpsi**2
            Vpsi2 = Vpsi2 + Bpsi2**2
        end do

        Aeta = Aeta/(Ntraj*sr)
        Veta = Veta/(Ntraj*sr**2)
        Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

        Apsi = Apsi/(Ntraj*sr**2)
        Vpsi = Vpsi/(Ntraj*sr**4)
        Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

        Apsi2 = Apsi2/(Ntraj*sr**2)
        Vpsi2 = Vpsi2/(Ntraj*sr**4)
        Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Aeta-0.63212 = ", (Aeta-0.63212D0), " +- ", Veta
        print *, "Apsi-0.52848 = ", (Apsi-0.52848D0), " +- ", Vpsi

        !Check that both Qavg and S are within acceptable range
        call assertEquals(0.63212D0, Aeta, Veta*2.D0, &
                          "eta != 0.63212 for sr=1 Hookean Dumbbell (semimp)")
        call assertEquals(0.52848D0, Apsi, Vpsi*2.D0, &
                          "psi_1 != 0.52848 for sr=1 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_Hookean_psi2_with_HI_semimp(dt, Nsteps, Ntraj)

        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, h, a
        real*8 :: Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2, Beta, Bpsi, Bpsi2
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running sr=1, hstar=0.15 hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        sr = 1.D0;

        h = 0.15D0
        a = h*sqrt(PI)

        k(:,:) = 0.D0
        k(1,2) = 1.D0
        k = sr*k

        delX = (/1.D0, 0.D0, 0.D0/)

        Aeta = 0.D0
        Apsi = 0.D0
        Apsi2 = 0.D0
        Veta = 0.D0
        Vpsi = 0.D0
        Vpsi2 = 0.D0

        tau = 0.D0

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, a, dW)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            F(:) = (Ql - 0.D0)/(1.0D0-(Ql-0.D0)**2/10000.D0)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            Beta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            Aeta = Aeta + Beta
            Apsi = Apsi + Bpsi
            Apsi2 = Apsi2 + Bpsi2
            Veta = Veta + Beta**2
            Vpsi = Vpsi + Bpsi**2
            Vpsi2 = Vpsi2 + Bpsi2**2
        end do

        Aeta = Aeta/(Ntraj*sr)
        Veta = Veta/(Ntraj*sr**2)
        Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

        Apsi = Apsi/(Ntraj*sr**2)
        Vpsi = Vpsi/(Ntraj*sr**4)
        Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

        Apsi2 = Apsi2/(Ntraj*sr**2)
        Vpsi2 = Vpsi2/(Ntraj*sr**4)
        Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Apsi2+0.01348 = ", (Apsi2+0.01348D0), " +- ", Vpsi2

        !Check that both Qavg and S are within acceptable range
        call assertEquals(-0.01348D0, Apsi2, Vpsi2*2.D0, &
                          "psi_2 != -0.01348 for sr=1, hstar=0.15 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_Hookean_psi2_with_HI_semimp_2nd_method(dt, Nsteps, Ntraj)

        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, h, a
        real*8 :: Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2, Beta, Bpsi, Bpsi2
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running sr=1, hstar=0.15 hookean dumbbell with semimp integrator test, 2nd method"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        sr = 1.D0;

        h = 0.15D0
        a = h*sqrt(PI)

        k(:,:) = 0.D0
        k(1,2) = 1.D0
        k = sr*k

        delX = (/1.D0, 0.D0, 0.D0/)

        Aeta = 0.D0
        Apsi = 0.D0
        Apsi2 = 0.D0
        Veta = 0.D0
        Vpsi = 0.D0
        Vpsi2 = 0.D0

        tau = 0.D0

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, a, dW, 5.D0)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            F(:) = (Ql - 0.D0)/(1.0D0-(Ql-0.D0)**2/10000.D0)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            Beta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            Aeta = Aeta + Beta
            Apsi = Apsi + Bpsi
            Apsi2 = Apsi2 + Bpsi2
            Veta = Veta + Beta**2
            Vpsi = Vpsi + Bpsi**2
            Vpsi2 = Vpsi2 + Bpsi2**2
        end do

        Aeta = Aeta/(Ntraj*sr)
        Veta = Veta/(Ntraj*sr**2)
        Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

        Apsi = Apsi/(Ntraj*sr**2)
        Vpsi = Vpsi/(Ntraj*sr**4)
        Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

        Apsi2 = Apsi2/(Ntraj*sr**2)
        Vpsi2 = Vpsi2/(Ntraj*sr**4)
        Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Apsi2+0.01348 = ", (Apsi2+0.01348D0), " +- ", Vpsi2

        !Check that both Qavg and S are within acceptable range
        call assertEquals(-0.01348D0, Apsi2, Vpsi2*2.D0, &
                          "psi_2 != -0.01348 for sr=1, hstar=0.15 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_eq_FENE_semimp(dt, Nsteps, Ntraj)
        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, Bq, Bs, Qavg, Vqavg, S, Serr, start_time, stop_time, b
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running zero-shear FENE Dumbbell (b=10) with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        k(:,:) = 0.D0

        delX = (/1.D0, 0.D0, 0.D0/)

        Qavg = 0.D0
        Vqavg = 0.D0
        S = 0.D0
        Serr = 0.D0

        b = 10.D0

        Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, b, 0.D0, dW)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product(delX, Q(:,i))**2/Ql2

            Qavg = Qavg + Ql2
            Vqavg = Vqavg + Ql
            S = S + 0.5D0*(3.D0*Bs - 1)
            Serr = Serr + 0.25D0*(9.D0*Bs**2 - 6.D0*Bs + 1.D0)
        end do

        Qavg = sqrt(Qavg/Ntraj)
        Vqavg = Vqavg/Ntraj
        Vqavg = sqrt((Qavg**2 - Vqavg**2)/(Ntraj-1))

        S = S/Ntraj
        Serr = Serr/Ntraj
        Serr = sqrt((Serr - S**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Q-(3*b)/(b+5) = ", (Qavg-sqrt(3*b/(b+5))), " +- ", Vqavg
        print *, "S = ", S, " +- ", Serr

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(3*b/(b+5)), Qavg, Vqavg*2.D0, &
                          "Qavg != sqrt(3*b/(b+5)) for sr=0, b=10 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, S, Serr*2.D0, &
                          "S != 0 for sr=0, b=10 FENE Dumbbell (semimp)")

        print *, ""
        print *, "Running zero-shear FENE Dumbbell (b=50) with semimp integrator test"

        b = 50.D0

        Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, b, 0.D0, dW)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product(delX, Q(:,i))**2/Ql2

            Qavg = Qavg + Ql2
            Vqavg = Vqavg + Ql
            S = S + 0.5D0*(3.D0*Bs - 1)
            Serr = Serr + 0.25D0*(9.D0*Bs**2 - 6.D0*Bs + 1.D0)
        end do

        Qavg = sqrt(Qavg/Ntraj)
        Vqavg = Vqavg/Ntraj
        Vqavg = sqrt((Qavg**2 - Vqavg**2)/(Ntraj-1))

        S = S/Ntraj
        Serr = Serr/Ntraj
        Serr = sqrt((Serr - S**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Q-(3*b)/(b+5) = ", (Qavg-sqrt(3*b/(b+5))), " +- ", Vqavg
        print *, "S = ", S, " +- ", Serr

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(3*b/(b+5)), Qavg, Vqavg*2.D0, &
                          "Qavg != sqrt(3*b/(b+5)) for sr=0, b=50 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, S, Serr*2.D0, &
                          "S != 0 for sr=0, b=50 FENE Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_semimp_euler_equal(dt, Nsteps, Ntraj)
        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, a, h, Beta, Bpsi, Bpsi2
        real*8, dimension(2) :: Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q_semimp, Q_euler

        print *, "Running semimp vs euler Hookean test, h*=0.15, sr=1"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q_semimp(3,1:Ntraj), Q_euler(3,1:Ntraj))

        sr = 1.D0;

        k(:,:) = 0.D0
        k(1,2) = 1.D0
        k = sr*k

        h = 0.15D0
        a = sqrt(PI)*h

        delX = (/1.D0, 0.D0, 0.D0/)

        Aeta = 0.D0
        Apsi = 0.D0
        Apsi2 = 0.D0
        Veta = 0.D0
        Vpsi = 0.D0
        Vpsi2 = 0.D0

        tau = 0.D0

        Q_semimp = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)
        Q_euler = Q_semimp

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q_euler, Q_semimp)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q_semimp(:,i) =  step(Q_semimp(:,i), k, dt, 0.D0, 10000.D0, a, dW)
                dW = Wiener_step(seed, dt)
                Q_euler(:,i) = step(Q_euler(:,i), k, dt, a, dW)
            end do
            !$OMP END DO

        end do
        !$OMP END parallel

        ! Measurement
        do i = 1,Ntraj
            Ql2 = Q_semimp(1,i)**2 + Q_semimp(2,i)**2 + Q_semimp(3,i)**2
            Ql = sqrt(Ql2)
            F(:) = (Ql - 0.D0)/(1.0D0-(Ql-0.D0)**2/10000.D0)*Q_semimp(:,i)/Ql
            tau(:,:) = dyadic_prod(Q_semimp(:,i), F)

            Beta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            Aeta(1) = Aeta(1) + Beta
            Apsi(1) = Apsi(1) + Bpsi
            Apsi2(1) = Apsi2(1) + Bpsi2
            Veta(1) = Veta(1) + Beta**2
            Vpsi(1) = Vpsi(1) + Bpsi**2
            Vpsi2(1) = Vpsi2(1) + Bpsi2**2

            Ql2 = Q_euler(1,i)**2 + Q_euler(2,i)**2 + Q_euler(3,i)**2
            Ql = sqrt(Ql2)
            F(:) = (Ql - 0.D0)/(1.0D0-(Ql-0.D0)**2/10000.D0)*Q_euler(:,i)/Ql
            tau(:,:) = dyadic_prod(Q_euler(:,i), F)

            Beta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            Aeta(2) = Aeta(2) + Beta
            Apsi(2) = Apsi(2) + Bpsi
            Apsi2(2) = Apsi2(2) + Bpsi2
            Veta(2) = Veta(2) + Beta**2
            Vpsi(2) = Vpsi(2) + Bpsi**2
            Vpsi2(2) = Vpsi2(2) + Bpsi2**2
        end do

        Aeta = Aeta/(Ntraj*sr)
        Veta = Veta/(Ntraj*sr**2)
        Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

        Apsi = Apsi/(Ntraj*sr**2)
        Vpsi = Vpsi/(Ntraj*sr**4)
        Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

        Apsi2 = Apsi2/(Ntraj*sr**2)
        Vpsi2 = Vpsi2/(Ntraj*sr**4)
        Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Aeta_semimp-Aeta_euler = ", (Aeta(2) - Aeta(1)), " +- ", (Veta(2) + Veta(1))
        print *, "Apsi_semimp-Apsi_euler = ", (Apsi(2) - Apsi(1)), " +- ", (Vpsi(2) + Vpsi(1))
        print *, "Apsi2_semimp-Apsi2_euler = ", (Apsi2(2) - Apsi2(1)), " +- ", (Vpsi2(2) + Vpsi2(1))

        !Check that both Qavg and S are within acceptable range
        call assertEquals(Aeta(1), Aeta(2), 2*sqrt(Veta(2)**2 + Veta(1)**2), &
                          "Aeta_semimp != Aeta_euler, h*=0.15, sr=1")
        call assertEquals(Apsi(1), Apsi(2), 2*sqrt(Vpsi(2)**2 + Vpsi(1)**2), &
                          "Apsi_semimp != Apsi_euler, h*=0.15, sr=1")
        call assertEquals(Apsi2(1), Apsi2(2), 2*sqrt(Vpsi2(2)**2 + Vpsi2(1)**2), &
                          "Apsi2_semimp != Apsi2_euler, h*=0.15, sr=1")

        deallocate(Q_euler, Q_semimp)
        print *, ""
    end subroutine

    subroutine test_FENE_HI_shear_semimp_vs_Kailash_code(dt, Nsteps, Ntraj)
        implicit none

        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i, shear_steps
        real*8 :: Ql, Ql2, Bq, Bs, Qavg, Vqavg, S, Serr, start_time, stop_time, b, sr_vals(5), h, a, sr, Q0
        real*8 :: Aeta(5), Veta(5), Apsi(5), Vpsi(5), Beta, Bpsi
        real*8 :: K_Aeta(5), K_Veta(5), K_Apsi(5), K_Vpsi(5)
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q_eq_VR, Q

        print *, "Comparison with Kailash's results for shear-flow FENE dumbbells with HI, VR"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj), Q_eq_VR(3,1:Ntraj))

        k(:,:) = 0.D0

        delX = (/1.D0, 0.D0, 0.D0/)

        Qavg = 0.D0
        Vqavg = 0.D0
        S = 0.D0
        Serr = 0.D0
        Aeta = 0.D0
        Apsi = 0.D0
        Veta = 0.D0
        Vpsi = 0.D0

        sr_vals = (/0.03D0, 0.1D0, 1.D0, 45.D0, 100.D0/)

        Q0 = 0.D0
        b = 100.D0
        h = 0.3D0
        a = sqrt(PI)*h

        Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)
        Q_eq_VR = Q

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()


        do shear_steps = 1,size(sr_vals)
            sr = sr_vals(shear_steps)
            Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)
            Q_eq_VR = Q

            !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR)
            !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8
            do steps = 1,Nsteps
                !$OMP DO
                do i = 1,Ntraj
                    dW = Wiener_step(seed, dt)
                    k(1,2) = sr
                    Q(:,i) =  step(Q(:,i), k, dt, Q0, b, a, dW)
                    k(1,2) = 0.D0
                    Q_eq_VR(:,i) =  step(Q_eq_VR(:,i), k, dt, Q0, b, a, dW)
                end do
                !$OMP END DO

            end do
            !$OMP END parallel

            ! Measurement
            do i = 1,Ntraj
                Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
                Ql = sqrt(Ql2)
                F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q(:,i)/Ql
                tau(:,:) = dyadic_prod(Q(:,i), F)

                Beta = tau(1,2)
                Bpsi = (tau(1,1) - tau(2,2))

                Ql2 = Q_eq_VR(1,i)**2 + Q_eq_VR(2,i)**2 + Q_eq_VR(3,i)**2
                Ql = sqrt(Ql2)
                F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q_eq_VR(:,i)/Ql
                tau(:,:) = dyadic_prod(Q_eq_VR(:,i), F)

                Beta = Beta - tau(1,2)
                Bpsi = Bpsi - (tau(1,1) - tau(2,2))
                Aeta(shear_steps) = Aeta(shear_steps) + Beta
                Apsi(shear_steps) = Apsi(shear_steps) + Bpsi
                Veta(shear_steps) = Veta(shear_steps) + Beta**2
                Vpsi(shear_steps) = Vpsi(shear_steps) + Bpsi**2
            end do

            Aeta(shear_steps) = Aeta(shear_steps)/(Ntraj*sr)
            Veta(shear_steps) = Veta(shear_steps)/(Ntraj*sr**2)
            Veta(shear_steps) = sqrt((Veta (shear_steps)- Aeta(shear_steps)**2)/(Ntraj-1))

            Apsi(shear_steps) = Apsi(shear_steps)/(Ntraj*sr**2)
            Vpsi(shear_steps) = Vpsi(shear_steps)/(Ntraj*sr**4)
            Vpsi(shear_steps) = sqrt((Vpsi(shear_steps) - Apsi(shear_steps)**2)/(Ntraj-1))

        end do


        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

!        print *, "Q-(3*b)/(b+5) = ", (Qavg-sqrt(3*b/(b+5))), " +- ", Vqavg
!        print *, "S = ", S, " +- ", Serr
!        print *, "Aeta values are", Aeta
!        print *, "Veta values are", Veta
!        print *, "Apsi values are", Apsi
!        print *, "Vpsi values are", Vpsi

        !Check that my results agree with Kailash's results
        K_Aeta = (/1.36355, 1.36426, 1.20931, 0.25452, 0.15777/)
        K_Veta = (/0.00616, 0.00482, 0.00107, 0.00022, 0.00014/)
        K_Apsi = (/3.53116, 3.64638, 2.79976, 0.09888, 0.03609/)
        K_Vpsi = (/0.03368, 0.0078, 0.00212, 0.00005, 0.00002/)
        do shear_steps = 1,size(sr_vals)
            print *, "sr=", sr_vals(shear_steps)
            print *, "My eta=", Aeta(shear_steps), "+-", Veta(shear_steps)
            print *, "Kailash eta=", K_Aeta(shear_steps), "+-", K_Veta(shear_steps)
            call assertEquals(K_Aeta(shear_steps), Aeta(shear_steps), &
                              2*sqrt(Veta(shear_steps)**2+K_Veta(shear_steps)**2), "Aeta")
            print *, ""
            print *, "My psi=", Apsi(shear_steps), "+-", Vpsi(shear_steps)
            print *, "Kailash psi=", K_Apsi(shear_steps), "+-", K_Vpsi(shear_steps)

            call assertEquals(K_Apsi(shear_steps), Apsi(shear_steps), &
                              2*sqrt(Vpsi(shear_steps)**2+K_Vpsi(shear_steps)**2), "Veta")
            print *, ""
        end do

    end subroutine


End Program
