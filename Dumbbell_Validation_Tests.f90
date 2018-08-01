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
        !call test_FF_zero_shear_viscosity(0.01D0, 2000, 100000)

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
        real*8 :: Ql, Ql2, start_time, stop_time, sr, dW(3), k(3,3)
        type(measured_variables) out_var
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running zero-shear hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        k(:,:) = 0.D0

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, 0.D0, dW)
            end do
            !$OMP END DO

        end do

        call measure_shear_no_VR(out_var, Q, 0.D0, 10000.D0, 0.D0, Ntraj)

        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Q-sqrt(3) = ", (out_var%Qavg-sqrt(3.D0)), " +- ", out_var%Vqavg
        print *, "S = ", out_var%S, " +- ", out_var%Serr

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(3.D0), out_var%Qavg, out_var%Vqavg*2.D0, &
                          "Qavg != sqrt(3) for sr=0 Hookean Dumbbell (semimp)")
        call assertEquals(0.D0, out_var%S, out_var%Serr*2.D0, &
                          "S != 0 for sr=0 Hookean Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_Hookean_viscosity_semimp(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, dW(3), k(3,3)
        type(measured_variables) out_var
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

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, 0.D0, dW)
            end do
            !$OMP END DO
        end do

        call measure_shear_no_VR(out_var, Q, 0.D0, 10000.D0, sr, Ntraj)

        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Aeta-0.63212 = ", (out_var%Aeta-0.63212D0), " +- ", out_var%Veta
        print *, "Apsi-0.52848 = ", (out_var%Apsi-0.52848D0), " +- ", out_var%Vpsi

        !Check that both Qavg and S are within acceptable range
        call assertEquals(0.63212D0, out_var%Aeta, out_var%Veta*2.D0, &
                          "eta != 0.63212 for sr=1 Hookean Dumbbell (semimp)")
        call assertEquals(0.52848D0, out_var%Apsi, out_var%Vpsi*2.D0, &
                          "psi_1 != 0.52848 for sr=1 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_Hookean_psi2_with_HI_semimp(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, dW(3), k(3,3), h, a
        type(measured_variables) out_var
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

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, a, dW)
            end do
            !$OMP END DO

        end do

        call measure_shear_no_VR(out_var, Q, 0.D0, 10000.D0, sr, Ntraj)

        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Apsi2+0.01348 = ", (out_var%Apsi2+0.01348D0), " +- ", out_var%Vpsi2

        !Check that both Qavg and S are within acceptable range
        call assertEquals(-0.01348D0, out_var%Apsi2, out_var%Vpsi2*2.D0, &
                          "psi_2 != -0.01348 for sr=1, hstar=0.15 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_Hookean_psi2_with_HI_semimp_2nd_method(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, dW(3), k(3,3), h, a
        type(measured_variables) out_var
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

        Q = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, 10000.D0, a, dW, 5.D0)
            end do
            !$OMP END DO

        end do

        call measure_shear_no_VR(out_var, Q, 0.D0, 10000.D0, sr, Ntraj)

        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Apsi2+0.01348 = ", (out_var%Apsi2+0.01348D0), " +- ", out_var%Vpsi2

        !Check that both Qavg and S are within acceptable range
        call assertEquals(-0.01348D0, out_var%Apsi2, out_var%Vpsi2*2.D0, &
                          "psi_2 != -0.01348 for sr=1, hstar=0.15 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_eq_FENE_semimp(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, dW(3), k(3,3), b
        type(measured_variables) out_var
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running zero-shear FENE Dumbbell (b=10) with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        k(:,:) = 0.D0

        b = 10.D0

        Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8
        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, b, 0.D0, dW)
            end do
            !$OMP END DO
        end do
        call measure_shear_no_VR(out_var, Q, 0.D0, b, 0.D0, Ntraj)
        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Q-(3*b)/(b+5) = ", (out_var%Qavg-sqrt(3*b/(b+5))), " +- ", out_var%Vqavg
        print *, "S = ", out_var%S, " +- ", out_var%Serr

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(3*b/(b+5)), out_var%Qavg, out_var%Vqavg*2.D0, &
                          "Qavg != sqrt(3*b/(b+5)) for sr=0, b=10 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, out_var%S, out_var%Serr*2.D0, &
                          "S != 0 for sr=0, b=10 FENE Dumbbell (semimp)")

        print *, ""
        print *, "Running zero-shear FENE Dumbbell (b=50) with semimp integrator test"

        b = 50.D0

        Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8
        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, b, 0.D0, dW)
            end do
            !$OMP END DO
        end do
        call measure_shear_no_VR(out_var, Q, 0.D0, b, 0.D0, Ntraj)
        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Q-(3*b)/(b+5) = ", (out_var%Qavg-sqrt(3*b/(b+5))), " +- ", out_var%Vqavg
        print *, "S = ", out_var%S, " +- ", out_var%Serr

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(3*b/(b+5)), out_var%Qavg, out_var%Vqavg*2.D0, &
                          "Qavg != sqrt(3*b/(b+5)) for sr=0, b=50 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, out_var%S, out_var%Serr*2.D0, &
                          "S != 0 for sr=0, b=50 FENE Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_semimp_euler_equal(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: Ql, Ql2, start_time, stop_time, sr, dW(3), k(3,3), h, a
        type(measured_variables) Vsemimp, Veuler
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

        Q_semimp = generate_Q_FF(0.D0, 10000.D0, Ntraj, seed, 10000)
        Q_euler = Q_semimp

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q_euler, Q_semimp, Veuler, Vsemimp)
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
        call measure_shear_no_VR(Veuler, Q_euler, 0.D0, 10000.D0, sr, Ntraj)
        call measure_shear_no_VR(Vsemimp, Q_semimp, 0.D0, 10000.D0, sr, Ntraj)
        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        print *, "Aeta_semimp-Aeta_euler = ", (Vsemimp%Aeta - Veuler%Aeta), &
                                    " +- ", (Vsemimp%Veta + Veuler%Veta)
        print *, "Apsi_semimp-Apsi_euler = ", (Vsemimp%Apsi - Veuler%Apsi), &
                                    " +- ", (Vsemimp%Vpsi + Veuler%Vpsi)
        print *, "Apsi2_semimp-Apsi2_euler = ", (Vsemimp%Apsi2 - Veuler%Apsi2), &
                                    " +- ", (Vsemimp%Vpsi2 + Veuler%Vpsi2)

        !Check that both Qavg and S are within acceptable range
        call assertEquals(Veuler%Aeta, Vsemimp%Aeta, 2*sqrt(Vsemimp%Veta**2 + Veuler%Veta**2), &
                          "Aeta_semimp != Aeta_euler, h*=0.15, sr=1")
        call assertEquals(Veuler%Apsi, Vsemimp%Apsi, 2*sqrt(Vsemimp%Vpsi**2 + Veuler%Vpsi**2), &
                          "Apsi_semimp != Apsi_euler, h*=0.15, sr=1")
        call assertEquals(Veuler%Apsi2, Vsemimp%Apsi2, 2*sqrt(Vsemimp%Vpsi2**2 + Veuler%Vpsi2**2), &
                          "Apsi2_semimp != Apsi2_euler, h*=0.15, sr=1")

        deallocate(Q_euler, Q_semimp)
        print *, ""
    end subroutine

    subroutine test_FENE_HI_shear_semimp_vs_Kailash_code(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i, s
        real*8 :: start_time, stop_time, b, sr_vals(5), h, a, sr, Q0
        type(measured_variables) out_var(5)
        real*8, dimension(3) :: dW
        real*8, dimension(3,3) :: k
        real*8, dimension(5) :: K_Aeta, K_Apsi, K_Veta, K_Vpsi
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q_eq_VR, Q

        print *, "Comparison with Kailash's results for shear-flow FENE dumbbells with HI, VR"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj), Q_eq_VR(3,1:Ntraj))

        k(:,:) = 0.D0

        sr_vals = (/0.03D0, 0.1D0, 1.D0, 45.D0, 100.D0/)

        Q0 = 0.D0
        b = 100.D0
        h = 0.3D0
        a = sqrt(PI)*h

        Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)
        Q_eq_VR = Q

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        do s = 1,size(sr_vals)
            sr = sr_vals(s)
            Q = generate_Q_FF(0.D0, b, Ntraj, seed, 10000)
            Q_eq_VR = Q

            !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, out_var)
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
            ! Measurement
            call measure_shear_with_VR(out_var(s), Q, Q_eq_VR, Q0, b, sr, Ntraj)
            !$OMP END parallel

        end do

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        !Check that my results agree with Kailash's results
        K_Aeta = (/1.36355, 1.36426, 1.20931, 0.25452, 0.15777/)
        K_Veta = (/0.00616, 0.00482, 0.00107, 0.00022, 0.00014/)
        K_Apsi = (/3.53116, 3.64638, 2.79976, 0.09888, 0.03609/)
        K_Vpsi = (/0.03368, 0.0078, 0.00212, 0.00005, 0.00002/)
        do s = 1,size(sr_vals)
            print *, "sr=", sr_vals(s)
            print *, "My eta=", out_var(s)%Aeta, "+-", out_var(s)%Veta
            print *, "Kailash eta=", K_Aeta(s), "+-", K_Veta(s)
            call assertEquals(K_Aeta(s), out_var(s)%Aeta, &
                              2*sqrt(out_var(s)%Veta**2+K_Veta(s)**2), "Aeta")
            print *, ""
            print *, "My psi=", out_var(s)%Apsi, "+-", out_var(s)%Vpsi
            print *, "Kailash psi=", K_Apsi(s), "+-", K_Vpsi(s)

            call assertEquals(K_Apsi(s), out_var(s)%Apsi, &
                              2*sqrt(out_var(s)%Vpsi**2+K_Vpsi(s)**2), "Veta")
            print *, ""
        end do

    end subroutine

    subroutine test_FF_zero_shear_viscosity(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i, s
        real*8 :: start_time, stop_time, alpha, h, a, sr, Q0, eta_ana, psi_ana, Q2_ana
        type(measured_variables) out_var
        real*8, dimension(3) :: dW
        real*8, dimension(3,3) :: k
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q, Q_eq_VR

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj), Q_eq_VR(3,1:Ntraj))

        k(:,:) = 0.D0

        sr = 0.001D0

        Q0 = 5.D0
        alpha = 0.1D0
        h = 0.D0
        a = sqrt(PI)*h

        print *, "Correct zero-shear (sr=", sr, ") viscosity for FF dumbbells"

        Q = generate_Q_FF(Q0, alpha, Ntraj, seed, 10000)
        Q_eq_VR = Q

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8
        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                k(1,2) = sr
                Q(:,i) =  step(Q(:,i), k, dt, Q0, alpha, a, dW)
                k(1,2) = 0.D0
                Q_eq_VR(:,i) =  step(Q_eq_VR(:,i), k, dt, Q0, alpha, a, dW)
            end do
            !$OMP END DO

        end do
        ! Measurement
        call measure_shear_with_VR(out_var, Q, Q_eq_VR, Q0, alpha, sr, Ntraj)
        !$OMP END parallel

        call cpu_time(stop_time)
        !$ stop_time = omp_get_wtime()
        print *, "Run time:", &
                stop_time - start_time, "seconds"

        Q2_ana = (((3.D0*alpha)/(alpha+5.D0)+6.D0*Q0**2)*(alpha/(alpha+3.D0))+Q0**4)/&
                  (alpha/(alpha+3.D0)+Q0**2)

        eta_ana = (1.D0/3.D0)*Q2_ana

        psi_ana = (1.D0/15.D0)*((((5.D0*alpha)/(alpha+7.D0)+15.D0*Q0**2)*&
                  (3.D0*alpha)/(alpha+5.D0)+15.D0*Q0**4)*(alpha/(alpha+3.D0))+Q0**6)/&
                  (alpha/(alpha+3.D0)+Q0**2)

        print *, "Q^2_analytical = ", Q2_ana
        print *, "eta_analytical = ", eta_ana
        print *, "psi_analytical = ", psi_ana
        print *, "Qavg = ", out_var%Qavg
        print *, "Aeta = ", out_var%Aeta
        print *, "Apsi = ", out_var%Apsi
        print *, "Qavg-sqrt(Q^2)_0 = " , out_var%Qavg-sqrt(Q2_ana), " +- ", out_var%Vqavg
        print *, "Aeta-eta_0 = ", (out_var%Aeta-eta_ana), " +- ", out_var%Veta
        print *, "Apsi-psi_0 = ", (out_var%Apsi-psi_ana), " +- ", out_var%Vpsi

        !Check that both Qavg and S are within acceptable range
        call assertEquals(sqrt(Q2_ana), out_var%Qavg, out_var%Vqavg*2.D0, &
                          "Qavg != sqrt(Q^2)_analytical for sr=0.001 FF dumbbell")
        call assertEquals(eta_ana, out_var%Aeta, out_var%Veta*2.D0, &
                          "Aeta != eta_analytical for sr=0.001 FF dumbbell")
        call assertEquals(psi_ana, out_var%Apsi, out_var%Vpsi*2.D0, &
                          "Apsi != psi_analytical for sr=0.001 FF dumbbell")

        print *, ""

    end subroutine

End Program
