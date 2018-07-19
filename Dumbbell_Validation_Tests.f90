Program Dumbbell_Validation_tests
    use Random_Numbers
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
        call test_eq_hookean_semimp
        call test_Hookean_viscosity
        call test_Hookean_psi2_with_HI
        call test_eq_FENE_semimp

    end subroutine

    subroutine test_eq_hookean_semimp()
        implicit none

        integer*8 :: Nsteps, steps, time(1:8), seed, i, Ntraj
        real*8 :: dt, Ql, Ql2, Bq, Bs, Qavg, Vqavg, S, Serr, start_time, stop_time
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, delT, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running zero-shear hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        dt = 0.01D0
        Nsteps = 1000
        Ntraj = 100000

        allocate(Q(3,1:Ntraj))

        k(:,:) = 0.D0

        delT(:,:) = 0.D0
        delT(1,1) = 1.D0
        delT(2,2) = 1.D0
        delT(3,3) = 1.D0

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
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, delT, 10000.D0, 0.D0, dW)
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
        call assertEquals(sqrt(3.D0), Qavg, Vqavg*3.D0, &
                          "Qavg != sqrt(3) for sr=0 Hookean Dumbbell (semimp)")
        call assertEquals(0.D0, S, Serr*3.D0, &
                          "S != 0 for sr=0 Hookean Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_Hookean_viscosity()
        implicit none

        integer*8 :: Nsteps, steps, time(1:8), seed, i, Ntraj
        real*8 :: dt, Ql, Ql2, start_time, stop_time, sr
        real*8 :: Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2, Beta, Bpsi, Bpsi2
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, delT, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running sr=1 hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        dt = 0.01D0
        Nsteps = 100
        Ntraj = 100000

        allocate(Q(3,1:Ntraj))

        sr = 1.D0;

        k(:,:) = 0.D0
        k(1,2) = 1.D0
        k = sr*k

        delT(:,:) = 0.D0
        delT(1,1) = 1.D0
        delT(2,2) = 1.D0
        delT(3,3) = 1.D0

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
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, delT, 10000.D0, 0.D0, dW)
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
        call assertEquals(0.63212D0, Aeta, Veta*3.D0, &
                          "eta != 0.63212 for sr=1 Hookean Dumbbell (semimp)")
        call assertEquals(0.52848D0, Apsi, Vpsi*3.D0, &
                          "psi_1 != 0.52848 for sr=1 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_Hookean_psi2_with_HI()
        implicit none

        integer*8 :: Nsteps, steps, time(1:8), seed, i, Ntraj
        real*8 :: dt, Ql, Ql2, start_time, stop_time, sr, h, a
        real*8 :: Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2, Beta, Bpsi, Bpsi2
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, delT, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running sr=1, hstar=0.15 hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        dt = 0.01D0
        Nsteps = 1000
        Ntraj = 100000

        allocate(Q(3,1:Ntraj))

        sr = 1.D0;

        h = 0.15D0
        a = h*sqrt(PI)

        k(:,:) = 0.D0
        k(1,2) = 1.D0
        k = sr*k

        delT(:,:) = 0.D0
        delT(1,1) = 1.D0
        delT(2,2) = 1.D0
        delT(3,3) = 1.D0

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
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, delT, 10000.D0, a, dW)
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
        call assertEquals(-0.01348D0, Apsi2, Vpsi2*3.D0, &
                          "psi_2 != -0.01348 for sr=1, hstar=0.15 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_eq_FENE_semimp()
        implicit none

        integer*8 :: Nsteps, steps, time(1:8), seed, i, Ntraj
        real*8 :: dt, Ql, Ql2, Bq, Bs, Qavg, Vqavg, S, Serr, start_time, stop_time, b
        real*8, dimension(3) :: F, dW, Qtemp, delX
        real*8, dimension(3,3) :: k, delT, tau
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running zero-shear FENE Dumbbell (b=10) with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        dt = 0.05D0
        Nsteps = 1000
        Ntraj = 100000

        allocate(Q(3,1:Ntraj))

        k(:,:) = 0.D0

        delT(:,:) = 0.D0
        delT(1,1) = 1.D0
        delT(2,2) = 1.D0
        delT(3,3) = 1.D0

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
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, delT, b, 0.D0, dW)
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
        call assertEquals(sqrt(3*b/(b+5)), Qavg, Vqavg*3.D0, &
                          "Qavg != sqrt(3*b/(b+5)) for sr=0, b=10 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, S, Serr*3.D0, &
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
                Q(:,i) =  step(Q(:,i), k, dt, 0.D0, delT, b, 0.D0, dW)
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
        call assertEquals(sqrt(3*b/(b+5)), Qavg, Vqavg*3.D0, &
                          "Qavg != sqrt(3*b/(b+5)) for sr=0, b=50 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, S, Serr*3.D0, &
                          "S != 0 for sr=0, b=50 FENE Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine
End Program
