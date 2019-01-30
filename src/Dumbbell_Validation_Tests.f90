Program Dumbbell_Validation_tests
    use Dumbbell_util
    use Generate_Initial_Conditions
    use Integrate
    use fruit
    use omp_lib
    implicit none

    call init_fruit
    call Validation_tests
    !call Unit_tests
    call fruit_summary
    call fruit_finalize

    contains

    subroutine Unit_tests()
        implicit none

        print *, "Running unit tests"

        call unit_test_locate()
        call unit_test_measure()
        call unit_test_measure_all_variables_no_VR()
        print *, ""

    end subroutine

    subroutine unit_test_locate()
        implicit none
        real*8 :: xvals(10), x, yvals(10)
        integer :: j

        xvals = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10/)
        yvals = xvals**2
        x = 4.6

        j = locate(xvals, x, 2)
        call assertEquals(4, j, "j!=4")

        j = locate(xvals, x, 3)
        call assertEquals(4, j, char(j))

        j = locate(xvals, 1.1D0, 4)
        call assertEquals(1, j)

        j = locate(xvals, 10.D0, 4)
        call assertEquals(6, j)

        j = locate(xvals, 5123.D0, 4)
        call assertEquals(6, j)

        j = locate(xvals, -423.D0, 4)
        call assertEquals(1, j)

        x = lin_interp_bs(xvals, yvals, 2.5D0)
        call assertEquals(6.5D0, x, 1D-8)

        x = poly_interp_bs(xvals, yvals, 2.1D0, 3)
        call assertEquals(2.1D0**2, x, 1D-8)

        x = poly_interp_bs(xvals, yvals, 5.6D0, 3)
        call assertEquals(5.6D0**2, x, 1D-8)

        x = poly_interp_bs(xvals, yvals, 11.D0, 3)
        call assertEquals(11.D0**2, x, 1D-8)

        x = poly_interp_bs(xvals, yvals, 0.5D0, 3)
        call assertEquals(0.5D0**2, x, 1D-8)

        x = lin_interp_bs(xvals, yvals, 5.5D0)
        call assertNotEquals(5.5D0**2, x, 1D-8)

        call assertEquals(lin_interp_bs(xvals,yvals, 5.41D0), poly_interp_bs(xvals, yvals, 5.41D0, 2), 1D-8)

        yvals = xvals**5

        x = poly_interp_bs(xvals, yvals, 1.1D0, 6)
        call assertEquals(1.1D0**5, x, 1D-8)

        yvals = xvals**3

        x = poly_interp_bs(xvals, yvals, 0.9D0, 4)
        call assertEquals(0.9D0**3, x, 1D-8)

    end subroutine

    subroutine unit_test_measure()
        implicit none
        real*8, dimension(3,30) :: Q
        integer*8 :: N
        real*8 :: avgA, avgB, varA, varB, covAB
        N = 30

        Q = reshape((/  5.6029,   14.5768,   -1.3101, &
                        6.6505,   13.7525,    3.4800, &
                        9.5981,   -5.8658,   10.8930, &
                       13.0620,    2.4046,    8.3193, &
                      -10.8488,   -6.1034,    9.4963, &
                       -6.1441,  -13.3344,    5.4456, &
                        6.2891,   12.1568,    7.6120, &
                       10.6699,   -0.5442,  -11.4562, &
                        6.7870,  -14.1031,   -0.2352, &
                       12.9328,    7.7425,    4.2637, &
                        1.2288,  -13.1725,   -8.3873, &
                        5.6379,   14.6187,   -0.4575, &
                       -9.1385,   -2.4393,   12.5012, &
                       -6.4316,    4.8980,   13.4052, &
                       11.5259,   10.0613,    3.3339, &
                      -15.0162,   -2.1274,    3.8903, &
                      -15.3353,    1.0163,   -2.9702, &
                       -5.6396,    6.2548,  -13.2135, &
                        5.0111,    5.7722,  -13.6645, &
                      -15.3403,   -1.5905,    2.8346, &
                       -0.8900,   14.6681,    5.4176, &
                       10.9586,    7.8865,   -7.9505, &
                        6.0137,   -7.7398,   12.2107, &
                        2.5731,   13.3188,   -7.8074, &
                       15.2916,   -1.7554,   -2.8481, &
                      -15.6276,    0.7459,   -0.7577, &
                       14.5073,   -4.8661,   -3.3425, &
                       -1.6233,   -6.3572,   14.2364, &
                       -5.1473,   -9.5874,   11.2803, &
                      -14.6866,   -2.0634,   -5.0104/), &
                                            shape(Q))

        print *, ""
        print *, "performing parallel measure"
        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(avgA, avgB, varA, varB, covAB)
        call measure(N, Q, avgA, avgB, varA, varB, covAB)
        !$OMP END PARALLEL
        print *, "writing parallel results"
        print *, "avgA, avgB, varA, varB, covAB"
        write(*,*) avgA, avgB, varA, varB, covAB

        call assertEquals(20.694007550333335D0, avgA, 1D-5, &
                "avgA != matlab average with parallelisation")
        call assertEquals(22.006523030666667D0, avgB, 1D-5, &
                "avgB != matlab average with parallelisation")
        call assertEquals(3.059346682411035D3, varA, 1D-5, &
                "varA != matlab average with parallelisation")
        call assertEquals(2.276684944252046D4, varB, 1D-2, &
                "varB != matlab average with parallelisation")
        call assertEquals(-0.152114655473403D4, covAB, 1D-2, &
                "covAB != matlab average with parallelisation")

        print *, ""
        print *, "performing serial measure"
        call measure(N, Q, avgA, avgB, varA, varB, covAB)
        print *, "writing serial results"
        print *, "avgA, avgB, varA, varB, covAB"
        write(*,*) avgA, avgB, varA, varB, covAB

        call assertEquals(20.694007550333335D0, avgA, 1D-5, &
                "avgA != matlab average without parallelisation")
        call assertEquals(22.006523030666667D0, avgB, 1D-5, &
                "avgB != matlab average without parallelisation")
        call assertEquals(3.059346682411035D3, varA, 1D-5, &
                "varA != matlab average without parallelisation")
        call assertEquals(2.276684944252046D4, varB, 1D-2, &
                "varB != matlab average without parallelisation")
        call assertEquals(-0.152114655473403D4, covAB, 1D-2, &
                "covAB != matlab average without parallelisation")

    end subroutine

    subroutine unit_test_measure_all_variables_no_VR()
        implicit none
        real*8, dimension(3,30) :: Q
        type(measured_variables) :: output_variables
        real*8 :: sigma, alpha, sr, time1, time2
        integer*8 :: N

        N = 30
        sigma = 15.665239548400001D0
        alpha = 2.4539973011D-4
        sr = 0.005D0

        Q = reshape((/ 5.543438785899431D0, -10.402170853586162D0,  10.320844780622027D0, &
                       9.702769059757708D0, -12.175639414571998D0,  -1.817227268324901D0, &
                     -11.893173302153242D0,  -7.511011366390604D0,  -6.905269745987831D0, &
                     -10.344087833814578D0,  11.754784230629834D0,   0.324780392859595D0, &
                      13.101854738650246D0,  -6.641515950155460D0,  -5.450300333615164D0, &
                      -7.101255874171087D0,  12.661871035957237D0,  -5.879055890488128D0, &
                       3.785195377055119D0,  11.351672697083293D0, -10.126255368223024D0, &
                     -11.162834599985779D0, -10.236077181086715D0,  -4.025269634663749D0, &
                      10.057392688859574D0,  -5.794391079197713D0, -10.543171863658531D0, &
                      -0.274177517771730D0,   0.524480049998568D0, -15.661613174645945D0, &
                       4.941799194604504D0, -11.612651233049332D0,  -9.296095567076499D0, &
                      10.611689358556761D0,   9.924893250892564D0,  -5.845792577448453D0, &
                      10.625826404246803D0,  -4.345548335041154D0, -10.644612279633108D0, &
                       2.312578859666597D0, -15.487584261965582D0,  -0.071327788705253D0, &
                     -12.503515453415714D0,   3.817103615995090D0,   8.609379306338216D0, &
                      -1.184721776828417D0,  13.877519207069220D0,  7.170121076041205D0, &
                      -2.499825893888357D0, -14.732985482146356D0, -4.733087231236928D0, &
                     -13.751486472138488D0,   7.370329339771016D0,  -1.569397791691287D0, &
                      13.323619968103518D0,   4.371114309155102D0,   7.016397013747762D0, &
                      10.658988187230912D0,   0.640127912302644D0,  11.480046826966507D0, &
                      -6.033269430176779D0,   8.233677403204039D0, -11.866373649718588D0, &
                       0.341037196382751D0, -12.896816362812922D0,   8.869204274660524D0, &
                     -13.403050184244714D0,   4.179634316157828D0,   6.918727900159459D0, &
                       1.179986624631155D0, -11.357887857229105D0, -10.703563177325268D0, &
                      -3.199471133801409D0,   8.608893947892446D0, -12.693348259935862D0, &
                       6.194902552624718D0,  11.688089424576185D0,  -8.363762075217574D0, &
                      13.549565203435876D0,   2.222923346655979D0,  -7.529195190194018D0, &
                     -15.539768712321557D0,  -1.465318688866009D0,  -1.252789137341028D0, &
                     -13.325359010157820D0,  -5.157877109647665D0,   6.453348518733802D0, &
                      -6.347529139979506D0,  10.254582791728707D0,  -9.988882460573381D0/), &
                    shape(Q))

        !$ time1 = omp_get_wtime()
        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(output_variables)
        call measure_shear_no_VR(output_variables, Q, sigma, alpha, sr, N)
        !$OMP END PARALLEL
        !$ time2 = omp_get_wtime()

        print *, ""
        print *, "time taken is ", time2-time1, "seconds"

        print *, "Average eta, error eta"
        print *, output_variables%AvgEta, output_variables%ErrEta

        print *, "Average psi, error psi"
        print *, output_variables%AvgPsi, output_variables%ErrPsi

        print *, "Average psi2, error psi2"
        print *, output_variables%AvgPsi2, output_variables%ErrPsi2

        print *, "Average Ql, error Ql"
        print *, output_variables%AvgQ, output_variables%ErrQ

        print *, "Average S, error S"
        print *, output_variables%S, output_variables%ErrS

        print *, "Average chiG, error ChiG"
        print *, output_variables%AvgChiG, output_variables%ErrChiG

        print *, "Average chiTau, error ChiTau"
        print *, output_variables%AvgChiTau, output_variables%ErrChiTau

        call assertEquals(-74.305185195104102D0, output_variables%AvgEta, 1D-5, &
                "Average eta != matlab average")
        call assertEquals(64.5448809703651D0, output_variables%ErrEta, 1D-5, &
                "Error eta != matlab average")
        call assertEquals(22940.7503883371D0, output_variables%AvgPsi, 1D-5, &
                "Average Psi1 != matlab average")
        call assertEquals(1.700507684704652D4, output_variables%ErrPsi, 1D-5, &
                "Error Psi1 != matlab average")
        call assertEquals(1.444778929686133D03, output_variables%AvgPsi2, 1D-5, &
                "Error Ql != matlab average")
        call assertEquals(7.814469868788609D03, output_variables%ErrPsi2, 1D-5, &
                "Error Ql != matlab average")
        call assertEquals(15.665802591181151D0, output_variables%AvgQ, 1D-5, &
                "Average Ql != matlab average")
        call assertEquals(0.001790858468768D0, output_variables%ErrQ, 1D-5, &
                "Error Ql != matlab average")



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
        print *, "Running Validation tests"
        print *, "We expect about 5% of tests to fail by chance, or ~1/20"
        print *, ""
        call test_eq_hookean_semimp(0.001D0, 5000, 10000)
        call test_Hookean_viscosity_semimp(0.001D0, 100000)
        !call test_Hookean_psi2_with_HI_semimp(0.001D0, 5000, 10000)
        call test_eq_FENE_semimp(0.001D0, 5000, 10000)
        !call test_semimp_euler_equal(0.001D0, 500, 10000)
        !call test_FENE_HI_shear_semimp_vs_Kailash_code(0.001D0, 10000, 10000)
        call test_FF_zero_shear_viscosity(0.002D0, 10000, 5000)
        !call test_lookup_method(0.01D0, 1000, 10000)

        !p-test on likelihood of k or more 'failures' in n trials
        call get_failed_count(k)
        call get_total_count(n)
        k_r = k
        n_r = n
        !Incomplete beta function gives cumulative binomial probability
        !https://en.wikipedia.org/wiki/Binomial_distribution#Cumulative_distribution_function
        prob = 1.D0 - betai(n_r-k_r+1.D0, k_r, 0.95D0)

        print "(A, F5.1, A, I2, A, I3, A)", &
         "There is a", prob*100.D0, "% chance of getting", k, " or more failures in", &
          n, " tests, assuming the underlying failure probability is 5% (2 sigma)"
        print *, ""

    end subroutine

    subroutine test_eq_hookean_semimp(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: start_time, stop_time, dW(3), k(3,3)
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

        print *, "Q-sqrt(3) = ", (out_var%AvgQ-sqrt(3.D0)), " +- ", out_var%ErrQ
        print *, "S = ", out_var%S, " +- ", out_var%ErrS
        print *, "chi_G = ", (out_var%AvgChiG), " +- ", out_var%ErrChiG
        print *, "chi_tau = ", (out_var%AvgChiTau), " +- ", out_var%ErrChiTau

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(sqrt(3.D0), out_var%AvgQ, out_var%ErrQ*2.D0, &
                          "AvgQ != sqrt(3) for sr=0 Hookean Dumbbell (semimp)")
        call assertEquals(0.D0, out_var%S, out_var%ErrS*2.D0, &
                          "S != 0 for sr=0 Hookean Dumbbell (semimp)")



        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_Hookean_viscosity_semimp(dt, Ntraj)
        implicit none
        integer*8, intent(in) :: Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i, Nsteps
        real*8 :: start_time, stop_time, sr, dW(3), k(3,3)
        type(measured_variables) out_var
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q

        print *, "Running sr=1 hookean dumbbell with semimp integrator test"

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj))

        Nsteps = nint(1.D0/dt)

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

        print *, "AvgEta-0.63212 = ", (out_var%AvgEta-0.63212D0), " +- ", out_var%ErrEta
        print *, "AvgPsi-0.52848 = ", (out_var%AvgPsi-0.52848D0), " +- ", out_var%ErrPsi
        print *, "chi_G = ", (out_var%AvgChiG), " +- ", out_var%ErrChiG
        print *, "chi_tau = ", (out_var%AvgChiTau), " +- ", out_var%ErrChiTau


        !Check that both AvgQ and S are within acceptable range
        call assertEquals(0.63212D0, out_var%AvgEta, out_var%ErrEta*2.D0, &
                          "eta != 0.63212 for sr=1 Hookean Dumbbell (semimp)")
        call assertEquals(0.52848D0, out_var%AvgPsi, out_var%ErrPsi*2.D0, &
                          "psi_1 != 0.52848 for sr=1 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_Hookean_psi2_with_HI_semimp(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: start_time, stop_time, sr, dW(3), k(3,3), h, a
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

        print *, "AvgPsi2+0.01348 = ", (out_var%AvgPsi2+0.01348D0), " +- ", out_var%ErrPsi2

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(-0.01348D0, out_var%AvgPsi2, out_var%ErrPsi2*2.D0, &
                          "psi_2 != -0.01348 for sr=1, hstar=0.15 Hookean Dumbbell (semimp)")

        deallocate(Q)
        print *, ""
    end subroutine

    subroutine test_eq_FENE_semimp(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: start_time, stop_time, dW(3), k(3,3), b
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

        print *, "Q-(3*b)/(b+5) = ", (out_var%AvgQ-sqrt(3*b/(b+5))), " +- ", out_var%ErrQ
        print *, "S = ", out_var%S, " +- ", out_var%ErrS
        print *, "chi_G = ", (out_var%AvgChiG), " +- ", out_var%ErrChiG
        print *, "chi_tau = ", (out_var%AvgChiTau), " +- ", out_var%ErrChiTau

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(sqrt(3*b/(b+5)), out_var%AvgQ, out_var%ErrQ*2.D0, &
                          "AvgQ != sqrt(3*b/(b+5)) for sr=0, b=10 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, out_var%S, out_var%ErrS*2.D0, &
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

        print *, "Q-(3*b)/(b+5) = ", (out_var%AvgQ-sqrt(3*b/(b+5))), " +- ", out_var%ErrQ
        print *, "S = ", out_var%S, " +- ", out_var%ErrS
        print *, "chi_G = ", (out_var%AvgChiG), " +- ", out_var%ErrChiG
        print *, "chi_tau = ", (out_var%AvgChiTau), " +- ", out_var%ErrChiTau

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(sqrt(3*b/(b+5)), out_var%AvgQ, out_var%ErrQ*2.D0, &
                          "AvgQ != sqrt(3*b/(b+5)) for sr=0, b=50 FENE Dumbbell (semimp)")
        call assertEquals(0.D0, out_var%S, out_var%ErrS*2.D0, &
                          "S != 0 for sr=0, b=50 FENE Dumbbell (semimp)")

        deallocate(Q)

        print *, ""
    end subroutine

    subroutine test_semimp_euler_equal(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: start_time, stop_time, sr, dW(3), k(3,3), h, a
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

        print *, "Aeta_semimp-Aeta_euler = ", (Vsemimp%AvgEta - Veuler%AvgEta), &
                                    " +- ", (Vsemimp%ErrEta + Veuler%ErrEta)
        print *, "Apsi_semimp-Apsi_euler = ", (Vsemimp%AvgPsi - Veuler%AvgPsi), &
                                    " +- ", (Vsemimp%ErrPsi + Veuler%ErrPsi)
        print *, "Apsi2_semimp-Apsi2_euler = ", (Vsemimp%AvgPsi2 - Veuler%AvgPsi2), &
                                    " +- ", (Vsemimp%ErrPsi2 + Veuler%ErrPsi2)

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(Veuler%AvgEta, Vsemimp%AvgEta, 2*sqrt(Vsemimp%ErrEta**2 + Veuler%ErrEta**2), &
                          "Aeta_semimp != Aeta_euler, h*=0.15, sr=1")
        call assertEquals(Veuler%AvgPsi, Vsemimp%AvgPsi, 2*sqrt(Vsemimp%ErrPsi**2 + Veuler%ErrPsi**2), &
                          "Apsi_semimp != Apsi_euler, h*=0.15, sr=1")
        call assertEquals(Veuler%AvgPsi2, Vsemimp%AvgPsi2, 2*sqrt(Vsemimp%ErrPsi2**2 + Veuler%ErrPsi2**2), &
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
            k(1,2) = sr

            !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, out_var)
            !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8
            do steps = 1,Nsteps
                !$OMP DO
                do i = 1,Ntraj
                    dW = Wiener_step(seed, dt)
                    Q(:,i) =  step(Q(:,i), k, dt, Q0, b, a, dW)
                    Q_eq_VR(:,i) =  step(Q_eq_VR(:,i), dt, Q0, b, a, dW)
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
            print *, "My eta=", out_var(s)%AvgEta, "+-", out_var(s)%ErrEta
            print *, "Kailash eta=", K_Aeta(s), "+-", K_Veta(s)
            call assertEquals(K_Aeta(s), out_var(s)%AvgEta, &
                              2*sqrt(out_var(s)%ErrEta**2+K_Veta(s)**2), "AvgEta")
            print *, ""
            print *, "My psi=", out_var(s)%AvgPsi, "+-", out_var(s)%ErrPsi
            print *, "Kailash psi=", K_Apsi(s), "+-", K_Vpsi(s)

            call assertEquals(K_Apsi(s), out_var(s)%AvgPsi, &
                              2*sqrt(out_var(s)%ErrPsi**2+K_Vpsi(s)**2), "AvgPsi")
            print *, ""
        end do

    end subroutine

    subroutine test_FF_zero_shear_viscosity(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
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
        sr = 0.0001D0
        k(1,2) = sr

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
                Q(:,i) =  step(Q(:,i), k, dt, Q0, alpha, a, dW)
                Q_eq_VR(:,i) =  step(Q_eq_VR(:,i), dt, Q0, alpha, a, dW)
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

        print *, "Q_analytical = ", sqrt(Q2_ana)
        print *, "eta_analytical = ", eta_ana
        print *, "psi_analytical = ", psi_ana
        print *, "AvgQ = ", out_var%AvgQ
        print *, "AvgEta = ", out_var%AvgEta
        print *, "AvgPsi = ", out_var%AvgPsi
        print *, "AvgQ-sqrt(Q^2)_0 = " , out_var%AvgQ-sqrt(Q2_ana), " +- ", out_var%ErrQ
        print *, "AvgEta-eta_0 = ", (out_var%AvgEta-eta_ana), " +- ", out_var%ErrEta
        print *, "AvgPsi-psi_0 = ", (out_var%AvgPsi-psi_ana), " +- ", out_var%ErrPsi
        print *, "chi_G = ", (out_var%AvgChiG), " +- ", out_var%ErrChiG
        print *, "chi_tau = ", (out_var%AvgChiTau), " +- ", out_var%ErrChiTau

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(sqrt(Q2_ana), out_var%AvgQ, out_var%ErrQ*2.D0, &
                          "AvgQ != sqrt(Q^2)_analytical for sr=0.001 FF dumbbell")
        call assertEquals(eta_ana, out_var%AvgEta, out_var%ErrEta*2.D0, &
                          "AvgEta != eta_analytical for sr=0.001 FF dumbbell")
        call assertEquals(psi_ana, out_var%AvgPsi, out_var%ErrPsi*2.D0, &
                          "AvgPsi != psi_analytical for sr=0.001 FF dumbbell")

        print *, ""

    end subroutine

    subroutine test_lookup_method(dt, Nsteps, Ntraj)
        implicit none
        integer*8, intent(in) :: Nsteps, Ntraj
        real*8, intent(in) :: dt
        integer*8 :: steps, time(1:8), seed, i
        real*8 :: start_time, stop_time, alpha, h, a, sr, Q0, eta_ana, psi_ana, Q2_ana
        type(measured_variables) out_var
        real*8, dimension(3) :: dW
        real*8, dimension(3,3) :: k
        !large arrays must be declared allocatable so they go into heap, otherwise
        !OpenMP threads run out of memory
        real*8, dimension(:, :), allocatable :: Q, Q_eq_VR, storage_vals
        real*8, dimension(:), allocatable :: Yvals, Qvals
        integer :: Nvals

        Nvals = 1001

        call date_and_time(values=time)
        seed = time(8)*100 + time(7)*10

        allocate(Q(3,1:Ntraj), Q_eq_VR(3,1:Ntraj))
        allocate(Yvals(2*Nvals), Qvals(2*Nvals), storage_vals(2*Nvals, 2))

        k(:,:) = 0.D0
        sr = 0.0001D0
        k(1,2) = sr

        Q0 = 5.D0
        alpha = 0.1D0
        h = 0.D0
        a = sqrt(PI)*h

        print *, "Correct zero-shear (sr=", sr, ") viscosity for FF dumbbells, lookup method"

        Q = generate_Q_FF(Q0, alpha, Ntraj, seed, 10000)
        Q_eq_VR = Q

        call cpu_time(start_time)
        !$ start_time = omp_get_wtime()

        storage_vals = generate_lookup_table(real(dt*alpha/4.D0,kind=quad), real(Q0, kind=quad), real(alpha, kind=quad), Nvals)
        Yvals = storage_vals(:,1)
        Qvals = storage_vals(:,2)

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8
        do steps = 1,Nsteps
            !$OMP DO
            do i = 1,Ntraj
                dW = Wiener_step(seed, dt)
                Q(:,i) =  step(Q(:,i), k, dt, Q0, alpha, a, dW, Yvals, Qvals)
                Q_eq_VR(:,i) =  step(Q_eq_VR(:,i), dt, Q0, alpha, a, dW, Yvals, Qvals)
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

        print *, "Q_ana = ", sqrt(Q2_ana)
        print *, "eta_analytical = ", eta_ana
        print *, "psi_analytical = ", psi_ana
        print *, "AvgQ = ", out_var%AvgQ
        print *, "AvgEta = ", out_var%AvgEta
        print *, "AvgPsi = ", out_var%AvgPsi
        print *, "AvgQ-sqrt(Q^2)_0 = " , out_var%AvgQ-sqrt(Q2_ana), " +- ", out_var%ErrQ
        print *, "AvgEta-eta_0 = ", (out_var%AvgEta-eta_ana), " +- ", out_var%ErrEta
        print *, "AvgPsi-psi_0 = ", (out_var%AvgPsi-psi_ana), " +- ", out_var%ErrPsi

        !Check that both AvgQ and S are within acceptable range
        call assertEquals(sqrt(Q2_ana), out_var%AvgQ, out_var%ErrQ*2.D0, &
                          "AvgQ != sqrt(Q^2)_analytical for sr=0.001 FF dumbbell")
        call assertEquals(eta_ana, out_var%AvgEta, out_var%ErrEta*2.D0, &
                          "AvgEta != eta_analytical for sr=0.001 FF dumbbell")
        call assertEquals(psi_ana, out_var%AvgPsi, out_var%ErrPsi*2.D0, &
                          "AvgPsi != psi_analytical for sr=0.001 FF dumbbell")

        print *, ""

    end subroutine

End Program
