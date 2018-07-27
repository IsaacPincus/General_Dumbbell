
program General_Dumbbell
    use Integrate
    use Generate_Initial_Conditions
    use Dumbbell_util
    use omp_lib
    implicit none

    integer*8 :: Nblocks, Nsteps, steps, block, time(1:8), seed, i, Ntraj, VR_opt, dist_opt
    integer :: NTimeSteps, timestep, trajectories, Ndtwidths
    real*8 :: sr, b, h, a, Q0, Nrelax_times, dt, Ql, Ql2, F(3), dW(3), Qtemp(3), Bq, Bs, delX(3)
    real*8 :: time1, time2, B_eta, Bpsi, Bpsi2, output_delay
    real*8 :: Qavg, Vqavg, S, Serr, Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2
    real*8, dimension(3,3) :: k, tau
    !large arrays must be declared allocatable so they go into heap, otherwise
    !OpenMP threads run out of memory
    real*8, dimension(:, :), allocatable :: Q, Q_eq_VR
    real*8, dimension(:), allocatable :: timestepwidths

    call date_and_time(values=time)
    seed = time(8)*100 + time(7)*10

    open(unit=31, file='inputparameters.inp')
    open(unit=30, file='timestepdata.inp')
    open(unit=32, file='options.inp')
    open(unit=25, file='Q_dist_output.dat')

    open(unit=20, file='Ql.dat')
    open(unit=21, file='S.dat')
    open(unit=22, file='eta.dat')
    open(unit=23, file='psi.dat')
    open(unit=24, file='psi2.dat')

    read (32, *) VR_opt, dist_opt
    read (31, *) sr, b, h, Q0
    read (30, *) Ntraj, Ndtwidths, Nrelax_times, output_delay
    allocate(timestepwidths(Ndtwidths))
    do i=1,Ndtwidths
        read(30, *) timestepwidths(i)
    end do

    allocate(Q(3,1:Ntraj))
    if (VR_opt.eq.1) then
        allocate(Q_eq_VR(3,1:Ntraj))
    end if

    a = h*sqrt(PI)

    !For shear flow
    k(:,:) = 0.D0
    k(1,2) = 1.D0
    k = sr*k

    delX = (/1.D0, 0.D0, 0.D0/)

    looptimesteps: do timestep=1,Ndtwidths

        dt = timestepwidths(timestep)
        NtimeSteps = int(Nrelax_times/dt)

        !Also works for Fraenkel, FENE and Hookean springs
        !Reduce final input for faster computation of equilibrium dist
        Q(:,:) = generate_Q_FF(Q0, b, Ntraj, seed, 10000)

        call initialise_output_files_timestep()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, VR_opt), &
        !$OMP& SHARED(Qavg, Vqavg, S, Serr, Aeta, Veta, Apsi, Vpsi, Apsi2, Vpsi2)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        if (VR_opt.eq.0) then
            do steps=1,Ntimesteps
                k(1,2) = sr
                !$OMP DO schedule(dynamic, 100)
                do i=1,Ntraj
                    dW = Wiener_step(seed, dt)
                    Q(:,i) =  step(Q(:,i), k, dt, Q0, b, a, dW)
                end do
                !$OMP END DO
                if (steps.eq.1.or.(mod(steps,nint(output_delay/dt)).eq.0).or.steps.eq.NTimeSteps) then
                    call write_data_to_files_no_VR()
                end if
            end do

        elseif (VR_opt.eq.1) then

            !$OMP single
            Q_eq_VR(:,:) = Q(:,:)
            !$OMP end single

            do steps=1,Ntimesteps

                !$OMP DO schedule(dynamic, 100)
                do i=1,Ntraj
                    dW = Wiener_step(seed, dt)
                    k(1,2) = sr
                    Q(:,i) =  step(Q(:,i), k, dt, Q0, b, a, dW)
                    k(1,2) = 0.D0
                    Q_eq_VR(:,i) = step(Q_eq_VR(:,i), k, dt, Q0, b, a, dW)
                end do
                !$OMP END DO
                if (steps.eq.1.or.(mod(steps,nint(output_delay/dt)).eq.0).or.steps.eq.NTimeSteps) then
                    call write_data_to_files_with_VR()
                end if
            end do

        else
            print *, "Variance Reduction option not set, place in options.inp file"
        end if

        !$OMP END PARALLEL

        if (dist_opt.eq.1) then
            25 format(2(E12.5,4X), E12.5)
            write(25, *) 'Timestep width is: ', dt
            do i=1,Ntraj
                write(25, 25) Q(:,i)
            end do
        end if

    end do looptimesteps

    !close(unit=21)
    close(unit=20)
    close(unit=21)
    close(unit=22)
    close(unit=23)
    close(unit=24)
    close(unit=25)

    deallocate(Q)
    if (VR_opt.eq.1) then
        deallocate(Q_eq_VR)
    end if

    contains

    subroutine write_data_to_files_no_VR()
        implicit none
        call measure_no_VR

        121 format(F7.2, 2X, E15.8, 2X, E15.8)
        !$OMP single
        write(20,121) steps*dt, Qavg, Vqavg
        write(21,121) steps*dt, S, Serr
        write(22,121) steps*dt, Aeta, Veta
        write(23,121) steps*dt, Apsi, Vpsi
        write(24,121) steps*dt, Apsi2, Vpsi2
        !$OMP end single

    end subroutine

    subroutine write_data_to_files_with_VR()
        implicit none
        call measure_with_VR

        121 format(F7.2, 2X, E15.8, 2X, E15.8)
        !$OMP single
        write(20,121) steps*dt, Qavg, Vqavg
        write(21,121) steps*dt, S, Serr
        write(22,121) steps*dt, Aeta, Veta
        write(23,121) steps*dt, Apsi, Vpsi
        write(24,121) steps*dt, Apsi2, Vpsi2
        !$OMP end single

    end subroutine

    subroutine initialise_output_files_timestep()
        implicit none

        write(20,*) "Timestep width is: ", dt
        write(21,*) "Timestep width is: ", dt
        write(22,*) "Timestep width is: ", dt
        write(23,*) "Timestep width is: ", dt
        write(24,*) "Timestep width is: ", dt

    end subroutine

    subroutine measure_with_VR()
        implicit none
        real*8 :: Qavg_tmp, Vqavg_tmp, S_tmp, Serr_tmp
        real*8 :: Aeta_tmp, Veta_tmp, Apsi_tmp, Vpsi_tmp, Apsi2_tmp, Vpsi2_tmp

        tau = 0.D0
        ! These variables are all global and shared between threads
        Aeta = 0.D0; Apsi = 0.D0; Apsi2 = 0.D0
        Veta = 0.D0; Vpsi = 0.D0; Vpsi2 = 0.D0
        Qavg = 0.D0; Vqavg = 0.D0; S = 0.D0; Serr = 0.D0

        Aeta_tmp = 0.D0; Apsi_tmp = 0.D0; Apsi2_tmp = 0.D0
        Veta_tmp = 0.D0; Vpsi_tmp = 0.D0; Vpsi2_tmp = 0.D0
        Qavg_tmp = 0.D0; Vqavg_tmp = 0.D0; S_tmp = 0.D0; Serr_tmp = 0.D0

        !$OMP DO
        do i=1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product(delX, Q(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            B_eta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))

            Qavg_tmp = Qavg_tmp + Ql2
            Vqavg_tmp = Vqavg_tmp + Ql
            S_tmp = S_tmp + 0.5*(3*Bs - 1)
            Serr_tmp = Serr_tmp + 0.25*(9*Bs**2 - 6*Bs + 1)

            !subtract equilibrium values from shear-flow values
            Ql2 = Q_eq_VR(1,i)**2 + Q_eq_VR(2,i)**2 + Q_eq_VR(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product(delX, Q_eq_VR(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q_eq_VR(:,i)/Ql
            tau(:,:) = dyadic_prod(Q_eq_VR(:,i), F)

            S_tmp = S_tmp - 0.5*(3*Bs - 1)
            Serr_tmp = Serr_tmp - 0.25*(9*Bs**2 - 6*Bs + 1)

            B_eta = B_eta - tau(1,2)
            Bpsi = Bpsi - (tau(1,1) - tau(2,2))
            Bpsi2 = Bpsi2 - (tau(2,2) - tau(3,3))

            Aeta_tmp = Aeta_tmp + B_eta
            Apsi_tmp = Apsi_tmp + Bpsi
            Apsi2_tmp = Apsi2_tmp + Bpsi2
            Veta_tmp = Veta_tmp + B_eta**2
            Vpsi_tmp = Vpsi_tmp + Bpsi**2
            Vpsi2_tmp = Vpsi2_tmp + Bpsi2**2
        end do
        !$OMP END DO

        !$OMP atomic
        Aeta = Aeta + Aeta_tmp
        !$OMP atomic
        Veta = Veta + Veta_tmp
        !$OMP atomic
        Apsi = Apsi + Apsi_tmp
        !$OMP atomic
        Vpsi = Vpsi + Vpsi_tmp
        !$OMP atomic
        Apsi2 = Apsi2 + Apsi2_tmp
        !$OMP atomic
        Vpsi2 = Vpsi2 + Vpsi2_tmp
        !$OMP atomic
        Qavg = Qavg + Qavg_tmp
        !$OMP atomic
        Vqavg = Vqavg + Vqavg_tmp
        !$OMP atomic
        S = S + S_tmp
        !$OMP atomic
        Serr = Serr + Serr_tmp

        !$OMP barrier

        !$OMP single
        Aeta = Aeta/(Ntraj*sr)
        Veta = Veta/(Ntraj*sr**2)
        Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

        Apsi = Apsi/(Ntraj*sr**2)
        Vpsi = Vpsi/(Ntraj*sr**4)
        Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

        Apsi2 = Apsi2/(Ntraj*sr**2)
        Vpsi2 = Vpsi2/(Ntraj*sr**4)
        Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

        Qavg = sqrt(Qavg/Ntraj)
        Vqavg = Vqavg/Ntraj
        Vqavg = sqrt((Qavg**2 - Vqavg**2)/(Ntraj-1))

        S = S/Ntraj
        Serr = Serr/Ntraj
        Serr = sqrt((Serr - S**2)/(Ntraj-1))
        !$OMP end single

    end subroutine

    subroutine measure_no_VR()
        implicit none
        real*8 :: Qavg_tmp, Vqavg_tmp, S_tmp, Serr_tmp
        real*8 :: Aeta_tmp, Veta_tmp, Apsi_tmp, Vpsi_tmp, Apsi2_tmp, Vpsi2_tmp

        ! These variables are all global and shared between threads
        Aeta = 0.D0; Apsi = 0.D0; Apsi2 = 0.D0
        Veta = 0.D0; Vpsi = 0.D0; Vpsi2 = 0.D0
        Qavg = 0.D0; Vqavg = 0.D0; S = 0.D0; Serr = 0.D0

        !These variables SHOULD(!) be private
        tau = 0.D0
        Aeta_tmp = 0.D0; Apsi_tmp = 0.D0; Apsi2_tmp = 0.D0
        Veta_tmp = 0.D0; Vpsi_tmp = 0.D0; Vpsi2_tmp = 0.D0
        Qavg_tmp = 0.D0; Vqavg_tmp = 0.D0; S_tmp = 0.D0; Serr_tmp = 0.D0

        !$OMP DO
        do i=1,Ntraj
            Ql2 = Q(1,i)**2 + Q(2,i)**2 + Q(3,i)**2
            Ql = sqrt(Ql2)
            Bs = dot_product(delX, Q(:,i))**2/Ql2
            F(:) = (Ql - Q0)/(1.0D0-(Ql-Q0)**2/b)*Q(:,i)/Ql
            tau(:,:) = dyadic_prod(Q(:,i), F)

            Qavg_tmp = Qavg_tmp + Ql2
            Vqavg_tmp = Vqavg_tmp + Ql
            S_tmp = S_tmp + 0.5*(3*Bs - 1)
            Serr_tmp = Serr_tmp + 0.25*(9*Bs**2 - 6*Bs + 1)

            B_eta = tau(1,2)
            Bpsi = (tau(1,1) - tau(2,2))
            Bpsi2 = (tau(2,2) - tau(3,3))
            Aeta_tmp = Aeta_tmp + B_eta
            Apsi_tmp = Apsi_tmp + Bpsi
            Apsi2_tmp = Apsi2_tmp + Bpsi2
            Veta_tmp = Veta_tmp + B_eta**2
            Vpsi_tmp = Vpsi_tmp + Bpsi**2
            Vpsi2_tmp = Vpsi2_tmp + Bpsi2**2
        end do
        !$OMP END DO

        !$OMP atomic
        Aeta = Aeta + Aeta_tmp
        !$OMP atomic
        Veta = Veta + Veta_tmp
        !$OMP atomic
        Apsi = Apsi + Apsi_tmp
        !$OMP atomic
        Vpsi = Vpsi + Vpsi_tmp
        !$OMP atomic
        Apsi2 = Apsi2 + Apsi2_tmp
        !$OMP atomic
        Vpsi2 = Vpsi2 + Vpsi2_tmp
        !$OMP atomic
        Qavg = Qavg + Qavg_tmp
        !$OMP atomic
        Vqavg = Vqavg + Vqavg_tmp
        !$OMP atomic
        S = S + S_tmp
        !$OMP atomic
        Serr = Serr + Serr_tmp

        !$OMP barrier

        !$OMP single
        Aeta = Aeta/(Ntraj*sr)
        Veta = Veta/(Ntraj*sr**2)
        Veta = sqrt((Veta - Aeta**2)/(Ntraj-1))

        Apsi = Apsi/(Ntraj*sr**2)
        Vpsi = Vpsi/(Ntraj*sr**4)
        Vpsi = sqrt((Vpsi - Apsi**2)/(Ntraj-1))

        Apsi2 = Apsi2/(Ntraj*sr**2)
        Vpsi2 = Vpsi2/(Ntraj*sr**4)
        Vpsi2 = sqrt((Vpsi2 - Apsi2**2)/(Ntraj-1))

        Qavg = sqrt(Qavg/Ntraj)
        Vqavg = Vqavg/Ntraj
        Vqavg = sqrt((Qavg**2 - Vqavg**2)/(Ntraj-1))

        S = S/Ntraj
        Serr = Serr/Ntraj
        Serr = sqrt((Serr - S**2)/(Ntraj-1))
        !$OMP end single

    end subroutine

!    subroutine Write_data()
!        implicit none
!        open(unit=20, file='Ql.dat')
!        open(unit=21, file='S.dat')
!        open(unit=22, file='eta.dat')
!        open(unit=23, file='psi.dat')
!        open(unit=24, file='psi2.dat')
!
!        !10 format(F4.2,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5,4X,F7.5,2X,F7.5)
!        10 format(F6.3,4X, 2(E12.5, 2X, E12.5, 4X))
!        write(*,*) 'dtw       Q             err             S             err'
!        do i = 1,Ndtwidths
!            write(*,10) timestepwidths(i), Qavg(i), Vqavg(i), S(i), Serr(i)
!        end do
!
!        12 format(F6.3,4X, 3(E12.5, 2X, E12.5, 4X))
!        write(*,*) 'dtw       eta           err             psi           err       psi2        err'
!        do i = 1,Ndtwidths
!            write(*,12) timestepwidths(i), Aeta(i), Veta(i), Apsi(i), Vpsi(i), Apsi2(i), Vpsi2(i)
!        end do
!
!
!        11 format(F6.3,4X, E15.8, 2X, E15.8, 4X)
!        do i=20,24
!            write(i,*) 'timestepwidth    avg    err'
!        end do
!        do i = 1,Ndtwidths
!            write(20,11) timestepwidths(i), Qavg(i), Vqavg(i)
!            write(21,11) timestepwidths(i), S(i), Serr(i)
!            write(22,11) timestepwidths(i), Aeta(i), Veta(i)
!            write(23,11) timestepwidths(i), Apsi(i), Vpsi(i)
!            write(24,11) timestepwidths(i), Apsi2(i), Vpsi2(i)
!        end do
!
!        close(unit=20)
!        close(unit=21)
!        close(unit=22)
!        close(unit=23)
!        close(unit=24)
!
!    end subroutine

end program

