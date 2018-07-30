
program General_Dumbbell
    use Integrate
    use Generate_Initial_Conditions
    use Dumbbell_util
    use omp_lib
    implicit none

    integer*8 :: Nblocks, Nsteps, steps, block, time(1:8), seed, i, Ntraj, VR_opt, dist_opt
    integer :: NTimeSteps, timestep, trajectories, Ndtwidths
    real*8 :: sr, alpha, h, a, Q0, Nrelax_times, dt, dW(3), Bq, Bs, time1, time2, output_delay, k(3,3)
    type(measured_variables) :: out_var
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
    read (31, *) sr, alpha, h, Q0
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

    looptimesteps: do timestep=1,Ndtwidths

        dt = timestepwidths(timestep)
        NtimeSteps = int(Nrelax_times/dt)

        !Also works for Fraenkel, FENE and Hookean springs
        !Reduce final input for faster computation of equilibrium dist
        Q(:,:) = generate_Q_FF(Q0, alpha, Ntraj, seed, 10000)

        call initialise_output_files_timestep()

        !$OMP PARALLEL DEFAULT(firstprivate) SHARED(Q, Q_eq_VR, VR_opt, out_var)
        !$ seed = seed + 932117 + OMP_get_thread_num()*2685821657736338717_8

        if (VR_opt.eq.0) then
            do steps=1,Ntimesteps
                k(1,2) = sr
                !$OMP DO schedule(dynamic, 100)
                do i=1,Ntraj
                    dW = Wiener_step(seed, dt)
                    Q(:,i) =  step(Q(:,i), k, dt, Q0, alpha, a, dW)
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
                    Q(:,i) =  step(Q(:,i), k, dt, Q0, alpha, a, dW)
                    k(1,2) = 0.D0
                    Q_eq_VR(:,i) = step(Q_eq_VR(:,i), k, dt, Q0, alpha, a, dW)
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
        call measure_shear_no_VR(out_var, Q, Q0, alpha, sr, Ntraj)

        121 format(F7.2, 2X, E15.8, 2X, E15.8)
        !$OMP single
        write(20,121) steps*dt, out_var%Qavg, out_var%Vqavg
        write(21,121) steps*dt, out_var%S, out_var%Serr
        write(22,121) steps*dt, out_var%Aeta, out_var%Veta
        write(23,121) steps*dt, out_var%Apsi, out_var%Vpsi
        write(24,121) steps*dt, out_var%Apsi2, out_var%Vpsi2
        !$OMP end single

    end subroutine

    subroutine write_data_to_files_with_VR()
        implicit none
        call measure_shear_with_VR(out_var, Q, Q_eq_VR, Q0, alpha, sr, Ntraj)

        121 format(F7.2, 2X, E15.8, 2X, E15.8)
        !$OMP single
        write(20,121) steps*dt, out_var%Qavg, out_var%Vqavg
        write(21,121) steps*dt, out_var%S, out_var%Serr
        write(22,121) steps*dt, out_var%Aeta, out_var%Veta
        write(23,121) steps*dt, out_var%Apsi, out_var%Vpsi
        write(24,121) steps*dt, out_var%Apsi2, out_var%Vpsi2
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

end program

