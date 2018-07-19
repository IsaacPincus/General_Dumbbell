module Random_Numbers

    implicit none

    contains

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

end module
