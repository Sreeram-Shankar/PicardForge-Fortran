module sdirk_module
    implicit none
contains

!defines the SDIRK step with Gauss-Seidel relaxation
subroutine step_sdirk(f, t, y, h, A, b, c, s, n, y_next, sweeps, tol)
    implicit none
    integer, intent(in) :: s, n, sweeps
    real(8), intent(in) :: t, h, tol
    real(8), intent(in) :: y(n), A(s,s), b(s), c(s)
    real(8), intent(out) :: y_next(n)

    real(8) :: Y_stages(s,n), Y_old(s,n), rhs(n), fval(n)
    real(8) :: diff_norm
    integer :: i, j, k

    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface

    !initializes all stages with the initial guess y
    do i = 1, s
        Y_stages(i,:) = y(:)
    end do

    !implements Gauss-Seidel relaxation
    do k = 1, sweeps
        Y_old = Y_stages
        do i = 1, s
            rhs = 0.0d0
            do j = 1, s
                call f(t + c(j)*h, Y_stages(j,:), fval)
                rhs = rhs + A(i,j)*fval
            end do
            Y_stages(i,:) = y + h*rhs
        end do

        diff_norm = sqrt(sum((Y_stages - Y_old)**2))
        if (diff_norm < tol) exit
    end do

    !computes the final state update
    y_next = y
    do i = 1, s
        call f(t + c(i)*h, Y_stages(i,:), fval)
        y_next = y_next + h*b(i)*fval
    end do

end subroutine step_sdirk


!solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK2
subroutine solve_sdirk2(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)

    real(8) :: gamma, A(2,2), b(2), c(2)
    real(8) :: y(n)
    integer :: nsteps, k

    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface

    gamma = 1.0d0 - 1.0d0/sqrt(2.0d0)
    A = reshape([gamma,0.0d0, 1.0d0-gamma, gamma],[2,2])
    b = [1.0d0-gamma, gamma]
    c = [gamma, 1.0d0]

    nsteps = ceiling((tf - t0)/h)

    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do

    Yout(1,:) = y0
    y = y0

    do k = 1, nsteps
        call step_sdirk(f, tgrid(k), y, h, A, b, c, 2, n, Yout(k+1,:), sweeps, tol)
        y = Yout(k+1,:)
    end do

end subroutine solve_sdirk2


!solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK3
subroutine solve_sdirk3(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out):: Yout(:,:), tgrid(:)

    real(8) :: gamma, A(3,3), b(3), c(3)
    real(8) :: y(n)
    integer :: nsteps, k

    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out):: dydt(:)
        end subroutine
    end interface

    gamma = 0.435866521508459d0
    A = reshape([ &
        gamma, 0.0d0, 0.0d0, &
        0.2820667395d0, gamma, 0.0d0, &
        1.208496649d0, -0.644363171d0, gamma], [3,3])

    b = [1.208496649d0, -0.644363171d0, gamma]
    c = [gamma, 0.7179332605d0, 1.0d0]

    nsteps = ceiling((tf - t0)/h)

    do k=1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do

    Yout(1,:) = y0
    y = y0

    do k=1, nsteps
        call step_sdirk(f, tgrid(k), y, h, A, b, c, 3, n, Yout(k+1,:), sweeps, tol)
        y = Yout(k+1,:)
    end do

end subroutine solve_sdirk3


!solves the nonlinear system of equations with a Gauss-Seidel relaxation SDIRK4
subroutine solve_sdirk4(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out):: Yout(:,:), tgrid(:)

    real(8) :: gamma, A(4,4), b(4), c(4)
    real(8) :: y(n)
    integer :: nsteps, k

    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out):: dydt(:)
        end subroutine
    end interface

    gamma = 0.572816062482135d0
    A = reshape([ &
        gamma, 0.0d0, 0.0d0, 0.0d0, &
       -0.6557110092d0, gamma, 0.0d0, 0.0d0, &
        0.757184241d0, 0.237758128d0, gamma, 0.0d0, &
        0.155416858d0, 0.701913790d0, 0.142669351d0, gamma ], [4,4])

    b = [0.155416858d0, 0.701913790d0, 0.142669351d0, gamma]
    c = [gamma, 0.344d0, 0.995d0, 1.0d0]

    nsteps = ceiling((tf - t0)/h)

    do k=1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do

    Yout(1,:) = y0
    y = y0

    do k=1, nsteps
        call step_sdirk(f, tgrid(k), y, h, A, b, c, 4, n, Yout(k+1,:), sweeps, tol)
        y = Yout(k+1,:)
    end do

end subroutine solve_sdirk4

end module sdirk_module