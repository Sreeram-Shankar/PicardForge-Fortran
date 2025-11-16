module rk_module
    implicit none
contains

!defines the step for RK1
subroutine step_rk1(f, t, y, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, h
    real(8), intent(in) :: y(n)
    real(8), intent(out) :: y_next(n)
    
    real(8) :: k1(n)
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    call f(t, y, k1)
    y_next = y + h * k1
end subroutine step_rk1

!main solver for RK1
subroutine solve_rk1(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    integer :: nsteps, k
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    
    do k = 1, nsteps
        call step_rk1(f, tgrid(k), Yout(k,:), h, n, Yout(k+1,:))
    end do
    
end subroutine solve_rk1

!defines the step for RK2
subroutine step_rk2(f, t, y, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, h
    real(8), intent(in) :: y(n)
    real(8), intent(out) :: y_next(n)
    
    real(8) :: k1(n), k2(n)
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    call f(t, y, k1)
    call f(t + h, y + h*k1, k2)
    y_next = y + (h / 2.0d0) * (k1 + k2)
end subroutine step_rk2

!main solver for RK2
subroutine solve_rk2(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    integer :: nsteps, k
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    
    do k = 1, nsteps
        call step_rk2(f, tgrid(k), Yout(k,:), h, n, Yout(k+1,:))
    end do
    
end subroutine solve_rk2

!defines the step for RK3
subroutine step_rk3(f, t, y, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, h
    real(8), intent(in) :: y(n)
    real(8), intent(out) :: y_next(n)
    
    real(8) :: k1(n), k2(n), k3(n)
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    call f(t, y, k1)
    call f(t + h/2.0d0, y + (h/2.0d0)*k1, k2)
    call f(t + h, y + h*(-k1 + 2.0d0*k2), k3)
    y_next = y + (h / 6.0d0) * (k1 + 4.0d0*k2 + k3)
end subroutine step_rk3

!main solver for RK3
subroutine solve_rk3(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    integer :: nsteps, k
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    
    do k = 1, nsteps
        call step_rk3(f, tgrid(k), Yout(k,:), h, n, Yout(k+1,:))
    end do
    
end subroutine solve_rk3

!defines the step for RK4
subroutine step_rk4(f, t, y, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, h
    real(8), intent(in) :: y(n)
    real(8), intent(out) :: y_next(n)
    
    real(8) :: k1(n), k2(n), k3(n), k4(n)
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    call f(t, y, k1)
    call f(t + h/2.0d0, y + (h/2.0d0)*k1, k2)
    call f(t + h/2.0d0, y + (h/2.0d0)*k2, k3)
    call f(t + h, y + h*k3, k4)
    y_next = y + (h / 6.0d0) * (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
end subroutine step_rk4

!main solver for RK4
subroutine solve_rk4(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    integer :: nsteps, k
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    
    do k = 1, nsteps
        call step_rk4(f, tgrid(k), Yout(k,:), h, n, Yout(k+1,:))
    end do
    
end subroutine solve_rk4

!defines the step for RK5
subroutine step_rk5(f, t, y, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t, h
    real(8), intent(in) :: y(n)
    real(8), intent(out) :: y_next(n)
    
    real(8) :: k1(n), k2(n), k3(n), k4(n), k5(n), k6(n)
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    call f(t, y, k1)
    call f(t + h/4.0d0, y + (h/4.0d0)*k1, k2)
    call f(t + h/4.0d0, y + (h/8.0d0)*(k1 + k2), k3)
    call f(t + h/2.0d0, y + (h/2.0d0)*k3, k4)
    call f(t + 3.0d0*h/4.0d0, y + (h/16.0d0)*(3.0d0*k1 + 9.0d0*k4), k5)
    call f(t + h, y + (h/7.0d0)*(2.0d0*k1 + 3.0d0*k2 + 4.0d0*k4 - 12.0d0*k3), k6)
    y_next = y + (h / 90.0d0) * (7.0d0*k1 + 32.0d0*k3 + 12.0d0*k4 + 32.0d0*k5 + 7.0d0*k6)
end subroutine step_rk5

!main solver for RK5
subroutine solve_rk5(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    integer :: nsteps, k
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    
    do k = 1, nsteps
        call step_rk5(f, tgrid(k), Yout(k,:), h, n, Yout(k+1,:))
    end do
    
end subroutine solve_rk5

end module rk_module