module ab_module
    implicit none
    
    abstract interface
        subroutine f_interface(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine f_interface
    end interface
    
contains

!defines the step for AB2
subroutine step_ab2(y_n, f_n, f_prev, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: h
    real(8), intent(in) :: y_n(n), f_n(n), f_prev(n)
    real(8), intent(out) :: y_next(n)
    
    y_next = y_n + h * ((3.0d0/2.0d0)*f_n - (1.0d0/2.0d0)*f_prev)
end subroutine step_ab2

!main solver for AB2
subroutine solve_ab2(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: k1(n), k2(n)
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
    call f(tgrid(1), Yout(1,:), F_vals(1,:))
    
    !uses RK2 for first step
    call f(tgrid(1), Yout(1,:), k1)
    call f(tgrid(1) + h, Yout(1,:) + h*k1, k2)
    Yout(2,:) = Yout(1,:) + 0.5d0 * h * (k1 + k2)
    call f(tgrid(2), Yout(2,:), F_vals(2,:))
    
    do k = 2, nsteps
        call step_ab2(Yout(k,:), F_vals(k,:), F_vals(k-1,:), h, n, Yout(k+1,:))
        call f(tgrid(k+1), Yout(k+1,:), F_vals(k+1,:))
    end do
    
end subroutine solve_ab2


!defines the step for AB3
subroutine step_ab3(y_n, f_n, f_prev, f_prev2, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: h
    real(8), intent(in) :: y_n(n), f_n(n), f_prev(n), f_prev2(n)
    real(8), intent(out) :: y_next(n)
    
    y_next = y_n + h * ((23.0d0/12.0d0)*f_n - (16.0d0/12.0d0)*f_prev + (5.0d0/12.0d0)*f_prev2)
end subroutine step_ab3

!main solver for AB3
subroutine solve_ab3(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: k1(n), k2(n), k3(n)
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
    call f(tgrid(1), Yout(1,:), F_vals(1,:))
    
    !uses RK3 for first two steps
    call f(tgrid(1), Yout(1,:), k1)
    call f(tgrid(1) + h/2.0d0, Yout(1,:) + (h/2.0d0)*k1, k2)
    call f(tgrid(1) + h, Yout(1,:) + h*(-k1 + 2.0d0*k2), k3)
    Yout(2,:) = Yout(1,:) + (h/6.0d0) * (k1 + 4.0d0*k2 + k3)
    call f(tgrid(2), Yout(2,:), F_vals(2,:))
    
    call f(tgrid(2), Yout(2,:), k1)
    call f(tgrid(2) + h/2.0d0, Yout(2,:) + (h/2.0d0)*k1, k2)
    call f(tgrid(2) + h, Yout(2,:) + h*(-k1 + 2.0d0*k2), k3)
    Yout(3,:) = Yout(2,:) + (h/6.0d0) * (k1 + 4.0d0*k2 + k3)
    call f(tgrid(3), Yout(3,:), F_vals(3,:))
    
    do k = 3, nsteps
        call step_ab3(Yout(k,:), F_vals(k,:), F_vals(k-1,:), F_vals(k-2,:), h, n, Yout(k+1,:))
        call f(tgrid(k+1), Yout(k+1,:), F_vals(k+1,:))
    end do
    
end subroutine solve_ab3


!defines the step for AB4
subroutine step_ab4(y_n, f_n, f_prev, f_prev2, f_prev3, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: h
    real(8), intent(in) :: y_n(n)
    real(8), intent(in) :: f_n(n), f_prev(n), f_prev2(n), f_prev3(n)
    real(8), intent(out) :: y_next(n)
    
    y_next = y_n + h * ((55.0d0/24.0d0)*f_n - (59.0d0/24.0d0)*f_prev + &
                        (37.0d0/24.0d0)*f_prev2 - (9.0d0/24.0d0)*f_prev3)
end subroutine step_ab4

!main solver for AB4
subroutine solve_ab4(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: k1(n), k2(n), k3(n), k4(n)
    integer :: nsteps, k, i
    
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
    call f(tgrid(1), Yout(1,:), F_vals(1,:))
    
    !uses RK4 for first three steps
    do i = 1, 3
        call f(tgrid(i), Yout(i,:), k1)
        call f(tgrid(i) + h/2.0d0, Yout(i,:) + (h/2.0d0)*k1, k2)
        call f(tgrid(i) + h/2.0d0, Yout(i,:) + (h/2.0d0)*k2, k3)
        call f(tgrid(i) + h, Yout(i,:) + h*k3, k4)
        Yout(i+1,:) = Yout(i,:) + (h/6.0d0) * (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
        call f(tgrid(i+1), Yout(i+1,:), F_vals(i+1,:))
    end do
    
    do k = 4, nsteps
        call step_ab4(Yout(k,:), F_vals(k,:), F_vals(k-1,:), F_vals(k-2,:), F_vals(k-3,:), h, n, Yout(k+1,:))
        call f(tgrid(k+1), Yout(k+1,:), F_vals(k+1,:))
    end do
    
end subroutine solve_ab4


!defines the step for AB5
subroutine step_ab5(y_n, f_n, f_prev, f_prev2, f_prev3, f_prev4, h, n, y_next)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: h
    real(8), intent(in) :: y_n(n)
    real(8), intent(in) :: f_n(n), f_prev(n), f_prev2(n), f_prev3(n), f_prev4(n)
    real(8), intent(out) :: y_next(n)
    
    y_next = y_n + h * ((1901.0d0/720.0d0)*f_n - (2774.0d0/720.0d0)*f_prev + &
                        (2616.0d0/720.0d0)*f_prev2 - (1274.0d0/720.0d0)*f_prev3 + &
                        (251.0d0/720.0d0)*f_prev4)
end subroutine step_ab5

!main solver for AB5
subroutine solve_ab5(f, t0, tf, y0, h, n, Yout, tgrid)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: t0, tf, h
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: k1(n), k2(n), k3(n), k4(n)
    integer :: nsteps, k, i
    
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
    call f(tgrid(1), Yout(1,:), F_vals(1,:))
    
    !uses RK4 for first four steps
    do i = 1, 4
        call f(tgrid(i), Yout(i,:), k1)
        call f(tgrid(i) + h/2.0d0, Yout(i,:) + (h/2.0d0)*k1, k2)
        call f(tgrid(i) + h/2.0d0, Yout(i,:) + (h/2.0d0)*k2, k3)
        call f(tgrid(i) + h, Yout(i,:) + h*k3, k4)
        Yout(i+1,:) = Yout(i,:) + (h/6.0d0) * (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
        call f(tgrid(i+1), Yout(i+1,:), F_vals(i+1,:))
    end do
    
    do k = 5, nsteps
        call step_ab5(Yout(k,:), F_vals(k,:), F_vals(k-1,:), F_vals(k-2,:), F_vals(k-3,:), F_vals(k-4,:), h, n, Yout(k+1,:))
        call f(tgrid(k+1), Yout(k+1,:), F_vals(k+1,:))
    end do
    
end subroutine solve_ab5

end module ab_module