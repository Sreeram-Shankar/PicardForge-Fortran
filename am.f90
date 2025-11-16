module am_module
    implicit none
    
    abstract interface
        subroutine f_interface(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine f_interface
    end interface
    
contains

!solves the nonlinear system of equations with a Gauss-Seidel relaxation AM2
subroutine solve_am2(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: y(n), y_old(n), f_next(n), diff_norm
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
    
    !uses BE bootstrap
    call f(tgrid(2), Yout(1,:), f_next)
    Yout(2,:) = Yout(1,:) + h * f_next
    call f(tgrid(2), Yout(2,:), F_vals(2,:))
    
    do k = 2, nsteps
        y = Yout(k,:)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_next)
            y = Yout(k,:) + h * (0.5d0*f_next + 0.5d0*F_vals(k,:))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
        F_vals(k+1,:) = f_next
    end do
    
end subroutine solve_am2

!solves the nonlinear system of equations with a Gauss-Seidel relaxation AM3
subroutine solve_am3(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_boot(size(Yout,1), n), tgrid_boot(size(tgrid))
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: y(n), y_old(n), f_next(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !uses AM2 to bootstrap 1 step
    call solve_am2(f, t0, t0+2.0d0*h, y0, h, n, Y_boot, tgrid_boot, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1:3,:) = Y_boot(1:3,:)
    
    !compute initial F
    do i = 1, 3
        call f(tgrid(i), Yout(i,:), F_vals(i,:))
    end do
    
    !defines the main AM3 solver
    do k = 3, nsteps
        y = Yout(k,:)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_next)
            y = Yout(k,:) + h * ((5.0d0/12.0d0)*f_next + (2.0d0/3.0d0)*F_vals(k,:) - (1.0d0/12.0d0)*F_vals(k-1,:))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
        F_vals(k+1,:) = f_next
    end do
    
end subroutine solve_am3

!solves the nonlinear system of equations with a Gauss-Seidel relaxation AM4
subroutine solve_am4(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_boot(size(Yout,1), n), tgrid_boot(size(tgrid))
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: y(n), y_old(n), f_next(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !uses AM3 to bootstrap 2 steps
    call solve_am3(f, t0, t0+3.0d0*h, y0, h, n, Y_boot, tgrid_boot, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1:4,:) = Y_boot(1:4,:)
    
    !compute initial F
    do i = 1, 4
        call f(tgrid(i), Yout(i,:), F_vals(i,:))
    end do
    
    !defines the main AM4 solver
    do k = 4, nsteps
        y = Yout(k,:)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_next)
            y = Yout(k,:) + h * ((3.0d0/8.0d0)*f_next + (19.0d0/24.0d0)*F_vals(k,:) - &
                                 (5.0d0/24.0d0)*F_vals(k-1,:) + (1.0d0/24.0d0)*F_vals(k-2,:))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
        F_vals(k+1,:) = f_next
    end do
    
end subroutine solve_am4

!solves the nonlinear system of equations with a Gauss-Seidel relaxation AM5
subroutine solve_am5(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_boot(size(Yout,1), n), tgrid_boot(size(tgrid))
    real(8) :: F_vals(size(Yout,1), n)
    real(8) :: y(n), y_old(n), f_next(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !uses AM4 to bootstrap 3 steps of history
    call solve_am4(f, t0, t0+4.0d0*h, y0, h, n, Y_boot, tgrid_boot, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1:5,:) = Y_boot(1:5,:)
    
    !computes initial F
    do i = 1, 5
        call f(tgrid(i), Yout(i,:), F_vals(i,:))
    end do
    
    !defines the main AM5 solver
    do k = 5, nsteps
        y = Yout(k,:)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_next)
            y = Yout(k,:) + h * ((251.0d0/720.0d0)*f_next + (646.0d0/720.0d0)*F_vals(k,:) - &
                                 (264.0d0/720.0d0)*F_vals(k-1,:) + (106.0d0/720.0d0)*F_vals(k-2,:) - &
                                 (19.0d0/720.0d0)*F_vals(k-3,:))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
        F_vals(k+1,:) = f_next
    end do
    
end subroutine solve_am5

end module am_module