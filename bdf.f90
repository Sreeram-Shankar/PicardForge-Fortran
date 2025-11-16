module bdf_module
    implicit none
contains

!solves the nonlinear system of equations with a Gauss-Seidel relaxation BE (BDF1)
subroutine solve_be(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: y(n), y_old(n), f_val(n), diff_norm
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
    
    !defines the main BE solver
    do k = 1, nsteps
        y = Yout(k,:)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_val)
            y = Yout(k,:) + h * f_val
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
    end do
    
end subroutine solve_be

!solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF2
subroutine solve_bdf2(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_BE(size(Yout,1), n), tgrid_BE(size(tgrid))
    real(8) :: y(n), y_old(n), f_val(n), rhs(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !uses backward euler to bootstrap
    call solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    Yout(2,:) = Y_BE(2,:)
    
    do k = 2, nsteps
        y = Yout(k,:)
        
        !defines the rhs
        rhs = (-4.0d0*Yout(k,:) + Yout(k-1,:)) / (2.0d0*h)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_val)
            y = (rhs + f_val) / (3.0d0/(2.0d0*h))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
    end do
    
end subroutine solve_bdf2

!solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF3
subroutine solve_bdf3(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_BE(size(Yout,1), n), tgrid_BE(size(tgrid))
    real(8) :: Y_BDF2(size(Yout,1), n), tgrid_BDF2(size(tgrid))
    real(8) :: y(n), y_old(n), f_val(n), rhs(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !bootstrap using BE then BDF2
    call solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol)
    call solve_bdf2(f, t0, t0+2.0d0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    Yout(2,:) = Y_BE(2,:)
    Yout(3,:) = Y_BDF2(3,:)
    
    do k = 3, nsteps
        y = Yout(k,:)
        
        !defines the rhs
        rhs = (-11.0d0*Yout(k,:) + 18.0d0*Yout(k-1,:) - 9.0d0*Yout(k-2,:) + 2.0d0*Yout(k-3,:)) / (6.0d0*h)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_val)
            y = (rhs + f_val) / (11.0d0/(6.0d0*h))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
    end do
    
end subroutine solve_bdf3

!solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF4
subroutine solve_bdf4(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_BE(size(Yout,1), n), tgrid_BE(size(tgrid))
    real(8) :: Y_BDF2(size(Yout,1), n), tgrid_BDF2(size(tgrid))
    real(8) :: Y_BDF3(size(Yout,1), n), tgrid_BDF3(size(tgrid))
    real(8) :: y(n), y_old(n), f_val(n), rhs(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !bootstrap using BDF1–3
    call solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol)
    call solve_bdf2(f, t0, t0+2.0d0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol)
    call solve_bdf3(f, t0, t0+3.0d0*h, y0, h, n, Y_BDF3, tgrid_BDF3, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    Yout(2,:) = Y_BE(2,:)
    Yout(3,:) = Y_BDF2(3,:)
    Yout(4,:) = Y_BDF3(4,:)
    
    do k = 4, nsteps
        y = Yout(k,:)
        
        !defines the rhs
        rhs = (-25.0d0*Yout(k,:) + 48.0d0*Yout(k-1,:) - 36.0d0*Yout(k-2,:) + &
              16.0d0*Yout(k-3,:) - 3.0d0*Yout(k-4,:)) / (12.0d0*h)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_val)
            y = (rhs + f_val) / (25.0d0/(12.0d0*h))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
    end do
    
end subroutine solve_bdf4

!solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF5
subroutine solve_bdf5(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_BE(size(Yout,1), n), tgrid_BE(size(tgrid))
    real(8) :: Y_BDF2(size(Yout,1), n), tgrid_BDF2(size(tgrid))
    real(8) :: Y_BDF3(size(Yout,1), n), tgrid_BDF3(size(tgrid))
    real(8) :: Y_BDF4(size(Yout,1), n), tgrid_BDF4(size(tgrid))
    real(8) :: y(n), y_old(n), f_val(n), rhs(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !bootstrap up to BDF4
    call solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol)
    call solve_bdf2(f, t0, t0+2.0d0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol)
    call solve_bdf3(f, t0, t0+3.0d0*h, y0, h, n, Y_BDF3, tgrid_BDF3, sweeps, tol)
    call solve_bdf4(f, t0, t0+4.0d0*h, y0, h, n, Y_BDF4, tgrid_BDF4, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    Yout(2,:) = Y_BE(2,:)
    Yout(3,:) = Y_BDF2(3,:)
    Yout(4,:) = Y_BDF3(4,:)
    Yout(5,:) = Y_BDF4(5,:)
    
    do k = 5, nsteps
        y = Yout(k,:)
        
        !defines the rhs
        rhs = (-137.0d0*Yout(k,:) + 300.0d0*Yout(k-1,:) - 300.0d0*Yout(k-2,:) + &
              200.0d0*Yout(k-3,:) - 75.0d0*Yout(k-4,:) + 12.0d0*Yout(k-5,:)) / (60.0d0*h)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_val)
            y = (rhs + f_val) / (137.0d0/(60.0d0*h))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
    end do
    
end subroutine solve_bdf5

!solves the nonlinear system of equations with a Gauss-Seidel relaxation BDF6
subroutine solve_bdf6(f, t0, tf, y0, h, n, Yout, tgrid, sweeps, tol)
    implicit none
    integer, intent(in) :: n, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: Y_BE(size(Yout,1), n), tgrid_BE(size(tgrid))
    real(8) :: Y_BDF2(size(Yout,1), n), tgrid_BDF2(size(tgrid))
    real(8) :: Y_BDF3(size(Yout,1), n), tgrid_BDF3(size(tgrid))
    real(8) :: Y_BDF4(size(Yout,1), n), tgrid_BDF4(size(tgrid))
    real(8) :: Y_BDF5(size(Yout,1), n), tgrid_BDF5(size(tgrid))
    real(8) :: y(n), y_old(n), f_val(n), rhs(n), diff_norm
    integer :: nsteps, k, i
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !bootstrap up to BDF5
    call solve_be(f, t0, t0+h, y0, h, n, Y_BE, tgrid_BE, sweeps, tol)
    call solve_bdf2(f, t0, t0+2.0d0*h, y0, h, n, Y_BDF2, tgrid_BDF2, sweeps, tol)
    call solve_bdf3(f, t0, t0+3.0d0*h, y0, h, n, Y_BDF3, tgrid_BDF3, sweeps, tol)
    call solve_bdf4(f, t0, t0+4.0d0*h, y0, h, n, Y_BDF4, tgrid_BDF4, sweeps, tol)
    call solve_bdf5(f, t0, t0+5.0d0*h, y0, h, n, Y_BDF5, tgrid_BDF5, sweeps, tol)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    Yout(2,:) = Y_BE(2,:)
    Yout(3,:) = Y_BDF2(3,:)
    Yout(4,:) = Y_BDF3(4,:)
    Yout(5,:) = Y_BDF4(5,:)
    Yout(6,:) = Y_BDF5(6,:)
    
    do k = 6, nsteps
        y = Yout(k,:)
        
        !defines the rhs
        rhs = (-147.0d0*Yout(k,:) + 360.0d0*Yout(k-1,:) - 450.0d0*Yout(k-2,:) + &
              400.0d0*Yout(k-3,:) - 225.0d0*Yout(k-4,:) + 72.0d0*Yout(k-5,:) - &
              10.0d0*Yout(k-6,:)) / (60.0d0*h)
        
        !implements Gauss-Seidel relaxation
        do i = 1, sweeps
            y_old = y
            call f(tgrid(k+1), y, f_val)
            y = (rhs + f_val) / (147.0d0/(60.0d0*h))
            diff_norm = sqrt(sum((y - y_old)**2))
            if (diff_norm < tol) exit
        end do
        
        Yout(k+1,:) = y
    end do
    
end subroutine solve_bdf6

end module bdf_module