module irk_module
    implicit none
    
    !defines the Butcher tableau for Gauss-Legendre s = 1
    real(8), parameter :: A_gl_1(1,1) = reshape([0.5d0], [1,1])
    real(8), parameter :: b_gl_1(1) = [1.0d0]
    real(8), parameter :: c_gl_1(1) = [0.5d0]
    
    !defines the Butcher tableau for Gauss-Legendre s = 2
    real(8), parameter :: A_gl_2(2,2) = reshape([ &
        0.25d0, -0.03867513459481288d0, &
        0.5386751345948129d0, 0.25d0], [2,2])
    real(8), parameter :: b_gl_2(2) = [0.5d0, 0.5d0]
    real(8), parameter :: c_gl_2(2) = [0.7886751345948129d0, 0.21132486540518712d0]
    
    !defines the Butcher tableau for Gauss-Legendre s = 3
    real(8), parameter :: A_gl_3(3,3) = reshape([ &
        0.1388888888888889d0, -0.022485417203086815d0, 0.009789444015308326d0, &
        0.48042111196938335d0, 0.2222222222222222d0, -0.03597666752493890d0, &
        0.26798833376246945d0, 0.3002631949808646d0, 0.1388888888888889d0], [3,3])
    real(8), parameter :: b_gl_3(3) = [0.2777777777777778d0, 0.4444444444444444d0, 0.2777777777777778d0]
    real(8), parameter :: c_gl_3(3) = [0.8872983346207417d0, 0.5d0, 0.1127016653792583d0]
    
    !defines the Butcher tableau for Gauss-Legendre s = 4
    real(8), parameter :: A_gl_4(4,4) = reshape([ &
        0.08696371128436346d0, -0.014190694931141143d0, 0.006735500594538155d0, -0.003555149685795683d0, &
        0.35267675751627186d0, 0.16303628871563654d0, -0.027880428602470895d0, 0.012627462689404725d0, &
        0.31344511474186835d0, 0.3539530060337440d0, 0.16303628871563654d0, -0.026604180084998793d0, &
        0.1774825722545226d0, 0.16719192197418877d0, 0.18811811749986807d0, 0.08696371128436346d0], [4,4])
    real(8), parameter :: b_gl_4(4) = [0.17392742256872693d0, 0.32607257743127307d0, 0.32607257743127307d0, 0.17392742256872693d0]
    real(8), parameter :: c_gl_4(4) = [0.9305681557970263d0, 0.6699905217924281d0, 0.33000947820757187d0, 0.06943184420297371d0]
    
    !defines the Butcher tableau for Gauss-Legendre s = 5
    real(8), parameter :: A_gl_5(5,5) = reshape([ &
        0.05923172126404727d0, -0.00968756314195074d0, 0.004687154523869941d0, -0.002768994398769603d0, 0.001588112967865998d0, &
        0.25888469960875927d0, 0.11965716762484162d0, -0.020690316430958285d0, 0.010318280670683357d0, -0.005593793660812185d0, &
        0.2731900436258015d0, 0.30903655906408665d0, 0.14222222222222222d0, -0.024592114619642200d0, 0.011254400818642956d0, &
        0.24490812891049542d0, 0.22899605457899988d0, 0.2600046516806415d0, 0.11965716762484162d0, -0.019570364359076037d0, &
        0.11687532956022855d0, 0.12123243692686415d0, 0.11377628800422460d0, 0.12815100567004528d0, 0.05923172126404727d0], [5,5])
    real(8), parameter :: b_gl_5(5) = [0.11846344252809454d0, 0.23931433524968323d0, 0.28444444444444444d0, &
                                       0.23931433524968323d0, 0.11846344252809454d0]
    real(8), parameter :: c_gl_5(5) = [0.9530899229693320d0, 0.7692346550528416d0, 0.5d0, &
                                       0.23076534494715845d0, 0.04691007703066800d0]
    
    !defines the Butcher tableau for Radau IIA s = 2
    real(8), parameter :: A_radau_2(2,2) = reshape([ &
        0.41666666666666667d0, 0.75d0, &
        -0.08333333333333333d0, 0.25d0], [2,2])
    real(8), parameter :: b_radau_2(2) = [0.75d0, 0.25d0]
    real(8), parameter :: c_radau_2(2) = [0.3333333333333333d0, 1.0d0]
    
    !defines the Butcher tableau for Radau IIA s = 3
    real(8), parameter :: A_radau_3(3,3) = reshape([ &
        0.19681547722366043d0, 0.39442431473908727d0, 0.37640306270046727d0, &
        -0.06553542585019839d0, 0.29207341166522846d0, 0.5124858261884216d0, &
        0.023770974348220152d0, -0.04154875212599793d0, 0.11111111111111111d0], [3,3])
    real(8), parameter :: b_radau_3(3) = [0.37640306270046727d0, 0.5124858261884216d0, 0.11111111111111111d0]
    real(8), parameter :: c_radau_3(3) = [0.15505102572168219d0, 0.6449489742783178d0, 1.0d0]
    
    !defines the Butcher tableau for Radau IIA s = 4
    real(8), parameter :: A_radau_4(4,4) = reshape([ &
        0.20689257393535890d0, -0.04030922072352221d0, 0.4061232638673733d0, 0.38819346884317188d0, &
        0.23438399574740026d0, 0.11299947932315619d0, 0.21668178462325034d0, 0.22046221117676838d0, &
        -0.04785712804854072d0, 0.02580237742033639d0, 0.18903651817005634d0, 0.32884431998005974d0, &
        0.01604742280651627d0, -0.009904676507266424d0, -0.02418210489983294d0, 0.06250000000000000d0], [4,4])
    real(8), parameter :: b_radau_4(4) = [0.38819346884317188d0, 0.22046221117676838d0, 0.32884431998005974d0, 0.0625d0]
    real(8), parameter :: c_radau_4(4) = [0.4094668644407347d0, 0.08858795951270395d0, 0.7876594617608471d0, 1.0d0]
    
    !defines the Butcher tableau for Radau IIA s = 5
    real(8), parameter :: A_radau_5(5,5) = reshape([ &
        0.14621486784749350d0, -0.02673533110794557d0, 0.29896712949128348d0, 0.27650006876015923d0, 0.28135601514946206d0, &
        0.15377523147918247d0, 0.07299886431790332d0, 0.14006304568480987d0, 0.14489430810953476d0, 0.14371356079122594d0, &
        -0.03644456890512809d0, 0.01867692976398435d0, 0.16758507013524896d0, 0.32579792291042103d0, 0.31182652297574125d0, &
        0.02123306311930472d0, -0.01287910609330644d0, -0.03396910168661775d0, 0.12875675325490976d0, 0.22310390108357074d0, &
        -0.007935579902728778d0, 0.005042839233882015d0, 0.010944288744192252d0, -0.015708917378805328d0, 0.04000000000000000d0], [5,5])
    real(8), parameter :: b_radau_5(5) = [0.28135601514946206d0, 0.14371356079122594d0, 0.31182652297574125d0, &
                                          0.22310390108357074d0, 0.04d0]
    real(8), parameter :: c_radau_5(5) = [0.27684301363812383d0, 0.05710419611451768d0, &
                                          0.5835904323689168d0, 0.8602401356562194d0, 1.0d0]
    
    !defines the Butcher tableau for Lobatto IIIC s = 2
    real(8), parameter :: A_lob_2(2,2) = reshape([ &
        0.0d0, 0.5d0, &
        0.0d0, 0.5d0], [2,2])
    real(8), parameter :: b_lob_2(2) = [0.5d0, 0.5d0]
    real(8), parameter :: c_lob_2(2) = [0.0d0, 1.0d0]
    
    !defines the Butcher tableau for Lobatto IIIC s = 3
    real(8), parameter :: A_lob_3(3,3) = reshape([ &
        0.0d0, 0.20833333333333333d0, 0.16666666666666667d0, &
        0.0d0, 0.3333333333333333d0, 0.6666666666666666d0, &
        0.0d0, -0.041666666666666664d0, 0.16666666666666666d0], [3,3])
    real(8), parameter :: b_lob_3(3) = [0.16666666666666666d0, 0.6666666666666666d0, 0.16666666666666666d0]
    real(8), parameter :: c_lob_3(3) = [0.0d0, 0.5d0, 1.0d0]
    
    !defines the Butcher tableau for Lobatto IIIC s = 4
    real(8), parameter :: A_lob_4(4,4) = reshape([ &
        0.0d0, 0.07303276685416842d0, 0.11030056647916491d0, 0.08333333333333333d0, &
        0.0d0, 0.22696723314583158d0, -0.03390736422914388d0, 0.4166666666666667d0, &
        0.0d0, 0.45057403089581055d0, 0.18969943352083508d0, 0.4166666666666667d0, &
        0.0d0, -0.02696723314583158d0, 0.01030056647916491d0, 0.08333333333333333d0], [4,4])
    real(8), parameter :: b_lob_4(4) = [0.08333333333333333d0, 0.4166666666666667d0, 0.4166666666666667d0, 0.08333333333333333d0]
    real(8), parameter :: c_lob_4(4) = [0.0d0, 0.7236067977499790d0, 0.2763932022500210d0, 1.0d0]
    
    !defines the Butcher tableau for Lobatto IIIC s = 5
    real(8), parameter :: A_lob_5(5,5) = reshape([ &
        0.0d0, 0.05370013924241453d0, 0.040625d0, 0.06772843218615690d0, 0.05d0, &
        0.0d0, 0.15247745287881054d0, -0.030961961100820556d0, 0.01063582422541549d0, 0.2722222222222222d0, &
        0.0d0, 0.37729127742211367d0, 0.17777777777777778d0, -0.021735721866558114d0, 0.35555555555555557d0, &
        0.0d0, 0.26158639799680673d0, 0.30318418332304278d0, 0.11974476934341168d0, 0.2722222222222222d0, &
        0.0d0, -0.017728432186156897d0, 0.009375d0, -0.003700139242414531d0, 0.05d0], [5,5])
    real(8), parameter :: b_lob_5(5) = [0.05d0, 0.2722222222222222d0, 0.35555555555555557d0, &
                                        0.2722222222222222d0, 0.05d0]
    real(8), parameter :: c_lob_5(5) = [0.0d0, 0.8273268353539886d0, 0.5d0, &
                                        0.1726731646460114d0, 1.0d0]
    
contains

!loads the Butcher tableaus for the given family and stage count
subroutine get_tableau(family, s, A, b, c)
    implicit none
    character(len=*), intent(in) :: family
    integer, intent(in) :: s
    real(8), intent(out) :: A(s,s), b(s), c(s)
    
    character(len=20) :: family_lower
    integer :: i
    
    !convert to lowercase
    family_lower = family
    do i = 1, len_trim(family_lower)
        if (family_lower(i:i) >= 'A' .and. family_lower(i:i) <= 'Z') then
            family_lower(i:i) = char(ichar(family_lower(i:i)) + 32)
        end if
    end do
    
    if (family_lower(1:5) == 'gauss') then
        if (s == 1) then
            A = A_gl_1
            b = b_gl_1
            c = c_gl_1
        else if (s == 2) then
            A = A_gl_2
            b = b_gl_2
            c = c_gl_2
        else if (s == 3) then
            A = A_gl_3
            b = b_gl_3
            c = c_gl_3
        else if (s == 4) then
            A = A_gl_4
            b = b_gl_4
            c = c_gl_4
        else if (s == 5) then
            A = A_gl_5
            b = b_gl_5
            c = c_gl_5
        else
            error stop 'Gauss-Legendre only implemented for s = 1...5'
        end if
        
    else if (family_lower(1:5) == 'radau') then
        if (s == 2) then
            A = A_radau_2
            b = b_radau_2
            c = c_radau_2
        else if (s == 3) then
            A = A_radau_3
            b = b_radau_3
            c = c_radau_3
        else if (s == 4) then
            A = A_radau_4
            b = b_radau_4
            c = c_radau_4
        else if (s == 5) then
            A = A_radau_5
            b = b_radau_5
            c = c_radau_5
        else
            error stop 'Radau IIA only implemented for s = 2...5'
        end if
        
    else if (family_lower(1:7) == 'lobatto') then
        if (s == 2) then
            A = A_lob_2
            b = b_lob_2
            c = c_lob_2
        else if (s == 3) then
            A = A_lob_3
            b = b_lob_3
            c = c_lob_3
        else if (s == 4) then
            A = A_lob_4
            b = b_lob_4
            c = c_lob_4
        else if (s == 5) then
            A = A_lob_5
            b = b_lob_5
            c = c_lob_5
        else
            error stop 'Lobatto IIIC only implemented for s = 2...5'
        end if
        
    else
        error stop 'Unknown IRK family. Must be "gauss", "radau", or "lobatto".'
    end if
    
end subroutine get_tableau


!defines the IRK step with Gauss-Seidel relaxation
subroutine step_collocation(f, t, y, h, A, b, c, s, n, y_next, sweeps, tol)
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
        
        !computes L2 norm of all stage differences (equivalent to Python's np.linalg.norm on 2D array)
        diff_norm = sqrt(sum((Y_stages - Y_old)**2))
        if (diff_norm < tol) exit
    end do
    
    !computes the final state update
    y_next = y
    do i = 1, s
        call f(t + c(i)*h, Y_stages(i,:), fval)
        y_next = y_next + h*b(i)*fval
    end do
    
end subroutine step_collocation


!main solver for any collocation IRK method using Gauss-Seidel relaxation
subroutine solve_collocation(f, t0, tf, y0, h, n, Yout, tgrid, family, s, sweeps, tol)
    implicit none
    integer, intent(in) :: n, s, sweeps
    real(8), intent(in) :: t0, tf, h, tol
    real(8), intent(in) :: y0(n)
    character(len=*), intent(in) :: family
    real(8), intent(out) :: Yout(:,:), tgrid(:)
    
    real(8) :: A(s,s), b(s), c(s)
    real(8) :: y(n)
    integer :: nsteps, k
    
    interface
        subroutine f(t_arg, y_arg, dydt)
            real(8), intent(in) :: t_arg
            real(8), intent(in) :: y_arg(:)
            real(8), intent(out) :: dydt(:)
        end subroutine
    end interface
    
    !gets the Butcher tableau
    call get_tableau(family, s, A, b, c)
    
    nsteps = ceiling((tf - t0)/h)
    
    do k = 1, nsteps+1
        tgrid(k) = t0 + (k-1)*h
    end do
    
    Yout(1,:) = y0
    y = y0
    
    do k = 1, nsteps
        call step_collocation(f, tgrid(k), y, h, A, b, c, s, n, Yout(k+1,:), sweeps, tol)
        y = Yout(k+1,:)
    end do
    
end subroutine solve_collocation

end module irk_module