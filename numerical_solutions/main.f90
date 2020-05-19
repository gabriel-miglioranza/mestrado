program main
    use types
    use simulation_parameters
    use simulation
    implicit none
    integer :: neq
    real(dp), allocatable, dimension(:) :: y
    real(dp) :: t0, tf
    real(dp), allocatable, dimension(:) :: atol, rtol
    integer, dimension(10) :: ipar
    real(dp), dimension(10) :: rpar 
    external jac, psol
    model_id = 4
    bi = 10.0_dp
    alfa = 1.0_dp
    cc1 = 1
    cc2 = 0
    s = 2 
    n = 300   
    dx = 1.0_dp/real(n + 1, dp)
    neq = n + 1 + cc1 + cc2
    
    allocate(y(neq), atol(neq), rtol(neq))

    t0 = 0.0_dp
    tf = 1.0_dp
    rtol = 1.0e-6_dp
    atol = 1.0e-6_dp
    rpar = 0.0_dp
    ipar = 0

    y = 0.0_dp
    y(neq) = 1.0_dp
    
    call ddaspk_integration(model, neq, y, t0, tf, rtol, atol, rpar, ipar, jac, psol)
end program main