module simulation
    use adsorption_models
    use types
    use simulation_parameters
    implicit none
    contains
        subroutine ddaspk_integration(model, neq, y, t0, tf, rtol, atol, rpar, ipar, jac, psol)
            integer, intent(in) :: neq !number of equations
            integer, intent(in), dimension(:) :: ipar !integer parameters used in res
            real(dp), intent(in), dimension(:) :: rpar !integer parameters used in res
            real(dp), intent(in) :: t0, tf
            real(dp), intent(inout), dimension(neq) :: rtol
            real(dp), intent(inout), dimension(neq) :: atol
            real(dp), intent(inout), dimension(neq) :: y

            integer, dimension(20) :: info 
            integer :: lrw
            integer :: liw
            integer :: i
            integer :: idid, ires
            real(dp), dimension(neq) :: delta
            real(dp), dimension(neq) :: yprime
            real(dp), allocatable, dimension(:) :: rwork
            integer, allocatable, dimension(:) :: iwork
            real(dp) :: t
            real(dp) :: cj
            real(dp) :: tout
            external jac, psol

            lrw = 50 + 9*neq + neq*neq
            liw = 40 + 3*neq

            allocate(rwork(lrw), iwork(liw))
            
                        ! inicializando par�metros
            iwork = 0
            rwork = 0.0_dp
            ! integrator time parameters
            t = t0          ! tempo adimensional inicial
            ! inicializando as derivadas e zerando os Deltas
            delta  = 0.0_dp
            yprime = 0.0_dp
    
            call model(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
            yprime = delta
            call model(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
           ! do i = 1, nsteps
                where(y < 1.e-20_dp) y = 0.0_dp
                ! define tempo final deste passo de integra��o
                tout = tf! real(i, dp)*dt
                ! passando par�metros da integra��o
                info = 0                ! zera os info
                info(2) = 1             ! para que atol e rtol sejam vetores
               ! info(10) = 1            ! se = 3, garante que as vari�veis sejam positivas
                iwork(40+1:40+neq) = 1  !if Y(I) must be .GE. 0,
                info(4) = 1           ! para n�o ultrapassar tout
                info(3) = 1
                RWORK(1) = tf        ! max tout
                ! chama rotina de integra��
                i=1
                write(*,*) 't, ', 'cb'
                do
                    call ddaspk(model,neq,t,y,yprime,tout,info,rtol,atol,idid,&
                    rwork,lrw,iwork,liw,rpar,ipar,jac,psol)
                    write(*,*) t, y(neq)
                    i = i + 1
                    if (IDID < -1) exit
                    if (t == tout) exit
                    ! se faltaram itera��es internas, volta e continua a integra��o
                    INFO(1) = 1
                end do
                ! verifca se houve erro de integra��o
                if (idid < 0) then
                    write(*,*) 'IDID = ', idid
                    !exit
                end if
                !print '(1f6.2,"% finished.")',float(i)/float(Ndado)*100
          !  end do
        end subroutine ddaspk_integration

        subroutine model(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
            integer, intent(in) :: neq
            real(dp), intent(in) :: t
            real(dp), intent(in), dimension(neq) :: y
            real(dp), intent(in), dimension(neq) :: yprime
            real(dp), intent(in) :: cj
            real(dp), intent(out), dimension(neq) :: delta
            integer, intent(in) :: ires
            real(dp), intent(in), dimension(:) :: rpar
            integer, intent(in), dimension(:) :: ipar
            
            select case(model_id)
                case(1)
                    call fd_eta_3pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
                case(2)
                    call fd_eta_5pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
                case(3)
                    call fd_mu_3pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
                case(4)
                    call fd_mu_5pts(neq, t, y, yprime, cj, delta, ires, rpar, ipar)
                case default
                    write(*,*) 'MODEL ID: BAD ARGUMENT' 
            end select
        end subroutine model
end module simulation