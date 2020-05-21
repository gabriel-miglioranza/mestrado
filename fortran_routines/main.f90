program main
    use orthogonal_collocation
    implicit none
    real(8) :: Raiz(20)
    integer :: points = 18
    real(8) :: x_mine(18)
    real(8) :: A(20), B(20), W(20)
    real(8) :: INF = 0.0
    real(8) :: SUP = 1.0
    integer :: i

    A = 0.0
    B = 0.0
    W = 0.0
    Raiz = 0.0
    call colocacao(points,INF,SUP,A,B,Raiz,W)
    !write(*,*) raiz(2:points+1)
    call find_collocation_points(x_mine, points)
    !write(*,*) x_mine
    write(*,*) abs(x_mine - raiz(2:points+1))
end program main