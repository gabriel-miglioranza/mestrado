module simulation_parameters
    use types
    implicit none
    real(dp) :: dx
    real(dp) :: bi
    real(dp) :: alfa
    real(dp), allocatable, dimension(:) :: x
    real(dp), allocatable, dimension(:,:) :: a
    real(dp), allocatable, dimension(:,:) :: b
    integer :: s
    integer :: n
    integer :: cc1
    integer :: cc2
    integer :: model_id
end module simulation_parameters