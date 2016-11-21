module com_vars
!Module with all global
implicit none
integer :: i_mon, j_mon, n_mon, i_time, n_time, i_dim, n_dim, mv_mon
real(kind=8) , dimension(:,:), allocatable :: r0
real(kind=8) , dimension(:), allocatable :: boundary, dr_mv
real(kind=8) :: Temp, a_box, Lx, Ly, Lz, energy, delta_energy, old_energy, k_spr
end module com_vars
 
