module com_vars
!Module with all global
implicit none
integer :: n_save, i_mon, init_mode, j_mon, n_mon, i_time, n_time, i_dim, n_dim, mv_mon
real(kind=8) , dimension(:,:), allocatable :: r0
real(kind=8) , dimension(:), allocatable :: boundary, dr_mv, r_cm, r_end
real(kind=8) :: Rend2, Temp, a_box, Lx, Ly, Lz, energy, delta_energy, old_energy, k_spr, inv_nmon
logical :: init_pos
end module com_vars
 
