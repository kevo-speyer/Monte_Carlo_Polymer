module com_vars
!Module with all global
#include "prepro.h"
implicit none
integer :: n_save, i_mon, init_mode, j_mon, n_mon, i_time, n_time, i_dim, n_dim, mv_mon
real(kind=8) , dimension(:,:), allocatable :: r0, r0_init
real(kind=8) , dimension(:), allocatable :: boundary, dr_mv, r_cm,r_cm_init, r_end
real(kind=8) :: Rend2, Temp, a_box, Lx, Ly, Lz, energy, delta_energy, old_energy, k_spr, inv_nmon
logical :: init_pos
#ifdef g3
real(kind=8) , dimension(:,:), allocatable :: r_cm_t
real(kind=8) , dimension(:), allocatable :: g3
integer :: j_save=1
#endif
end module com_vars
 
