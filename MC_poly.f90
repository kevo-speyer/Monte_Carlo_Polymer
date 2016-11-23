program MC_chain
use com_vars
use ziggurat

implicit none

! Read input parameters
call read_input()

! Initialize system
call init_system()

! Meassure Energy
call get_energy()

!begin time loop
do i_time = 1, n_time
    
    old_energy = energy

!   Try displacement
    call move_attempt()
    
    call get_energy()

! Accept move with proba = min(1,exp( -delta_E/kT ) )
    call MC_acceptance()

    call save_data()

end do

end program


subroutine save_data()
use com_vars
use ziggurat

implicit none

write(60,*) i_time, energy

end subroutine

subroutine MC_acceptance()
use com_vars
use ziggurat

implicit none
real(kind=8) :: prob_rej

!Difference energy between old and new config
delta_energy = energy - old_energy

! If energy lowers, accept move
!if (delta_energy .gt. 0.) then
    prob_rej = 1. - exp(-delta_energy / Temp)
    if( uni() .le. prob_rej) then !If change is rejected
        !Undo move
        r0(:,mv_mon) = r0(:,mv_mon) - dr_mv(:)
        energy = old_energy
    end if    
!end if
end subroutine


subroutine move_attempt()
use com_vars
use ziggurat

implicit none

!Pick random monomer
mv_mon = int( uni() * float(n_mon) ) + 1 

!Set random move
do i_dim = 1, n_dim
    dr_mv(i_dim) = a_box * uni() 
end do

!Perform move

r0(:,mv_mon) = r0(:,mv_mon) + dr_mv(:)

end subroutine

subroutine get_energy()
use com_vars
implicit none

energy = 0

!loop over monomers
do i_mon = 2, n_mon 
    energy = energy + k_spr * sqrt( sum( ( r0(:,i_mon) - r0(:,i_mon-1) ) ** 2 ) )    
end do

!write(61,*) energy
end subroutine

subroutine init_system()
use com_vars
use ziggurat
implicit none
real(kind=8) signo, dr
real(kind=8),dimension(n_dim) :: r_center

allocate( r0(n_dim,n_mon), boundary(n_dim), dr_mv(n_dim) ) 

!Initiate random seed
call init_rand_seed()

!Set box length
boundary(:) = Lx

!Set position of initial monomer
r_center(:) =  boundary(:) / 2.
r0(:,1) = r_center(:)

!Set position of remaining monomers
do i_mon = 2, n_mon
    do i_dim = 1, n_dim
        !Set distance to next monomer
        dr = rnor() + a_box
        !Set sign to move
        signo = uni()
        if (signo .lt. 0.5 ) then
            signo = -1.
        else
            signo = 1.
        end if
        r0( i_dim, i_mon ) = r0( i_dim, i_mon - 1 ) + dr * signo  
    end do
end do
 
end subroutine


subroutine read_input()
use com_vars
implicit none

open (unit = 53, file = "input.dat", status= "old")

read(53,*) n_dim
read(53,*) Lx !, Ly, Lz Square box
read(53,*) n_time
read(53,*) n_mon
read(53,*) a_box
read(53,*) k_spr
read(53,*) Temp
close(53)

end subroutine

      

subroutine init_rand_seed()
use ziggurat
logical :: file_ex
integer :: rand_seed
inquire(file='random_seed.dat',exist=file_ex)
if (file_ex) then
    open(unit = 68, file = "random_seed.dat", status= "old")
    read(68,*) rand_seed
    close(68)
else
    rand_seed = 465178
end if

call zigset( rand_seed )

open(unit = 68, file = "random_seed.dat", status="unknown")
write(68,*) shr3()
close(68)

end subroutine

