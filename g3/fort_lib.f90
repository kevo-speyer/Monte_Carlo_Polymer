program test
implicit none
integer :: i, n
real (kind=8) :: x
real (kind=8), allocatable,dimension(:) :: a, c

n=10
allocate (a(n),c(n))

do i=1,n
    x = (i-1)*2.14/(n-1) +1.
    a(i) = sin(x)
    print*,x,a(i)
end do

call f_routine(a,n,c)

print*,""
print*, c

end program

subroutine f_routine(x,n,c)
implicit none
real (kind=8), dimension(n), intent(in) :: x
real (kind=8), dimension(n) :: y
integer, intent(in) :: n
real (kind=8), dimension(n), intent(out)  :: c
integer ,dimension(n) :: n_corr
integer :: i,j
real (kind=8) :: x_mean = 0, x_var = 0

!Initialize counters to  0
c = 0. 
n_corr = 0


!First get mean value of data
do i=1,n
    x_mean = x_mean + x(i)
end do

x_mean = x_mean / float(n)


!Translate data so that it has 0 mean
y = x - x_mean

!Now calculate the variance and autocorrelation
do i=1,n
    x_var = x_var + y(i)**2 ! Accumulate deviations from mean
    do j=i,n
        c(j-i+1) = c(j-i+1) + y(i)*y(j) ! Accumulate covariance for each lag
        n_corr(j-i+1) = n_corr(j-i+1) + 1 !Count cases for each lag
    end do
end do

x_var = x_var / float(n) ! Normalize Variance

do i=1,n
    c(i) = c(i) / float(n_corr(i)) ! Get mean Covariance
    c(i) = c(i)  / x_var !Normalize with Variance to obtain correlation coefficient
end do

return

end subroutine f_routine







subroutine g_3(x,n,n_dim,c)
implicit none
real (kind=8), dimension(n_dim, n), intent(in) :: x
real (kind=8), dimension(n) :: y
integer, intent(in) :: n, n_dim
real (kind=8), dimension(n), intent(out)  :: c
integer ,dimension(n) :: n_corr
integer :: i,j
real (kind=8) :: x_mean = 0, x_var = 0

!Initialize counters to  0
c = 0. 
n_corr = 0


!First get mean value of data
!do i=1,n
!    x_mean = x_mean + x(i)
!end do
!
!x_mean = x_mean / float(n)


!Translate data so that it has 0 mean
!y = x - x_mean

!Now calculate the variance and autocorrelation
do i=1,n-1
    !x_var = x_var + y(i)**2 ! Accumulate deviations from mean
    do j=i+1,n
        c(j-i) =  sum( (x(:,i) - x(:,j) )**2 ) !c(j-i+1) + y(i)*y(j) ! Accumulate covariance for each lag
        n_corr(j-i) = n_corr(j-i) + 1 !Count cases for each lag
    end do
end do

!x_var = x_var / float(n) ! Normalize Variance

do i=1,n-1
    c(i) = c(i) / float(n_corr(i)) ! Get mean Covariance
    !c(i) = c(i)  / x_var !Normalize with Variance to obtain correlation coefficient
end do

return

end subroutine g_3
