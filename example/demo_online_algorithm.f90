program test_online_variance
  use, intrinsic :: iso_fortran_env
  use online_variance_m
  use online_covariance_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_var = 1.0_real64 / 12, exact_covar = 0.0_real64
  integer(int32), parameter :: n = 10000000
  real(real64), allocatable :: arr(:), arr2(:)
  type(online_variance) :: o_var
  type(online_covariance) :: o_cov
  integer(int32) :: i
  allocate(arr(n), arr2(n))
  call random_number(arr)
  call random_number(arr2)
  o_var = online_variance(arr)
  do i = 1, n - 1
     call o_cov%add_pair(arr(i), arr2(i))
     ! write(output_unit, '(a15, 2(g0, 1x))') "sample_covar: ", i, o_cov%sample_covar()
  end do
  write(error_unit, '(a)') "[test_online_variance]"
  write(error_unit, '(g0)') abs(o_var%mean() - exact_mean)
  write(error_unit, '(g0)') abs(sum(arr) / n - exact_mean)
  write(error_unit, '(g0)') abs(o_var%var() - exact_var)
  write(error_unit, '(g0)') abs(o_var%sample_var() - exact_var)
  write(error_unit, '(g0)') abs(sum((arr(:) - sum(arr) / n) ** 2) / n - exact_var)
  write(error_unit, '(a)') "[test_online_covariance]"
  write(error_unit, '(a15, g0)') "mean1: ", abs(o_cov%mean1() - exact_mean)
  write(error_unit, '(a15, g0)') "mean2: ", abs(o_cov%mean2() - exact_mean)
  write(error_unit, '(a15, g0)') "var1: ", abs(o_cov%var1() - exact_var)
  write(error_unit, '(a15, g0)') "sample_var1: ", abs(o_cov%sample_var1() - exact_var)
  write(error_unit, '(a15, g0)') "var2: ", abs(o_cov%var2() - exact_var)
  write(error_unit, '(a15, g0)') "sample_var2: ", abs(o_cov%sample_var2() - exact_var)
  write(error_unit, '(a15, g0)') "covar: ", abs(o_cov%covar() - exact_covar)
  write(error_unit, '(a15, g0)') "sample_covar: ", abs(o_cov%sample_covar() - exact_covar)
end program test_online_variance
