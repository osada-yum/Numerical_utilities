program test_online_variance
  use, intrinsic :: iso_fortran_env
  use online_variance_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_var = 1.0_real64 / 12
  integer(int32), parameter :: n = 10000000
  real(real64), allocatable :: arr(:)
  type(online_variance) :: ov
  integer(int32) :: i
  allocate(arr(n))
  call random_number(arr)
  ov = online_variance(arr)
  write(error_unit, '(a)') "[test_online_variance]"
  write(error_unit, '(g0)') abs(ov%mean() - exact_mean)
  write(error_unit, '(g0)') abs(sum(arr) / n - exact_mean)
  write(error_unit, '(g0)') abs(ov%var() - exact_var)
  write(error_unit, '(g0)') abs(ov%sample_var() - exact_var)
  write(error_unit, '(g0)') abs(sum((arr(:) - sum(arr) / n) ** 2) / n - exact_var)
end program test_online_variance
