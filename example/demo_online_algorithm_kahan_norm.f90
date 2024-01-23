program test_online_variance
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  use online_variance_m
  use online_variance_kahan_m
  implicit none
  real(real64), parameter :: exact_mean = 0.0_real64, exact_var = 1.0_real64
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int32), parameter :: n = 10000000
  real(real64), allocatable :: arr(:), y(:), x(:)
  type(online_variance) :: ov
  type(online_variance_kahan) :: ovk
  integer(int32) :: i
  allocate(y(n/2), x(n/2))
  call random_number(y)
  call random_number(x)
  allocate(arr(n))
  arr(1:n/2) = sqrt(- 2 * log(x)) * cos(2 * pi * y)
  arr(n/2+1:n) = sqrt(- 2 * log(x)) * sin(2 * pi * y)
  ov = online_variance(arr)
  ovk = online_variance_kahan(arr)
  write(error_unit, '(a)') "[test_online_variance_kahan_norm]"
  write(error_unit, '(g0)') abs(ov%mean() - exact_mean)
  write(error_unit, '(g0)') abs(ovk%mean() - exact_mean)
  write(error_unit, '(g0)') abs(sum(arr) / n - exact_mean)
  write(error_unit, '(g0)') abs(sum((arr(:) - sum(arr) / n) ** 2) / n - exact_var)
  write(error_unit, '(g0)') abs(ov%var() - exact_var)
  write(error_unit, '(g0)') abs(ovk%var() - exact_var)
  write(error_unit, '(g0)') abs(sum((arr(:) - sum(arr) / n) ** 2) / (n - 1) - exact_var)
  write(error_unit, '(g0)') abs(ov%sample_var() - exact_var)
  write(error_unit, '(g0)') abs(ovk%sample_var() - exact_var)
end program test_online_variance
