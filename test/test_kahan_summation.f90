program test_kahan_summation
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  integer(int64), parameter :: n = 10000000_int64
  real(real64), parameter :: v = 0.01_real64
  real(real64), parameter :: exact_sum = n * v
  real(real64) :: summ
  type(kahan_summation) :: summ_k
  integer(int32) :: i
  summ = 0.0_real64
  summ_k = kahan_summation(0.0_real64)
  do i = 1, n
     summ = summ + v
     summ_k = summ_k%add(v)
  end do
  write(error_unit, '(a)') "[test_kahan_summation]"
  write(error_unit, '(g0)') abs(summ - exact_sum)
  write(error_unit, '(g0)') abs(summ_k%val() - exact_sum)
end program test_kahan_summation
