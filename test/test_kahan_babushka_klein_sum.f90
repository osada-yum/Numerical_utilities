program test_kahan_babushka_klein_sum
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  integer(int64), parameter :: n = 10000000_int64
  real(real64), parameter :: v = 0.001_real64
  real(real64), parameter :: exact_sum = n * v
  real(real64), allocatable :: arr(:)
  real(real64) :: summ, summ_kahan_babushka_klein_sum
  type(kahan_summation) :: summ_k
  integer(int32) :: i
  allocate(arr(n))
  arr(:) = v
  ! call random_number(arr)
  summ = sum(arr)
  summ_kahan_babushka_klein_sum = kahan_babushka_klein_sum(arr)
  summ_k = kahan_summation(0.0_real64)
  do i = 1, n
     summ_k = summ_k%add(arr(i))
  end do
  write(error_unit, '(a)') "[test_kahan_babushka_klein_sum]"
  write(error_unit, '(g0)') exact_sum
  write(error_unit, '(g0)') abs(summ - exact_sum)
  write(error_unit, '(g0)') abs(summ_kahan_babushka_klein_sum - exact_sum)
  write(error_unit, '(g0)') abs(summ_k%val() - exact_sum)
end program test_kahan_babushka_klein_sum
