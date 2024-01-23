program test_kahan_summation_variance_norm
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  real(real64), parameter :: exact_mean = 0.0_real64, exact_var = 1.0_real64
  real(real64), parameter :: pi = 4 * atan(1.0_real64)
  integer(int64), parameter :: n = 100000_int64
  real(real64), allocatable :: arr(:), y(:), x(:)
  real(real64) :: summ, summ2
  real(real64) :: summ_kahan_sum, summ_kahan_sum2
  real(real64) :: summ_diff, summ_diff_kahan_sum
  type(kahan_summation) :: summ_diff_k
  type(kahan_summation) :: summ_k, summ_k2
  integer(int32) :: i
  allocate(y(n/2), x(n/2))
  call random_number(y)
  call random_number(x)
  allocate(arr(n))
  arr(1:n/2) = sqrt(- 2 * log(x)) * cos(2 * pi * y)
  arr(n/2+1:n) = sqrt(- 2 * log(x)) * sin(2 * pi * y)
  summ = 0.0_real64
  summ2 = 0.0_real64
  summ_k = kahan_summation(0.0_real64)
  summ_k2 = kahan_summation(0.0_real64)
  do i = 1, n
     summ = summ + arr(i)
     summ2 = summ2 + arr(i) ** 2
     summ_k = summ_k%add(arr(i))
     summ_k2 = summ_k2%add(arr(i) ** 2)
  end do
  summ_kahan_sum = kahan_babushka_sum(arr)
  summ_kahan_sum2 = kahan_babushka_sum(arr ** 2)

  summ_diff = 0.0_real64
  summ_diff_k = kahan_summation(0.0_real64)
  associate(mean => summ / n, mean_k => summ_k%val() / n, mean_summ_kahan_sum => summ_kahan_sum / n)
    do i = 1, n
       summ_diff = summ_diff + (arr(i) - mean) ** 2
       summ_diff_k = summ_diff_k%add((arr(i) - mean_k) ** 2)
    end do
    summ_diff_kahan_sum = kahan_babushka_sum((arr(:) - mean_summ_kahan_sum) ** 2)
    ! write(error_unit, '(a)') "[test_kahan_summation_variance_norm]"
    ! write(error_unit, '(a)') " [simple variance]"
    ! write(error_unit, '(g0)') abs(summ2 / n - (summ / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ2 * n - summ ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs(summ_diff / n - exact_var)
    ! write(error_unit, '(a)') " [kahan variance]"
    ! write(error_unit, '(g0)') abs(summ_k2%val() / n - (summ_k%val() / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ_k2%val() * n - summ_k%val() ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs(summ_diff_k%val() / n - exact_var)
    ! write(error_unit, '(a)') " [kahan_babushka_sum variance]"
    ! write(error_unit, '(g0)') abs(summ_kahan_sum2 / n - (summ_kahan_sum / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ_kahan_sum2 * n - summ_kahan_sum ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs(summ_diff_kahan_sum / n - exact_var)
    ! write(error_unit, '(a)') "[test_kahan_summation_mean_norm]"
    ! write(error_unit, '(g0)') abs(mean - exact_mean)
    ! write(error_unit, '(g0)') abs(mean_k - exact_mean)
    ! write(error_unit, '(g0)') abs(mean_summ_kahan_sum - exact_mean)
  end associate
end program test_kahan_summation_variance_norm
