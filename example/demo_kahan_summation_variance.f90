program test_kahan_summation_variance
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_var = 1.0_real64 / 12
  integer(int64), parameter :: n = 1000000_int64
  real(real64), allocatable :: arr(:)
  real(real64) :: summ, summ2
  real(real64) :: summ_kahan_babushka_sum, summ_kahan_babushka_sum2
  real(real64) :: summ_kahan_babushka_klein_sum, summ_kahan_babushka_klein_sum2
  real(real64) :: summ_diff, summ_diff_kahan_babushka_sum, summ_diff_kahan_babushka_klein_sum
  type(kahan_summation) :: summ_diff_k
  type(kahan_summation) :: summ_k, summ_k2
  integer(int32) :: i
  allocate(arr(n))
  call random_number(arr)
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
  summ_kahan_babushka_sum = kahan_babushka_sum(arr)
  summ_kahan_babushka_sum2 = kahan_babushka_sum(arr ** 2)
  summ_kahan_babushka_klein_sum = kahan_babushka_klein_sum(arr)
  summ_kahan_babushka_klein_sum2 = kahan_babushka_klein_sum(arr ** 2)

  summ_diff = 0.0_real64
  summ_diff_k = kahan_summation(0.0_real64)
  associate(mean => summ / n, mean_k => summ_k%val() / n, &
       & mean_summ_kahan_babushka_sum => summ_kahan_babushka_sum / n, &
       & mean_summ_kahan_babushka_klein_sum => summ_kahan_babushka_klein_sum / n)
    do i = 1, n
       summ_diff = summ_diff + (arr(i) - mean) ** 2
       summ_diff_k = summ_diff_k%add((arr(i) - mean_k) ** 2)
    end do
    summ_diff_kahan_babushka_sum = kahan_babushka_sum((arr(:) - mean_summ_kahan_babushka_sum) ** 2)
    summ_diff_kahan_babushka_klein_sum = kahan_babushka_klein_sum((arr(:) - mean_summ_kahan_babushka_sum) ** 2)
    ! write(error_unit, '(a)') "[test_kahan_summation_variance]"
    ! write(error_unit, '(a)') " [simple variance]"
    ! write(error_unit, '(g0)') abs(summ2 / n - (summ / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ2 * n - summ ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs(summ_diff / n - exact_var)
    ! write(error_unit, '(a)') " [kahan variance]"
    ! write(error_unit, '(g0)') abs(summ_diff_k%val() / n - exact_var)
    ! write(error_unit, '(g0)') abs(summ_k2%val() / n - (summ_k%val() / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ_k2%val() * n - summ_k%val() ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(a)') " [kahan_sum variance]"
    ! write(error_unit, '(g0)') abs(summ_diff_kahan_babushka_sum / n - exact_var)
    ! write(error_unit, '(g0)') abs(summ_diff_kahan_babushka_klein_sum / n - exact_var)
    ! write(error_unit, '(g0)') abs(summ_kahan_babushka_sum2 / n - (summ_kahan_babushka_sum / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ_kahan_babushka_sum2 * n - summ_kahan_babushka_sum ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs(summ_kahan_babushka_klein_sum2 / n - (summ_kahan_babushka_klein_sum / n) ** 2 - exact_var)
    ! write(error_unit, '(g0)') abs((summ_kahan_babushka_klein_sum2 * n - summ_kahan_babushka_klein_sum ** 2) / n ** 2 - exact_var)
    ! write(error_unit, '(a)') "[test_kahan_summation_mean]"
    ! write(error_unit, '(g0)') abs(mean - exact_mean)
    ! write(error_unit, '(g0)') abs(mean_k - exact_mean)
    ! write(error_unit, '(g0)') abs(mean_summ_kahan_babushka_sum - exact_mean)
    ! write(error_unit, '(g0)') abs(mean_summ_kahan_babushka_klein_sum - exact_mean)
  end associate
end program test_kahan_summation_variance
