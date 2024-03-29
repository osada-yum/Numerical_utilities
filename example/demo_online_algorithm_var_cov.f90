program demo_online_variance
  use, intrinsic :: iso_fortran_env
  use online_variance_covariance_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_var = 1.0_real64 / 12, exact_covar = 1.0_real64 / 12
  integer(int32), parameter :: n = 10000000
  real(real64), allocatable :: arr(:), arr2(:)
  type(online_variance_covariance) :: var_cov, var_cov_square
  integer(int32) :: i
  allocate(arr(n), arr2(n))
  call random_number(arr)
  call random_number(arr2)
  do i = 1, n
     arr(i) = real(i, real64) / n
  end do
  arr2(:) = arr(:)
  var_cov = online_variance_covariance(n, arr, arr2)
  var_cov_square = online_variance_covariance(n, arr ** 2, arr2 ** 2)
  write(error_unit, '(a)') "[test_online_variance_covariance]"
  associate(mean => sum(arr) / n)
    write(error_unit, '(a)') "mean of arr: "
    write(error_unit, '(g0)') abs(var_cov%mean1() - exact_mean)
    write(error_unit, '(*(g0, 1x))') abs(mean - exact_mean), abs(var_cov%mean1() - mean)
    associate(var => sum((arr(:) - sum(arr) / n) ** 2) / n, square_mean => sum(arr(:) ** 2) / n)
      write(error_unit, '(a)') "var of arr: "
      write(error_unit, '(g0)') abs(var_cov%var1() - exact_var)
      write(error_unit, '(*(g0, 1x))') abs(var - exact_var), abs(var_cov%var1() - var)
      write(error_unit, '(*(g0, 1x))') abs(square_mean - mean ** 2 - exact_var), &
           & abs(var_cov%var1() - (square_mean - mean ** 2))
    end associate
  end associate
  write(error_unit, '(a15, g0)') "sample_var1: ", abs(var_cov%sample_var1() - exact_var)

  associate(mean => sum(arr2) / n)
    write(error_unit, '(a)') "mean of arr2: "
    write(error_unit, '(g0)') abs(var_cov%mean2() - exact_mean)
    write(error_unit, '(*(g0, 1x))') abs(mean - exact_mean), abs(var_cov%mean2() - mean)
    associate(var => sum((arr2(:) - sum(arr2) / n) ** 2) / n, square_mean => sum(arr2(:) ** 2) / n)
      write(error_unit, '(a)') "var of arr2: "
      write(error_unit, '(g0)') abs(var_cov%var2() - exact_var)
      write(error_unit, '(*(g0, 1x))') abs(var - exact_var), abs(var_cov%var2() - var)
      write(error_unit, '(*(g0, 1x))') abs(square_mean - mean ** 2 - exact_var), &
           & abs(var_cov%var2() - (square_mean - mean ** 2))
    end associate
  end associate
  write(error_unit, '(a15, g0)') "sample_var2: ", abs(var_cov%sample_var2() - exact_var)

  associate(cov => sum((arr(:) - sum(arr) / n) * (arr2(:) - sum(arr2) / n))  / n, &
       & cov_another => sum(arr(:) * arr2(:)) / n - sum(arr) / n * sum(arr2) / n)
    write(error_unit, '(a)') "cov of arr, arr2: "
    write(error_unit, '(g0)') abs(var_cov%cov() - exact_covar)
    write(error_unit, '(g0)') var_cov%cov()
    write(error_unit, '(g0)') cov
    write(error_unit, '(g0)') cov_another
    write(error_unit, '(*(g0, 1x))') abs(cov - exact_covar), abs(var_cov%cov() - cov)
    write(error_unit, '(*(g0, 1x))') abs(cov_another - exact_covar), abs(var_cov%cov() - cov_another)
    write(error_unit, '(a15, g0)') "sample_cov: ", abs(var_cov%sample_cov() - exact_covar)
  end associate
end program demo_online_variance
