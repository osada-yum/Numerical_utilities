program demo_online_variance_kahan_var_cov
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  use online_variance_m
  use online_variance_covariance_kahan_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_var = 1.0_real64 / 12, exact_covar = 0.0_real64
  integer(int32), parameter :: n = 10000000
  real(real64), allocatable :: arr(:), arr2(:)
  type(online_variance) :: o_var
  type(online_variance_covariance_kahan) :: var_cov, var_cov_square
  integer(int32) :: i
  allocate(arr(n), arr2(n))
  call random_number(arr)
  call random_number(arr2)
  arr2(:) = arr
  var_cov = online_variance_covariance_kahan(n, arr, arr2)
  o_var = online_variance(arr(:) * arr2(:))
  var_cov_square = online_variance_covariance_kahan(n, arr ** 2, arr2 ** 2)
  write(error_unit, '(a)') "[demo_online_variance_covariance_kahan]"
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
       & cov_another => sum(arr(:) * arr2(:)) / n - sum(arr) / n * sum(arr2) / n, &
       & cov_on => o_var%mean() - var_cov%mean1() * var_cov%mean2())
    write(error_unit, '(a)') "cov of arr, arr2: "
    write(error_unit, '(g0)') abs(var_cov%cov() - exact_covar)
    write(error_unit, '(g0)') var_cov%cov()
    write(error_unit, '(g0)') cov_on
    write(error_unit, '(g0)') cov
    write(error_unit, '(g0)') cov_another
    write(error_unit, '(*(g0, 1x))') abs(cov - exact_covar), abs(var_cov%cov() - cov)
    write(error_unit, '(*(g0, 1x))') abs(cov_on - exact_covar), abs(var_cov%cov() - cov_on)
    write(error_unit, '(*(g0, 1x))') abs(cov_another - exact_covar), abs(var_cov%cov() - cov_another)
    write(error_unit, '(a15, g0)') "sample_cov: ", abs(var_cov%sample_cov() - exact_covar)
  end associate

  block
    type(kahan_summation) :: v1, v2, v1_square, v2_square, v12
    type(kahan_summation) :: var1_twopass_kahan, var2_twopass_kahan
    v1 = kahan_summation(0.0_real64); v2 = kahan_summation(0.0_real64)
    v1_square = kahan_summation(0.0_real64); v2_square = kahan_summation(0.0_real64)
    v12 = kahan_summation(0.0_real64)
    do i = 1, n
       v1 = v1 + arr(i)
       v2 = v2 + arr2(i)
       v1_square = v1_square + arr(i) ** 2
       v2_square = v2_square + arr2(i) ** 2
       v12 = v12 + (arr(i) * arr2(i))
    end do
    var1_twopass_kahan = kahan_summation(0.0_real64)
    var2_twopass_kahan = kahan_summation(0.0_real64)
    do i = 1, n
       var1_twopass_kahan = var1_twopass_kahan + (arr(i) - v1%val() / n) ** 2
       var2_twopass_kahan = var2_twopass_kahan + (arr2(i) - v2%val() / n) ** 2
    end do
    write(error_unit, '(a)') "[using kahan_summation]"
    associate(mean1 => v1%val() / n, mean2 => v2%val() / n)
      write(error_unit, '(a)') "mean of arr, arr2: "
      write(error_unit, '(*(g0, 1x))') abs(mean1 - var_cov%mean1()), abs(mean2 - var_cov%mean2())
    end associate
    associate(var1 => v1_square%val() / n - (v1%val() / n) ** 2, var2 => v2_square%val() / n - (v2%val() / n) ** 2, &
         & var1_twopass => sum((arr(:) - sum(arr(:)) / n) ** 2) / n, var2_twopass => sum((arr2(:) - sum(arr2(:)) / n) ** 2) / n)
      write(error_unit, '(a)') "var of arr, arr2: "
      write(error_unit, '(*(g0, 1x))') abs(var1 - var_cov%var1()), abs(var2 - var_cov%var2())
      write(error_unit, '(*(g0, 1x))') abs(var1_twopass - var_cov%var1()), abs(var2_twopass - var_cov%var2())
      write(error_unit, '(*(g0, 1x))') abs(var1_twopass - var1_twopass_kahan%val() / n), &
            abs(var2_twopass - var2_twopass_kahan%val() / n)
    end associate
    associate(cov => sum((arr(:) - v1%val() / n) * (arr2(:) - v2%val() / n))  / n, &
         & cov_another => v12%val() / n - (v1%val() / n) * (v2%val() / n))
      write(error_unit, '(a)') "cov of arr, arr2: "
      write(error_unit, '(g0)') abs(var_cov%cov() - exact_covar)
      write(error_unit, '(g0)') var_cov%cov()
      write(error_unit, '(g0)') cov
      write(error_unit, '(g0)') cov_another
      write(error_unit, '(*(g0, 1x))') abs(cov - exact_covar), abs(var_cov%cov() - cov)
      write(error_unit, '(*(g0, 1x))') abs(cov_another - exact_covar), abs(var_cov%cov() - cov_another)
      write(error_unit, '(*(g0, 1x))') abs(var_cov%cov() - cov), abs(var_cov%cov() - cov_another), abs(cov_another - cov)
      write(error_unit, '(a15, g0)') "sample_cov: ", abs(var_cov%sample_cov() - exact_covar)
    end associate
  end block
end program demo_online_variance_kahan_var_cov
