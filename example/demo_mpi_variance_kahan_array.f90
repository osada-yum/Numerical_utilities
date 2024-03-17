program demo_variance_covariance_kahan_m
  use, intrinsic :: iso_fortran_env
  use mpi
  use kahan_summation_m
  use mpi_variance_kahan_array_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_square_mean = 1.0_real64 / 3, &
       & exact_var = 1.0_real64 / 12, exact_covar = 0.0_real64
  integer(int32), parameter :: n = 1000000, m = 10
  real(real64), allocatable :: arr(:, :)
  integer(int32) :: myrank, num_proc, ierr
  type(mpi_variance_kahan_array) :: var_array
  integer(int32) :: i, j
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)
  write(error_unit, '(i0, a, i0)') myrank, " / ", num_proc
  allocate(arr(n, m))
  call random_number(arr)
  var_array = mpi_variance_kahan_array(m)
  do j = 1, m
     do i = 1, n
        call var_array%add_data(j, arr(i, j))
     end do
  end do
  write(error_unit, '(a, i0)') "[demo_mpi_variance_kahan_array]", myrank
  do j = 1, m
     associate(mean => sum(arr(:, j)) / n, square_mean => sum(arr(:, j) ** 2) / n)
       write(error_unit, '(a, i0)') "mean of arr: ", j
       write(error_unit, '(*(g0, 1x))') var_array%mean(j), abs(var_array%mean(j) - exact_mean)
       write(error_unit, '(*(g0, 1x))') abs(mean - exact_mean), abs(var_array%mean(j) - mean)
       write(error_unit, '(*(g0, 1x))') var_array%square_mean(j), abs(var_array%square_mean(j) - exact_square_mean)
       write(error_unit, '(*(g0, 1x))') abs(square_mean - exact_square_mean), abs(var_array%square_mean(j) - square_mean)
       associate(var => sum((arr(:, j) - mean) ** 2) / n)
         write(error_unit, '(a, i0)') "var of arr: ", j
         write(error_unit, '(g0)') abs(var_array%var(j) - exact_var)
         write(error_unit, '(*(g0, 1x))') abs(var - exact_var), abs(var_array%var(j) - var)
         write(error_unit, '(*(g0, 1x))') abs(square_mean - mean ** 2 - exact_var), &
              & abs(var_array%var(j) - (square_mean - mean ** 2))
       end associate
       write(error_unit, '(a15, g0)') "sample_var: ", abs(var_array%sample_var(j) - exact_var)
     end associate
  end do
  block
    type(mpi_variance_kahan_array) :: other, other2
    other = mpi_variance_kahan_array(var_array%num_sample([(i, i = 1, m)]), var_array%get_vs())
    other2 = other
    write(error_unit, '(a)') "[demo_copy]"
    write(error_unit, '(a)') "other"
    do j = 1, m
       write(error_unit, '(*(g0, 1x))') myrank, abs(var_array%v1(j) - other%v1(j)), &
            & abs(var_array%var(j) - other%var(j))
    end do
    call other2%merge_data(var_array)
    write(error_unit, '(a)') "other2"
    do j = 1, m
       write(error_unit, '(*(g0, 1x))') myrank, abs(2 * var_array%v1(j) - other2%v1(j)), &
            & abs(var_array%var(j) - other2%var(j))
    end do
  end block
  block
    type(mpi_variance_kahan_array) :: res
    call var_array%mpi_sum(myrank, num_proc, 0, res, ierr)
    if (myrank == 0) then
       do j = 1, m
          write(error_unit, '(*(g0, 1x))') "res: ", var_array%v1(j), res%v1(j), &
               & var_array%v1_square(j), res%v1_square(j)
          write(error_unit, '(*(g0, 1x))') abs(res%mean(j) - exact_mean), &
               & abs(res%var(j) - exact_var)
       end do
    end if
  end block
  call MPI_Finalize(ierr)
end program demo_variance_covariance_kahan_m
