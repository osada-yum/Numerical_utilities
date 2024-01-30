program demo_mpi_vck
  use, intrinsic :: iso_fortran_env
  use mpi
  use variance_covariance_kahan_m
  use mpi_variance_covariance_kahan_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_square_mean = 1.0_real64 / 3, &
       & exact_var = 1.0_real64 / 12, exact_covar = 0.0_real64
  integer(int32), parameter :: n = 10000000
  integer(int32) :: myrank, num_proc, ierr
  real(real64) :: arr(n), arr2(n)
  type(variance_covariance_kahan) :: var_cov, recv_vck
  integer(int32) :: i
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)
  write(error_unit, '(i0, a, i0)') myrank, " / ", num_proc
  call random_number(arr)
  call random_number(arr2)
  var_cov = variance_covariance_kahan()
  do i = 1, n
     call var_cov%add_data(arr(i), arr2(i))
  end do
  write(error_unit, '(a, i0, a)') "[demo_mpi_vck: ", myrank, "]"
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%num_sample()
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%mean1(), abs(var_cov%mean1() - exact_mean)
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%mean2(), abs(var_cov%mean2() - exact_mean)
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%square_mean1(), abs(var_cov%square_mean1() - exact_square_mean)
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%square_mean2(), abs(var_cov%square_mean2() - exact_square_mean)
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%var1(), abs(var_cov%var1() - exact_var)
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%var2(), abs(var_cov%var2() - exact_var)
  write(error_unit, '(*(g0, 1x))') myrank, var_cov%cov(), abs(var_cov%cov() - exact_covar)
  call vck_mpi_gather(var_cov, recv_vck, 0, myrank, num_proc, ierr)
  if (myrank == 0) then
     write(error_unit, '(a)') "[demo_mpi_vck]"
     write(error_unit, '(*(g0, 1x))') myrank, recv_vck%num_sample()
     write(error_unit, '(*(g0, 1x))') recv_vck%mean1(), abs(recv_vck%mean1() - exact_mean)
     write(error_unit, '(*(g0, 1x))') recv_vck%mean2(), abs(recv_vck%mean2() - exact_mean)
     write(error_unit, '(*(g0, 1x))') recv_vck%square_mean1(), abs(recv_vck%square_mean1() - exact_square_mean)
     write(error_unit, '(*(g0, 1x))') recv_vck%square_mean2(), abs(recv_vck%square_mean2() - exact_square_mean)
     write(error_unit, '(*(g0, 1x))') recv_vck%var1(), abs(recv_vck%var1() - exact_var)
     write(error_unit, '(*(g0, 1x))') recv_vck%var2(), abs(recv_vck%var2() - exact_var)
     write(error_unit, '(*(g0, 1x))') recv_vck%cov(), abs(recv_vck%cov() - exact_covar)
  end if
  call MPI_Finalize(ierr)
end program demo_mpi_vck
