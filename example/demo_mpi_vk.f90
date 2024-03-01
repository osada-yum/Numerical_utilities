program demo_mpi_vk
  use, intrinsic :: iso_fortran_env
  use mpi
  use variance_kahan_m
  use mpi_variance_kahan_m
  implicit none
  real(real64), parameter :: exact_mean = 0.5_real64, exact_square_mean = 1.0_real64 / 3, &
       & exact_var = 1.0_real64 / 12, exact_covar = 0.0_real64
  integer(int32), parameter :: n = 10000000
  integer(int32) :: myrank, num_proc, ierr
  real(real64) :: arr(n)
  type(variance_kahan) :: var, recv_vk
  integer(int32) :: i
  call MPI_Init(ierr)
  call MPI_Comm_Rank(MPI_COMM_WORLD, myrank, ierr)
  call MPI_Comm_Size(MPI_COMM_WORLD, num_proc, ierr)
  write(error_unit, '(i0, a, i0)') myrank, " / ", num_proc
  call random_number(arr)
  var = variance_kahan()
  do i = 1, n
     call var%add_data(arr(i))
  end do
  write(error_unit, '(a, i0, a)') "[demo_mpi_vk: ", myrank, "]"
  write(error_unit, '(*(g0, 1x))') myrank, var%num_sample()
  write(error_unit, '(*(g0, 1x))') myrank, var%mean(), abs(var%mean() - exact_mean)
  write(error_unit, '(*(g0, 1x))') myrank, var%square_mean(), abs(var%square_mean() - exact_square_mean)
  write(error_unit, '(*(g0, 1x))') myrank, var%var(), abs(var%var() - exact_var)
  call vk_mpi_gather(var, recv_vk, 0, myrank, num_proc, ierr)
  if (myrank == 0) then
     write(error_unit, '(a)') "[demo_mpi_vk]"
     write(error_unit, '(*(g0, 1x))') myrank, recv_vk%num_sample()
     write(error_unit, '(*(g0, 1x))') recv_vk%mean(), abs(recv_vk%mean() - exact_mean)
     write(error_unit, '(*(g0, 1x))') recv_vk%square_mean(), abs(recv_vk%square_mean() - exact_square_mean)
     write(error_unit, '(*(g0, 1x))') recv_vk%var(), abs(recv_vk%var() - exact_var)
  end if
  call MPI_Finalize(ierr)
end program demo_mpi_vk
