module mpi_variance_covariance_kahan_m
  use, intrinsic :: iso_fortran_env
  use mpi_f08
  use kahan_summation_m
  use variance_covariance_kahan_m
  implicit none
contains
  impure subroutine vck_mpi_gather(send_vck, res_vck, root, myrank, num_proc, ierr)
    type(variance_covariance_kahan), intent(in) :: send_vck
    type(variance_covariance_kahan), intent(inout) :: res_vck
    integer(int32), intent(in) :: root, myrank, num_proc
    integer(int32), intent(out), optional :: ierr
    integer(int64) :: n(num_proc)
    real(real64) :: v1(num_proc), v2(num_proc), v1_square(num_proc), v2_square(num_proc), v1v2(num_proc)
    integer(int32) :: i
    integer(int32) :: err(6)
    call MPI_Gather(send_vck%num_sample(), 1, MPI_INTEGER8, n(1),         1, MPI_INTEGER8, root, MPI_COMM_WORLD, err(1))
    call MPI_Gather(send_vck%v1(),         1, MPI_REAL8,    v1(1),        1, MPI_REAL8,    root, MPI_COMM_WORLD, err(2))
    call MPI_Gather(send_vck%v2(),         1, MPI_REAL8,    v2(1),        1, MPI_REAL8,    root, MPI_COMM_WORLD, err(3))
    call MPI_Gather(send_vck%v1_square(),  1, MPI_REAL8,    v1_square(1), 1, MPI_REAL8,    root, MPI_COMM_WORLD, err(4))
    call MPI_Gather(send_vck%v2_square(),  1, MPI_REAL8,    v2_square(1), 1, MPI_REAL8,    root, MPI_COMM_WORLD, err(5))
    call MPI_Gather(send_vck%v1v2(),       1, MPI_REAL8,    v1v2(1),      1, MPI_REAL8,    root, MPI_COMM_WORLD, err(6))
    if (any(err(:) /= 0)) then
       if (present(ierr)) then
          ierr = 0
          do i = 1, size(err)
             if (err(i) /= 0) &
                  ierr = ibset(ierr, i - 1)
          end do
       else
          write(error_unit, '(a)') "Error in `vck_mpi_gather`"
          write(error_unit, '(*(i0, 1x))') err(:)
          error stop 1
       end if
       return
    end if
    if (myrank /= root) return
    res_vck = send_vck
    do i = 1, num_proc
       if (i - 1 == root) cycle
       call res_vck%merge_data(variance_covariance_kahan(n(i), v1(i), v2(i), v1_square(i), v2_square(i), v1v2(i)))
    end do
  end subroutine vck_mpi_gather
end module mpi_variance_covariance_kahan_m
