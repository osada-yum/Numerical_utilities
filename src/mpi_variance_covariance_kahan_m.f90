module mpi_variance_covariance_kahan_m
  use, intrinsic :: iso_fortran_env
  use mpi
  use kahan_summation_m
  use variance_covariance_kahan_m
  implicit none
contains
  impure subroutine vck_mpi_gather(send_vck, res_vck, root, myrank, num_proc, ierr)
    type(variance_covariance_kahan), intent(in) :: send_vck
    type(variance_covariance_kahan), intent(inout) :: res_vck
    integer(int32), intent(in) :: root, myrank, num_proc
    integer(int32), intent(out), optional :: ierr
    integer(int64) :: num_sample(num_proc)
    real(real64) :: v1(num_proc), v2(num_proc), v1_square(num_proc), v2_square(num_proc), v1v2(num_proc)
    integer(int32) :: i
    integer(int32) :: err(6)
    call MPI_Gather(send_vck%num_sample(), 1, MPI_INTEGER8, num_sample(1), 1, MPI_INTEGER8, root, MPI_COMM_WORLD, err(1))
    call MPI_Gather(send_vck%v1(),         1, MPI_REAL8,    v1(1),         1, MPI_REAL8,    root, MPI_COMM_WORLD, err(2))
    call MPI_Gather(send_vck%v2(),         1, MPI_REAL8,    v2(1),         1, MPI_REAL8,    root, MPI_COMM_WORLD, err(3))
    call MPI_Gather(send_vck%v1_square(),  1, MPI_REAL8,    v1_square(1),  1, MPI_REAL8,    root, MPI_COMM_WORLD, err(4))
    call MPI_Gather(send_vck%v2_square(),  1, MPI_REAL8,    v2_square(1),  1, MPI_REAL8,    root, MPI_COMM_WORLD, err(5))
    call MPI_Gather(send_vck%v1v2(),       1, MPI_REAL8,    v1v2(1),       1, MPI_REAL8,    root, MPI_COMM_WORLD, err(6))
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
       call res_vck%merge_data(variance_covariance_kahan(num_sample(i), v1(i), v2(i), v1_square(i), v2_square(i), v1v2(i)))
    end do
  end subroutine vck_mpi_gather
  impure subroutine vck_mpi_multi_gather(n, send_vck, res_vck, root, myrank, num_proc, ierr)
    integer(int32), intent(in) :: n
    type(variance_covariance_kahan), intent(in) :: send_vck(n)
    type(variance_covariance_kahan), intent(inout) :: res_vck(n)
    integer(int32), intent(in) :: root, myrank, num_proc
    integer(int32), intent(out), optional :: ierr
    integer(int64) :: num_sample(n)
    real(real64) :: v1(n), v2(n), &
         & v1_square(n), v2_square(n), v1v2(n)
    integer(int64), allocatable :: recv_num_sample(:)
    real(real64), allocatable :: recv_v1(:), recv_v2(:), &
         & recv_v1_square(:), recv_v2_square(:), recv_v1v2(:)
    integer(int32) :: i, j
    integer(int32) :: err(6)
    do j = 1, n
       num_sample(j) = send_vck(j)%num_sample()
       v1(j) = send_vck(j)%v1()
       v2(j) = send_vck(j)%v2()
       v1_square(j) = send_vck(j)%v1_square()
       v2_square(j) = send_vck(j)%v2_square()
       v1v2(j) = send_vck(j)%v1v2()
    end do
    if (myrank == root) then
       allocate(recv_num_sample(num_proc * n), &
            & recv_v1(num_proc * n), recv_v2(num_proc * n), &
            & recv_v1_square(num_proc * n), recv_v2_square(num_proc * n), &
            & recv_v1v2(num_proc * n))
    else
       allocate(recv_num_sample(1), &
            & recv_v1(1), recv_v2(1), &
            & recv_v1_square(1), recv_v2_square(1), &
            & recv_v1v2(1))
    end if
    call MPI_Gather(num_sample(1), n, MPI_INTEGER8, recv_num_sample(1), n, MPI_INTEGER8, root, MPI_COMM_WORLD, err(1))
    call MPI_Gather(v1(1),         n, MPI_REAL8,    recv_v1(1),         n, MPI_REAL8,    root, MPI_COMM_WORLD, err(2))
    call MPI_Gather(v2(1),         n, MPI_REAL8,    recv_v2(1),         n, MPI_REAL8,    root, MPI_COMM_WORLD, err(3))
    call MPI_Gather(v1_square(1),  n, MPI_REAL8,    recv_v1_square(1),  n, MPI_REAL8,    root, MPI_COMM_WORLD, err(4))
    call MPI_Gather(v2_square(1),  n, MPI_REAL8,    recv_v2_square(1),  n, MPI_REAL8,    root, MPI_COMM_WORLD, err(5))
    call MPI_Gather(v1v2(1),       n, MPI_REAL8,    recv_v1v2(1),       n, MPI_REAL8,    root, MPI_COMM_WORLD, err(6))
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
    do j = 1, n
       res_vck(j) = send_vck(j)
    end do
    do i = 1, num_proc
       if (i - 1 == root) cycle
       do j = 1, n
          associate(idx => (i - 1) * n + j)
            call res_vck(j)%merge_data(variance_covariance_kahan(&
                 & recv_num_sample(idx), recv_v1(idx), recv_v2(idx), &
                 & recv_v1_square(idx), recv_v2_square(idx), recv_v1v2(idx)))
          end associate
       end do
    end do
  end subroutine vck_mpi_multi_gather
end module mpi_variance_covariance_kahan_m
