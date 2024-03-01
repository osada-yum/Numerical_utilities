module mpi_variance_kahan_m
  use, intrinsic :: iso_fortran_env
  use mpi
  use kahan_summation_m
  use variance_kahan_m
  implicit none
contains
  impure subroutine vk_mpi_gather(send_vk, res_vk, root, myrank, num_proc, ierr)
    type(variance_kahan), intent(in) :: send_vk
    type(variance_kahan), intent(inout) :: res_vk
    integer(int32), intent(in) :: root, myrank, num_proc
    integer(int32), intent(out), optional :: ierr
    integer(int64) :: num_sample(num_proc)
    real(real64) :: v1(num_proc), v2(num_proc), v1_square(num_proc), v2_square(num_proc), v1v2(num_proc)
    integer(int32) :: i
    integer(int32) :: err(3)
    call MPI_Gather(send_vk%num_sample(), 1, MPI_INTEGER8, num_sample(1), 1, MPI_INTEGER8, root, MPI_COMM_WORLD, err(1))
    call MPI_Gather(send_vk%v1(),         1, MPI_REAL8,    v1(1),         1, MPI_REAL8,    root, MPI_COMM_WORLD, err(2))
    call MPI_Gather(send_vk%v1_square(),  1, MPI_REAL8,    v1_square(1),  1, MPI_REAL8,    root, MPI_COMM_WORLD, err(3))
    if (any(err(:) /= 0)) then
       if (present(ierr)) then
          ierr = 0
          do i = 1, size(err)
             if (err(i) /= 0) &
                  ierr = ibset(ierr, i - 1)
          end do
       else
          write(error_unit, '(a)') "Error in `vk_mpi_gather`"
          write(error_unit, '(*(i0, 1x))') err(:)
          error stop 1
       end if
       return
    end if
    write(error_unit, *) "send: ", myrank, send_vk%var()
    if (myrank /= root) return
    res_vk = send_vk
    do i = 1, num_proc
       if (i - 1 == root) cycle
       call res_vk%merge_data(variance_kahan(num_sample(i), v1(i), v1_square(i)))
    end do
    write(error_unit, *) "recv: ", myrank, res_vk%var()
  end subroutine vk_mpi_gather
  impure subroutine vk_mpi_multi_gather(n, send_vk, res_vk, root, myrank, num_proc, ierr)
    integer(int32), intent(in) :: n
    type(variance_kahan), intent(in) :: send_vk(n)
    type(variance_kahan), intent(inout) :: res_vk(n)
    integer(int32), intent(in) :: root, myrank, num_proc
    integer(int32), intent(out), optional :: ierr
    integer(int64) :: num_sample(n)
    real(real64) :: v1(n), v2(n), &
         & v1_square(n), v2_square(n), v1v2(n)
    integer(int64), allocatable :: recv_num_sample(:)
    real(real64), allocatable :: recv_v1(:), recv_v2(:), &
         & recv_v1_square(:), recv_v2_square(:), recv_v1v2(:)
    integer(int32) :: i, j
    integer(int32) :: err(3)
    do j = 1, n
       num_sample(j) = send_vk(j)%num_sample()
       v1(j) = send_vk(j)%v1()
       v1_square(j) = send_vk(j)%v1_square()
    end do
    if (myrank == root) then
       allocate(recv_num_sample(num_proc * n), recv_v1(num_proc * n), recv_v1_square(num_proc * n))
    else
       allocate(recv_num_sample(1), recv_v1(1), recv_v1_square(1))
    end if
    call MPI_Gather(num_sample(1), n, MPI_INTEGER8, recv_num_sample(1), n, MPI_INTEGER8, root, MPI_COMM_WORLD, err(1))
    call MPI_Gather(v1(1),         n, MPI_REAL8,    recv_v1(1),         n, MPI_REAL8,    root, MPI_COMM_WORLD, err(2))
    call MPI_Gather(v1_square(1),  n, MPI_REAL8,    recv_v1_square(1),  n, MPI_REAL8,    root, MPI_COMM_WORLD, err(3))
    if (any(err(:) /= 0)) then
       if (present(ierr)) then
          ierr = 0
          do i = 1, size(err)
             if (err(i) /= 0) &
                  ierr = ibset(ierr, i - 1)
          end do
       else
          write(error_unit, '(a)') "Error in `vk_mpi_gather`"
          write(error_unit, '(*(i0, 1x))') err(:)
          error stop 1
       end if
       return
    end if
    if (myrank /= root) return
    do j = 1, n
       res_vk(j) = send_vk(j)
    end do
    do i = 1, num_proc
       if (i - 1 == root) cycle
       do j = 1, n
          associate(idx => (i - 1) * n + j)
            call res_vk(j)%merge_data(variance_kahan(&
                 & recv_num_sample(idx), recv_v1(idx), recv_v1_square(idx)))
          end associate
       end do
    end do
  end subroutine vk_mpi_multi_gather
end module mpi_variance_kahan_m
