module mpi_variance_kahan_array_m
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  use mpi
  implicit none
  private
  public mpi_variance_kahan_array
  type :: mpi_variance_kahan_array
     private
     integer(int32) :: array_size_
     integer(int64), allocatable :: n_(:)
     type(kahan_summation), allocatable :: vs_(:, :)
     !> n_: length(vs)
     !> vs_(1, i): sum(vs)
     !> vs_(2, i): sum(vs ** 2)
   contains
     procedure, pass :: add_data => add_data_mpi_vka
     procedure, pass :: merge_data => merge_data_mpi_vka
     !> getter
     procedure, pass :: num_sample => num_sample_mpi_vka
     procedure, pass :: v1 => v1_mpi_vka
     procedure, pass :: v1_square => v1_square_mpi_vka
     procedure, pass :: get_vs => get_vs_mpi_vka
     !> calculate
     procedure, pass :: mean => mean_mpi_vka
     procedure, pass :: square_mean => square_mean_mpi_vka
     procedure, pass :: var => var_mpi_vka
     procedure, pass :: sample_var => sample_var_mpi_vka
     !> mpi
     procedure, pass :: mpi_sum => mpi_sum_mpi_vka
  end type mpi_variance_kahan_array
  interface mpi_variance_kahan_array
     module procedure :: generate_mpi_variance_kahan_array
     module procedure :: generate_mpi_variance_kahan_array_by_data
  end interface mpi_variance_kahan_array
contains
  pure type(mpi_variance_kahan_array) function generate_mpi_variance_kahan_array(array_size) result(res)
    integer(int32), intent(in) :: array_size
    res%array_size_ = array_size
    allocate(res%n_(array_size), source = 0_int64)
    allocate(res%vs_(2, array_size), source = kahan_summation(0.0_real64))
  end function generate_mpi_variance_kahan_array
  pure subroutine add_data_mpi_vka(this, idx, v1)
    class(mpi_variance_kahan_array), intent(inout) :: this
    integer(int32), intent(in) :: idx
    real(real64), intent(in) :: v1
    if (idx < 1 .or. idx > this%array_size_) &
         & error stop "Error in `add_data(idx, v1)`."
    this%n_(idx) = this%n_(idx) + 1_int64
    this%vs_(1, idx) = this%vs_(1, idx) + v1
    this%vs_(2, idx) = this%vs_(2, idx) + v1 ** 2
  end subroutine add_data_mpi_vka
  pure type(mpi_variance_kahan_array) function generate_mpi_variance_kahan_array_by_data &
       & (ns, vs) result(res)
    integer(int64), intent(in) :: ns(:)
    real(real64), intent(in) :: vs(:, :)
    integer(int32) :: s1, s2
    integer(int32) :: i
    s1 = size(vs, dim = 1)
    if (s1 /= 2) &
         & error stop "Error in `mpi_variance_kahan_array(n, vs)`: shape of vs(:, :) is illeagal."
    s2 = size(vs, dim = 2)
    if (size(ns) /= s2) &
         & error stop "Error in `mpi_variance_kahan_array(n, vs)`: shape of ns(:) is illeagal."
    allocate(res%n_(s2))
    res%array_size_ = s2
    res%n_(:) = ns(:)
    allocate(res%vs_(2, s2))
    do i = 1, s2
       res%vs_(1, i) = kahan_summation(vs(1, i))
       res%vs_(2, i) = kahan_summation(vs(2, i))
    end do
  end function generate_mpi_variance_kahan_array_by_data
  pure subroutine merge_data_mpi_vka(this, other)
    class(mpi_variance_kahan_array), intent(inout) :: this
    type(mpi_variance_kahan_array), intent(in) :: other
    integer(int32) :: i
    do i = 1, this%array_size_
       this%vs_(1, i) = this%vs_(1, i) + other%vs_(1, i)%val()
       this%vs_(2, i) = this%vs_(2, i) + other%vs_(2, i)%val()
    end do
    this%n_(:) = this%n_(:) + other%n_(:)
  end subroutine merge_data_mpi_vka
  !> getter
  pure elemental integer(int64) function num_sample_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = this%n_(idx)
  end function num_sample_mpi_vka
  pure elemental real(real64) function v1_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = this%vs_(1, idx)%val()
  end function v1_mpi_vka
  pure elemental real(real64) function v1_square_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = this%vs_(2, idx)%val()
  end function v1_square_mpi_vka
  pure function get_vs_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    real(real64), allocatable :: res(:, :)
    allocate(res, source = this%vs_%val())
  end function get_vs_mpi_vka

  pure elemental real(real64) function mean_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = this%vs_(1, idx)%val() / this%n_(idx)
  end function mean_mpi_vka
  pure elemental real(real64) function square_mean_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = this%vs_(2, idx)%val() / this%n_(idx)
  end function square_mean_mpi_vka
  !> <v^2> - <v>^2
  !> == sum(v^2) / n - (sum(v) / n)^2
  !> == sum(v^2) / n - sum(v)^2 / n^2
  !> == (sum(v^2) - sum(v)^2 / n) / n
  pure elemental real(real64) function var_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = (this%vs_(2, idx)%val() - this%vs_(1, idx)%val() ** 2 / this%n_(idx)) / this%n_(idx)
  end function var_mpi_vka
  pure elemental real(real64) function sample_var_mpi_vka(this, idx) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    integer(int32), intent(in) :: idx
    res = this%var(idx) * this%n_(idx) / (this%n_(idx) - 1_int64)
  end function sample_var_mpi_vka
  impure subroutine mpi_sum_mpi_vka(this, myrank, num_proc, root, res, ierr)
    class(mpi_variance_kahan_array), intent(inout) :: this
    integer(int32), intent(in) :: myrank, num_proc, root
    type(mpi_variance_kahan_array), intent(inout) :: res
    integer(int32), intent(inout) :: ierr
    integer(int64), allocatable :: n(:)
    real(real64), allocatable :: vs(:, :)
    if (myrank == root) then
       allocate(n(this%array_size_))
    else
       allocate(n(1))
    end if
    call MPI_Reduce(this%n_(1:this%array_size_), n(1), this%array_size_, &
         & MPI_INTEGER8, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) return
    if (myrank == root) then
       allocate(vs(2, this%array_size_))
    else
       allocate(vs(2, 1))
    end if
    call MPI_Reduce(this%vs_(1:2, 1:this%array_size_)%val(), vs(1, 1), 2 * this%array_size_, &
         & MPI_REAL8, MPI_SUM, root, MPI_COMM_WORLD, ierr)
    if (ierr /= 0) return
    res = mpi_variance_kahan_array(n, vs)
  end subroutine mpi_sum_mpi_vka
end module mpi_variance_kahan_array_m
