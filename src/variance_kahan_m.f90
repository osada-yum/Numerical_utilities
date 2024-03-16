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
     integer(int64) :: n_ = 0_int64
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
     !> calculate
     procedure, pass :: mean => mean_mpi_vka
     procedure, pass :: square_mean => square_mean_mpi_vka
     procedure, pass :: var => var_mpi_vka
     procedure, pass :: sample_var => sample_var_mpi_vka
     !> mpi
     procedure, pass :: mpi_sum => mpi_sum_mpi_vka
  end type mpi_variance_kahan_array
  interface mpi_variance_kahan_array
     module procedure :: generate_variance_kahan
     module procedure :: generate_variance_kahan_by_data
  end interface mpi_variance_kahan_array
contains
  pure type(mpi_variance_kahan_array) function generate_variance_kahan(array_size) result(res)
    integer(int32), intent(in) :: array_size
    res%n_ = 0_int64
    res%array_size_ = array_size
    allocate(res%vs_(2, array_size), source = kahan_summation(0.0_real64))
  end function generate_variance_kahan
  pure subroutine add_data_mpi_vka(this, idx, v1)
    class(mpi_variance_kahan_array), intent(inout) :: this
    integer(int32), intent(in) :: idx
    real(real64), intent(in) :: v1
    if (idx < 1 .or. idx > this%array_size_) &
         & error stop "Error in `add_data(idx, v1)`."
    this%n_ = this%n_ + 1_int64
    this%vs_(1, idx) = this%vs_(1, idx) + v1
    this%vs_(2, idx) = this%vs_(2, idx) + v1 ** 2
  end subroutine add_data_mpi_vka
  pure type(mpi_variance_kahan_array) function generate_variance_kahan_by_data &
       & (n, vs) result(res)
    integer(int64), intent(in) :: n
    real(real64), intent(in) :: vs(:, :)
    integer(int32) :: s1, s2
    integer(int32) :: i
    s1 = size(vs, dim = 1)
    if (s1 /= 2) &
         & error stop "Error in `mpi_variance_kahan_array(n, vs)`"
    s2 = size(vs, dim = 2)
    res%n_ = n
    allocate(res%vs_(2, s2))
    do i = 1, s2
       res%vs_(1, i) = kahan_summation(vs(1, i))
       res%vs_(2, i) = kahan_summation(vs(2, i))
    end do
  end function generate_variance_kahan_by_data
  pure subroutine merge_data_mpi_vka(this, other)
    class(mpi_variance_kahan_array), intent(inout) :: this
    type(mpi_variance_kahan_array), intent(in) :: other
    this%v1_ = this%v1_ + other%v1_%val()
    this%v1_square_ = this%v1_square_ + other%v1_square_%val()
    this%n_ = this%n_ + other%n_
  end subroutine merge_data_mpi_vka
  !> getter
  pure integer(int64) function num_sample_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = this%n_
  end function num_sample_mpi_vka
  pure real(real64) function v1_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = this%v1_%val()
  end function v1_mpi_vka
  pure real(real64) function v1_square_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = this%v1_square_%val()
  end function v1_square_mpi_vka

  pure real(real64) function mean_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = this%v1_%val() / this%n_
  end function mean_mpi_vka
  pure real(real64) function square_mean_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = this%v1_square_%val() / this%n_
  end function square_mean_mpi_vka
  !> <v^2> - <v>^2
  !> == sum(v^2) / n - (sum(v) / n)^2
  !> == sum(v^2) / n - sum(v)^2 / n^2
  !> == (sum(v^2) - sum(v)^2 / n) / n
  pure real(real64) function var_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = (this%v1_square_%val() - this%v1_%val() ** 2 / this%n_) / this%n_
  end function var_mpi_vka
  pure real(real64) function sample_var_mpi_vka(this) result(res)
    class(mpi_variance_kahan_array), intent(in) :: this
    res = this%var() * this%n_ / (this%n_ - 1_int64)
  end function sample_var_mpi_vka
end module mpi_variance_kahan_array_m
