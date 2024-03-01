module variance_kahan_m
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  private
  public variance_kahan
  type :: variance_kahan
     private
     integer(int64) :: n_ = 0_int64
     type(kahan_summation) :: v1_
     type(kahan_summation) :: v1_square_
   contains
     procedure, pass :: add_data => add_data_vk
     procedure, pass :: merge_data => merge_data_vk
     !> getter
     procedure, pass :: num_sample => num_sample_vk
     procedure, pass :: v1 => v1_vk
     procedure, pass :: v1_square => v1_square_vk
     !> calculate
     procedure, pass :: mean => mean_vk
     procedure, pass :: square_mean => square_mean_vk
     procedure, pass :: var => var_vk
     procedure, pass :: sample_var => sample_var_vk
  end type variance_kahan
  interface variance_kahan
     module procedure :: generate_variance_covariance_kahan
     module procedure :: generate_variance_covariance_kahan_by_data
  end interface variance_kahan
contains
  pure type(variance_kahan) function generate_variance_covariance_kahan() result(res)
    res%n_ = 0_int64
    res%v1_ = kahan_summation(0.0_real64)
    res%v1_square_ = kahan_summation(0.0_real64)
  end function generate_variance_covariance_kahan
  pure subroutine add_data_vk(this, v1)
    class(variance_kahan), intent(inout) :: this
    real(real64), intent(in) :: v1
    this%n_ = this%n_ + 1_int64
    this%v1_ = this%v1_ + v1
    this%v1_square_ = this%v1_square_ + v1 ** 2
  end subroutine add_data_vk
  pure type(variance_kahan) function generate_variance_covariance_kahan_by_data &
       & (n, v1, v1_square) result(res)
    integer(int64), intent(in) :: n
    real(real64), intent(in) :: v1, v1_square
    res = variance_kahan()
    res%n_ = n
    res%v1_ = kahan_summation(v1)
    res%v1_square_ = kahan_summation(v1_square)
  end function generate_variance_covariance_kahan_by_data
  pure subroutine merge_data_vk(this, other)
    class(variance_kahan), intent(inout) :: this
    type(variance_kahan), intent(in) :: other
    this%v1_ = this%v1_ + other%v1_%val()
    this%v1_square_ = this%v1_square_ + other%v1_square_%val()
    this%n_ = this%n_ + other%n_
  end subroutine merge_data_vk
  !> getter
  pure integer(int64) function num_sample_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = this%n_
  end function num_sample_vk
  pure real(real64) function v1_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = this%v1_%val()
  end function v1_vk
  pure real(real64) function v1_square_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = this%v1_square_%val()
  end function v1_square_vk

  pure real(real64) function mean_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = this%v1_%val() / this%n_
  end function mean_vk
  pure real(real64) function square_mean_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = this%v1_square_%val() / this%n_
  end function square_mean_vk
  !> <v^2> - <v>^2
  !> == sum(v^2) / n - (sum(v) / n)^2
  !> == sum(v^2) / n - sum(v)^2 / n^2
  !> == (sum(v^2) - sum(v)^2 / n) / n
  pure real(real64) function var_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = (this%v1_square_%val() - this%v1_%val() ** 2 / this%n_) / this%n_
  end function var_vk
  pure real(real64) function sample_var_vk(this) result(res)
    class(variance_kahan), intent(in) :: this
    res = this%var() * this%n_ / (this%n_ - 1_int64)
  end function sample_var_vk
end module variance_kahan_m
