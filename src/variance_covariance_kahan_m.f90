module variance_covariance_kahan_m
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  private
  public variance_covariance_kahan
  type :: variance_covariance_kahan
     private
     integer(int64) :: n_ = 0_int64
     type(kahan_summation) :: v1_, v2_
     type(kahan_summation) :: v1_square_, v2_square_
     type(kahan_summation) :: v1v2_
   contains
     procedure, pass :: add_data => add_data_vck
     procedure, pass :: mean1 => mean1_vck
     procedure, pass :: mean2 => mean2_vck
     procedure, pass :: var1 => var1_vck
     procedure, pass :: sample_var1 => sample_var1_vck
     procedure, pass :: var2 => var2_vck
     procedure, pass :: sample_var2 => sample_var2_vck
     procedure, pass :: cov => cov_vck
     procedure, pass :: sample_cov => sample_cov_vck
  end type variance_covariance_kahan
  interface variance_covariance_kahan
     module procedure :: generate_variance_covariance_kahan
  end interface variance_covariance_kahan
contains
  pure type(variance_covariance_kahan) function generate_variance_covariance_kahan() result(res)
    res%n_ = 0_int64
    res%v1_ = kahan_summation(0.0_real64)
    res%v2_ = kahan_summation(0.0_real64)
    res%v1_square_ = kahan_summation(0.0_real64)
    res%v2_square_ = kahan_summation(0.0_real64)
    res%v1v2_ = kahan_summation(0.0_real64)
  end function generate_variance_covariance_kahan
  pure subroutine add_data_vck(this, v1, v2)
    class(variance_covariance_kahan), intent(inout) :: this
    real(real64), intent(in) :: v1, v2
    this%n_ = this%n_ + 1_int64
    this%v1_ = this%v1_ + v1
    this%v2_ = this%v2_ + v2
    this%v1_square_ = this%v1_square_ + v1 ** 2
    this%v2_square_ = this%v2_square_ + v2 ** 2
    this%v1v2_ = this%v1v2_ + v1 * v2
  end subroutine add_data_vck
  pure real(real64) function mean1_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = this%v1_%val() / this%n_
  end function mean1_vck
  pure real(real64) function mean2_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = this%v2_%val() / this%n_
  end function mean2_vck
  !> <v^2> - <v>^2
  !> == sum(v^2) / n - (sum(v) / n)^2
  !> == sum(v^2) / n - sum(v)^2 / n^2
  !> == (sum(v^2) - sum(v)^2 / n) / n
  pure real(real64) function var1_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = (this%v1_square_%val() - this%v1_%val() ** 2 / this%n_) / this%n_
  end function var1_vck
  pure real(real64) function sample_var1_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = this%var1() * this%n_ / (this%n_ - 1_int64)
  end function sample_var1_vck
  pure real(real64) function var2_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = (this%v2_square_%val() - this%v2_%val() ** 2 / this%n_) / this%n_
  end function var2_vck
  pure real(real64) function sample_var2_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = this%var2() * this%n_ / (this%n_ - 1_int64)
  end function sample_var2_vck
  pure real(real64) function cov_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = (this%v1v2_%val() - (this%v1_%val() * this%v2_%val()) / this%n_) / this%n_
  end function cov_vck
  pure real(real64) function sample_cov_vck(this) result(res)
    class(variance_covariance_kahan), intent(in) :: this
    res = this%cov() * this%n_ / (this%n_ - 1_int64)
  end function sample_cov_vck
end module variance_covariance_kahan_m
