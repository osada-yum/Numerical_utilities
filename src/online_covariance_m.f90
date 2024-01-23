module online_covariance_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: online_covariance
  type :: online_covariance
     private
     integer(int64) :: cnts_ = 0_int64
     real(real64) :: mean1_ = 0.0_real64, mean2_ = 0.0_real64
     real(real64) :: m2_1_ = 0.0_real64, m2_2_ = 0.0_real64
     real(real64) :: c_ = 0.0_real64
   contains
     procedure, pass :: add_pair => add_pair_online_variance
     procedure, pass :: mean1 => mean1_online_variance
     procedure, pass :: mean2 => mean2_online_variance
     procedure, pass :: var1 => var1_online_variance
     procedure, pass :: sample_var1 => sample_var1_online_variance
     procedure, pass :: var2 => var2_online_variance
     procedure, pass :: sample_var2 => sample_var2_online_variance
     procedure, pass :: covar => covar_online_variance
     procedure, pass :: sample_covar => sample_covar_online_variance
  end type online_covariance
  interface online_covariance
     module procedure :: generate_online_covariance
  end interface online_covariance
contains
  pure type(online_covariance) function generate_online_covariance(n, arr1, arr2) result(res)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: arr1(n), arr2(n)
    integer(int32) :: i
    do i = 1, n
       call res%add_pair(arr1(i), arr2(i))
    end do
  end function generate_online_covariance
  pure subroutine add_pair_online_variance(this, v1, v2)
    class(online_covariance), intent(inout) :: this
    real(real64), intent(in) :: v1, v2
    this%cnts_ = this%cnts_ + 1
    associate(d1 => v1 - this%mean1_, d2 => v2 - this%mean2_)
      this%mean1_ = this%mean1_ + d1 / this%cnts_
      this%mean2_ = this%mean2_ + d2 / this%cnts_
      this%c_ = this%c_ + d1 * d2
      associate(d1_2 => v1 - this%mean1_, d2_2 => v2 - this%mean2_)
        this%m2_1_ = this%m2_1_ + d1 * d1_2
        this%m2_2_ = this%m2_2_ + d2 * d2_2
      end associate
    end associate
  end subroutine add_pair_online_variance
  pure real(real64) function mean1_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%mean1_
  end function mean1_online_variance
  pure real(real64) function mean2_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%mean2_
  end function mean2_online_variance
  pure real(real64) function var1_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%m2_1_ / this%cnts_
  end function var1_online_variance
  pure real(real64) function sample_var1_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%m2_1_ / (this%cnts_ - 1)
  end function sample_var1_online_variance
  pure real(real64) function var2_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%m2_2_ / this%cnts_
  end function var2_online_variance
  pure real(real64) function sample_var2_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%m2_2_ / (this%cnts_ - 1)
  end function sample_var2_online_variance
  pure real(real64) function covar_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%c_ / this%cnts_
  end function covar_online_variance
  pure real(real64) function sample_covar_online_variance(this) result(res)
    class(online_covariance), intent(in) :: this
    res = this%c_ / (this%cnts_ - 1)
  end function sample_covar_online_variance
end module online_covariance_m
