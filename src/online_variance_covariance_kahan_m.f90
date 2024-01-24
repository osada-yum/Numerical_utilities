module online_variance_covariance_kahan_m
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  private
  public :: online_variance_covariance_kahan
  type :: online_variance_covariance_kahan
     private
     integer(int64) :: cnts_ = 0_int64
     type(kahan_summation) :: mean1_, mean2_
     type(kahan_summation) :: m2_1_, m2_2_
     type(kahan_summation) :: c_
   contains
     procedure, pass :: add_pair => add_pair_ovc
     procedure, pass :: mean1 => mean1_ovc
     procedure, pass :: mean2 => mean2_ovc
     procedure, pass :: var1 => var1_ovc
     procedure, pass :: sample_var1 => sample_var1_ovc
     procedure, pass :: var2 => var2_ovc
     procedure, pass :: sample_var2 => sample_var2_ovc
     procedure, pass :: cov => cov_ovc
     procedure, pass :: sample_cov => sample_cov_ovc
  end type online_variance_covariance_kahan
  interface online_variance_covariance_kahan
     module procedure :: generate_online_variance_covariance
  end interface online_variance_covariance_kahan
contains
  pure type(online_variance_covariance_kahan) function generate_online_variance_covariance(n, arr1, arr2) result(res)
    integer(int32), intent(in) :: n
    real(real64), intent(in) :: arr1(n), arr2(n)
    integer(int32) :: i
    do i = 1, n
       call res%add_pair(arr1(i), arr2(i))
    end do
  end function generate_online_variance_covariance
  pure subroutine add_pair_ovc(this, v1, v2)
    class(online_variance_covariance_kahan), intent(inout) :: this
    real(real64), intent(in) :: v1, v2
    type(kahan_summation) :: d1, d2, d1_2, d2_2
    this%cnts_ = this%cnts_ + 1
    d1 = v1 - this%mean1_
    d2 = v2 - this%mean2_
    this%mean1_ = this%mean1_ + d1%val() / this%cnts_
    this%mean2_ = this%mean2_ + d2%val() / this%cnts_
    this%c_ = this%c_ + d1%val() * d2%val()
    d1_2 = v1 - this%mean1_
    d2_2 = v2 - this%mean2_
    this%m2_1_ = this%m2_1_ + d1%val() * d1_2%val()
    this%m2_2_ = this%m2_2_ + d2%val() * d2_2%val()
  end subroutine add_pair_ovc
  pure real(real64) function mean1_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%mean1_%val()
  end function mean1_ovc
  pure real(real64) function mean2_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%mean2_%val()
  end function mean2_ovc
  pure real(real64) function var1_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%m2_1_%val() / this%cnts_
  end function var1_ovc
  pure real(real64) function sample_var1_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%m2_1_%val() / (this%cnts_ - 1)
  end function sample_var1_ovc
  pure real(real64) function var2_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%m2_2_%val() / this%cnts_
  end function var2_ovc
  pure real(real64) function sample_var2_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%m2_2_%val() / (this%cnts_ - 1)
  end function sample_var2_ovc
  pure real(real64) function cov_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%c_%val() / this%cnts_
  end function cov_ovc
  pure real(real64) function sample_cov_ovc(this) result(res)
    class(online_variance_covariance_kahan), intent(in) :: this
    res = this%c_%val() / (this%cnts_ - 1)
  end function sample_cov_ovc
end module online_variance_covariance_kahan_m
