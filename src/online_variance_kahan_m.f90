module online_variance_kahan_m
  use, intrinsic :: iso_fortran_env
  use kahan_summation_m
  implicit none
  private
  public :: online_variance_kahan
  type :: online_variance_kahan
     private
     integer(int64) :: cnts_ = 0_int64
     type(kahan_summation) :: mean_, m2_
   contains
     procedure, pass :: add => add_online_variance_kahan
     procedure, pass :: mean => mean_online_variance_kahan
     procedure, pass :: var => var_online_variance_kahan
     procedure, pass :: sample_var => sample_var_online_variance_kahan
  end type online_variance_kahan
  interface online_variance_kahan
     module procedure :: generate_online_variance_kahan
  end interface online_variance_kahan
contains
  pure type(online_variance_kahan) function generate_online_variance_kahan(arr) result(res)
    real(real64), intent(in) :: arr(:)
    integer(int32) :: n
    integer(int32) :: i
    n = size(arr)
    do i = 1, n
       call res%add(arr(i))
    end do
  end function generate_online_variance_kahan
  pure subroutine add_online_variance_kahan(this, new_val)
    class(online_variance_kahan), intent(inout) :: this
    real(real64), intent(in) :: new_val
    type(kahan_summation) :: delta, delta2
    this%cnts_ = this%cnts_ + 1
    delta = new_val - this%mean_
    this%mean_ = this%mean_ + delta%val() / this%cnts_
    delta2 = new_val - this%mean_
    this%m2_ = this%m2_ + delta%val() * delta2%val()
  end subroutine add_online_variance_kahan
  pure real(real64) function mean_online_variance_kahan(this) result(res)
    class(online_variance_kahan), intent(in) :: this
    if (this%cnts_ < 1) then
       error stop "Error in `mean_online_variance_kahan`: count of the numbers is too small."
    end if
    res = this%mean_%val()
  end function mean_online_variance_kahan
  pure real(real64) function var_online_variance_kahan(this) result(res)
    class(online_variance_kahan), intent(in) :: this
    if (this%cnts_ < 2) then
       error stop "Error in `var_online_variance_kahan`: count of the numbers is too small."
    end if
    res = this%m2_%val() / this%cnts_
  end function var_online_variance_kahan
  pure real(real64) function sample_var_online_variance_kahan(this) result(res)
    class(online_variance_kahan), intent(in) :: this
    if (this%cnts_ < 2) then
       error stop "Error in `sample_var_online_variance_kahan`: count of the numbers is too small."
    end if
    res = this%m2_%val() / (this%cnts_ - 1)
  end function sample_var_online_variance_kahan
end module online_variance_kahan_m
