module online_variance_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: online_variance
  type :: online_variance
     private
     integer(int64) :: cnts_ = 0_int64
     real(real64) :: mean_ = 0.0_real64, m2_ = 0.0_real64
   contains
     procedure, pass :: add => add_online_variance
     procedure, pass :: mean => mean_online_variance
     procedure, pass :: var => var_online_variance
     procedure, pass :: sample_var => sample_var_online_variance
  end type online_variance
  interface online_variance
     module procedure :: generate_online_variance
  end interface online_variance
contains
  pure type(online_variance) function generate_online_variance(arr) result(res)
    real(real64), intent(in) :: arr(:)
    integer(int32) :: n
    integer(int32) :: i
    n = size(arr)
    do i = 1, n
       call res%add(arr(i))
    end do
  end function generate_online_variance
  pure subroutine add_online_variance(this, new_val)
    class(online_variance), intent(inout) :: this
    real(real64), intent(in) :: new_val
    this%cnts_ = this%cnts_ + 1
    associate(delta => new_val - this%mean_)
      this%mean_ = this%mean_ + delta / this%cnts_
      associate(delta2 => new_val - this%mean_)
        this%m2_ = this%m2_ + delta * delta2
      end associate
    end associate
  end subroutine add_online_variance
  pure real(real64) function mean_online_variance(this) result(res)
    class(online_variance), intent(in) :: this
    if (this%cnts_ < 1) then
       error stop "Error in `mean_online_variance`: count of the numbers is too small."
    end if
    res = this%mean_
  end function mean_online_variance
  pure real(real64) function var_online_variance(this) result(res)
    class(online_variance), intent(in) :: this
    if (this%cnts_ < 2) then
       error stop "Error in `var_online_variance`: count of the numbers is too small."
    end if
    res = this%m2_ / this%cnts_
  end function var_online_variance
  pure real(real64) function sample_var_online_variance(this) result(res)
    class(online_variance), intent(in) :: this
    if (this%cnts_ < 2) then
       error stop "Error in `sample_var_online_variance`: count of the numbers is too small."
    end if
    res = this%m2_ / (this%cnts_ - 1)
  end function sample_var_online_variance
end module online_variance_m
