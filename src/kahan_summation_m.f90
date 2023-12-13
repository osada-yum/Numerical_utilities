module kahan_summation_m
  use, intrinsic :: iso_fortran_env
  implicit none
  type :: kahan_summation
     private
     real(real64) :: s_ = 0.0_real64, c_ = 0.0_real64
   contains
     procedure, pass :: val => val_kahan_summation
     procedure, pass :: add => add_kahan_summation
  end type kahan_summation
  interface kahan_summation
     module procedure :: generate_kahan_summation
  end interface kahan_summation
contains
  pure type(kahan_summation) function generate_kahan_summation(val) result(res)
    real(real64), intent(in) :: val
    res%s_ = val
    res%c_ = 0.0_real64
  end function generate_kahan_summation
  pure real(real64) function val_kahan_summation(this) result(res)
    class(kahan_summation), intent(in) :: this
    res = this%s_
  end function val_kahan_summation
  pure type(kahan_summation) function add_kahan_summation(this, val) result(res)
    class(kahan_summation), intent(in) :: this
    real(real64), intent(in) :: val
    real(real64) :: y, t
    y = val - this%c_
    t = this%s_ + y
    res%c_ = (t - this%s_) - y
    res%s_ = t
  end function add_kahan_summation

  pure real(real64) function kahan_sum(arr) result(res)
    real(real64), intent(in) :: arr(:)
    real(real64) :: c
    integer(int32) :: i
    res = 0.0_real64
    c = 0.0_real64
    associate(n => size(arr))
      do i = 1, n
         associate(t => res + arr(i))
           if (abs(res) >= abs(arr(i))) then
              c = c + (res - t) + arr(i)
           else
              c = c + (arr(i) - t) + res
           end if
           res = t
         end associate
      end do
    end associate
    res = res + c
  end function kahan_sum
end module kahan_summation_m
