module kahan_summation_m
  use, intrinsic :: iso_fortran_env
  implicit none
  private
  public :: kahan_babushka_sum, kahan_babushka_klein_sum
  public :: kahan_summation
  type :: kahan_summation
     private
     real(real64) :: s_ = 0.0_real64, c_ = 0.0_real64
   contains
     procedure, pass :: val => val_kahan_summation
     procedure, pass :: add => add_kahan_summation
     procedure, pass :: subtract => subtract_kahan_summation
     procedure, pass :: negate => negate_kahan_summation
     generic :: operator(+) => add
     generic :: operator(-) => subtract, negate
  end type kahan_summation
  interface kahan_summation
     module procedure :: generate_kahan_summation
  end interface kahan_summation
  public :: operator(+), operator(-)
  interface operator(+)
     module procedure :: add_r64_ks
  end interface operator(+)
  interface operator(-)
     module procedure :: subtract_r64_ks
  end interface operator(-)
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
  pure type(kahan_summation) function add_r64_ks(val, ks) result(res)
    real(real64), intent(in) :: val
    type(kahan_summation), intent(in) :: ks
    res = ks + val
  end function add_r64_ks
  pure type(kahan_summation) function subtract_kahan_summation(this, val) result(res)
    class(kahan_summation), intent(in) :: this
    real(real64), intent(in) :: val
    res = this + (- val)
  end function subtract_kahan_summation
  pure type(kahan_summation) function subtract_r64_ks(val, ks) result(res)
    real(real64), intent(in) :: val
    type(kahan_summation), intent(in) :: ks
    res = - ks + val
  end function subtract_r64_ks
  pure type(kahan_summation) function negate_kahan_summation(this) result(res)
    class(kahan_summation), intent(in) :: this
    res%s_ = -this%s_
    res%c_ = -this%c_
  end function negate_kahan_summation

  pure real(real64) function kahan_babushka_sum(arr) result(res)
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
  end function kahan_babushka_sum
  impure real(real64) function kahan_babushka_klein_sum(arr) result(res)
    real(real64), intent(in) :: arr(:)
    real(real64) :: cs, ccs, c, cc
    integer(int32) :: i
    write(error_unit, '(a)') "The implement of `kahan_babushka_klein_sum` may be not complete."
    res = 0.0_real64
    cs = 0.0_real64
    ccs = 0.0_real64
    c = 0.0_real64
    cc = 0.0_real64
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
         associate(t => cs + c)
           if (abs(cs) >= abs(c)) then
              cc = (cs - t) + c
           else
              cc = (c - t) + cs
           end if
           cs = t
           ccs = ccs + cc
         end associate
      end do
    end associate
    res = res + c + ccs
  end function kahan_babushka_klein_sum
end module kahan_summation_m
