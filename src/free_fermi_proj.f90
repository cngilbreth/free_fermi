! free_fermi_proj.f90: Code for calculating canonical thermodynamic quantities
! for free fermions in a harmonic trap of arbitrary dimension.
! This version uses particle number projection, and is therefore most suitable for
! low temperatures.
! http://infty.us/free_fermi/free_fermi.html
! v1.0, March 2013
!
! Copyright (c) 2013 Christopher N. Gilbreth
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! Notes:
!
! This version originally written just as a check against free_fermi.f90.
program free_fermi_proj
  implicit none

  integer, parameter :: rk = kind(1d0)
  integer, parameter :: d = 3 ! dimension
  character(len=*), parameter :: fmt = '(es15.8)'
  real(rk), parameter :: pi=3.141592653589793d0

  ! Beta: inverse temperature (1/T), with T in units of hbar * ω.
  real(rk) :: beta
  ! Number of particles for first and second species
  integer :: A
  ! Misc. variables
  character(len=256) :: buf
  character(len=32)  :: obs
  real(rk) :: val, Z

  if (command_argument_count() .lt. 3) then
     call print_help()
     stop
  end if
  call getarg(1,obs)
  if (obs == "help") then
     call print_help()
     stop
  end if
  call getarg(2,buf); read(buf,*) A
  call getarg(3,buf); read(buf,*) beta

  select case (obs)
  case ('mu')
     call find_mu(beta,A,val)
  case ('Z')
     call calc_Z(beta,A,val)
  case ('F')
     call calc_Z(beta,A,val)
     val = -log(val)/beta
  case ('E')
     call calc_Z(beta,A,Z)
     call calc_E(beta,A,Z,val)
  case default
     write (0,'(a)') "Error: invalid observable "//trim(obs)//"."
     write (0,'(a)') "Type ""free_fermi help"" for more info."
     stop
  end select

  write (6,fmt) val

contains


  subroutine print_help()
    implicit none
    write (6,'(a)') 'ffree: Calculate canonical-ensemble thermodynamic observables for'
    write (6,'(a)') '       noninteracting fermions in a harmonic trap.'
    write (6,'(a)') ''
  end subroutine print_help


  subroutine find_mu(beta,N,mu)
    implicit none
    real(rk),   intent(in)  :: beta
    integer,    intent(in)  :: N
    real(rk),   intent(out) :: mu

    ! the function may not have an exact zero in floating-point arithmetic, so
    ! a y_accuracy parameter is useful
    real(8), parameter :: y_accuracy = 1d-10
    real(8), parameter :: x_accuracy = 1d-10

    real(8) :: xlb, xub, fmid, xmid, flb, fub
    real*8  :: s
    integer :: i, max_iterations

    ! Bounds: μ ∈ [xlb,xub]
    xlb = -200.d0
    xub = 200.d0

    ! "exponent" gives logarithm base 2
    max_iterations = exponent((xub - xlb)/x_accuracy) + 1
    max_iterations = max_iterations * 2

    flb = calc_Nmu(beta,xlb) - N
    fub = calc_Nmu(beta,xub) - N
    if (.not. flb*fub < 0.d0) stop "find_mu: Zero not properly bracketed"

    ! Flip sign of function if necessary to make it increasing
    s = merge(1.d0,-1.d0,fub > 0.d0)

    do i=1,max_iterations
       ! Try the midpoint in [xlb,xub]
       xmid = (xlb + xub) * 0.5d0
       fmid = s*(calc_Nmu(beta,xmid) - N)
       ! Check for solution
       if (abs(fmid) <= y_accuracy .and. &
            abs(xub - xlb) < x_accuracy) then
          mu = xmid
          goto 10
       end if
       ! Decrease interval
       if (fmid <= 0.d0) then
          xlb = xmid
       else
          xub = xmid
       end if
    end do
    stop "Error in find_mu: too many iterations."
10  return
    ! Postconditions:
    !   1. abs(fmid) <= y_accuracy
    !   2. xmid = (xlb + xub)/2
    !   2. abs(xub - xlb) < x_accuracy
    !   3. sign(f(xlb)) .ne. sign(f(xub))
  end subroutine find_mu
       

  function calc_Nmu(beta,mu) result(nmu)
    implicit none
    real(rk), intent(in) :: beta,mu
    real(rk) :: nmu

    real(rk) :: p, term
    integer :: k

    nmu = 0.d0
    k = 0
    do
       p = exp(beta * (k + 1.5d0 - mu))
       term = 1.d0/(1.d0 + p) * (((k+1)*(k+2))/2)
       if (abs(term) .lt. epsilon(1d0)) exit
       nmu = nmu + term
       k = k + 1
    end do
  end function calc_Nmu


  subroutine calc_Z(beta,A,Z)
    ! Calculate the canonical partition functions for N free fermions in
    ! a d-dimensional harmonic trap.
    ! Notes:
    !   Z = Tr_N exp(-β h)
    !     = [1/2π] * ∫dφ exp(-i φ N) det(1 + exp(-β h) exp(i φ))
    implicit none
    real(rk), intent(in) :: beta
    integer,  intent(in) :: A
    real(rk), intent(out) :: Z

    real(rk), parameter :: lb = 0.d0, ub = 2.d0*pi
    integer,  parameter :: nreg = 128

    complex(rk) :: sum
    real(rk) :: phi, dphi, mu
    integer :: i

    call find_mu(beta,A,mu)

    ! Numerical integration -- trapezoidal
    dphi = (ub - lb)/nreg
    phi = lb
    sum = trgc_phi(beta,mu,A,phi) * dphi/2
    do i=1,nreg-1
       phi = phi + dphi
       sum = sum + trgc_phi(beta,mu,A,phi)*dphi
    end do
    sum = sum + trgc_phi(beta,mu,A,ub) * dphi/2

    Z = real(sum) / (2*pi)     
  end subroutine calc_Z



  function trgc_phi(beta,mu,A,phi) result(tr)
    ! Input:
    !   beta:  Inverse temperature
    !   mu:    Chemical potential
    !   A:     Number of particles
    ! Output:
    !   Return value:
    !     trgc = exp(-iφA - βμA) det(1 + exp(-βh) exp(iφ + βμ))
    implicit none
    real(rk), intent(in) :: beta,mu,phi
    integer, intent(in) :: A
    complex(rk) :: tr

    complex(rk) :: r, term
    integer :: m

    term = exp(-beta * (1.5d0 - mu)) * exp((0.d0,1.d0) * phi)
    r = exp(-beta)
    tr = 1.d0
    m = 0
    do
       if (abs(term) .lt. epsilon(1d0)) exit
       tr = tr * (1.d0 + term)**(((m+1)*(m+2))/2)
       term = term * r
       m = m + 1
    end do
    tr = exp(-(0.d0,1.d0) * phi * A - beta * mu * A) * tr
  end function trgc_phi


  subroutine calc_E(beta,A,Z,E)
    ! Calculate the canonical energy for N free fermions in a d-dimensional
    ! harmonic trap.
    ! Notes:
    !   E = Tr_N [exp(-β h) h]/Z
    !     = 1/(2π Z) * ∫dφ exp(-iφA - βμ) Tr[exp(-βh) h exp(iφ + βμ)]
    implicit none
    real(rk), intent(in) :: beta
    integer,  intent(in) :: A
    real(rk), intent(in) :: Z
    real(rk), intent(out) :: E

    real(rk), parameter :: lb = 0.d0, ub = 2.d0*pi
    integer,  parameter :: nreg = 128

    complex(rk) :: sum
    real(rk) :: phi, dphi, mu
    integer :: i

    call find_mu(beta,A,mu)

    ! Numerical integration -- trapezoidal
    dphi = (ub - lb)/nreg
    phi = lb
    sum = trh_phi(beta,mu,A,phi) * dphi/2
    do i=1,nreg-1
       phi = phi + dphi
       sum = sum + trh_phi(beta,mu,A,phi)*dphi
    end do
    sum = sum + trh_phi(beta,mu,A,ub) * dphi/2

    E = real(sum) / (2*pi) / Z     
  end subroutine calc_E



  function trh_phi(beta,mu,A,phi) result(tr)
    ! Compute
    !   <h>(φ) = exp(-iφA - βμA) Tr(exp[-β(h - μ) + iφ] h)
    !          = exp(-iφA - βμA) 
    !              * {Tr(exp[-β(h - μ) + iφ] h) / Tr(exp[-β(h - μ) + iφ])} 
    !              * Tr(exp[-β(h - μ) + iφ]
    !          = exp(-iφA - βμA) 
    !              * {∑_k g_k ϵ_k/(1 + exp[β(h - μ) - iφ])} 
    !              * Tr(exp[-β(h - μ) + iφ].
    implicit none
    real(rk), intent(in) :: beta,mu,phi
    integer,  intent(in) :: A

    integer :: k
    complex(rk) :: tr, fac, term
    real(rk) :: r

    fac = exp(beta * (1.5d0 - mu) - (0.d0,1.d0)*phi)
    k = 0; tr = 0; r = exp(beta)
    do
       term = 1.d0/(1.d0 + fac) * (((k+1)*(k+2))/2) * (k + 1.5d0)
       if (abs(term) .lt. epsilon(1d0)) exit
       tr = tr + term
       fac = fac * r
       k = k + 1
    end do
    tr = tr * trgc_phi(beta,mu,A,phi)
  end function trh_phi

end program free_fermi_proj
