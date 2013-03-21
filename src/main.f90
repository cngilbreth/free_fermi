! free_fermi.f90: Code for calculating canonical thermodynamic quantities for
! free fermions in a harmonic trap.
! http://infty.us/free_fermi/free_fermi.html
! v1.0
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
! (1) Calculations are performed using a recursion relation method from J. Chem
! Phys *98*, 2484 (1993).
!
! (2) This code uses the MPFUN90 arbitrary-precision arithmetic library, which
! is necessary for low-temperature calculations using this method. It is
! included in the files mpfun90.f90 and mpmod90.f90. See also
! http://crd-legacy.lbl.gov/~dhbailey/mpdist/. The internal_precision parameter below
! specifies the precision used.
!
! (3) The MPFUN90 library is not covered by the above copyright. See instead the
! files mpfun90.f90 and mpmod90.f90 for copyright information about those files.

program free_fermi
  use mpmodule
  implicit none

  integer, parameter :: rk = selected_real_kind(p=15)
  integer, parameter :: internal_precision = 200
  integer, parameter :: d = 3 ! dimension
  character(len=*), parameter :: fmt = '(es15.8)'

  ! Beta: inverse temperature (1/T), with T in units of hbar * ω.
  real(rk) :: beta
  ! Energy, partition function, and heat capacities
  type(mp_real), allocatable :: E1(:),Z1(:),E2(:),Z2(:),CC1(:),CC2(:)
  ! Number of particles for first and second species
  integer :: A1, A2
  ! Misc. variables
  character(len=256) :: buf
  character(len=32)  :: obs
  real(rk) :: val

  if (command_argument_count() .lt. 4) then
     call print_help()
     stop
  end if
  call getarg(1,obs)
  if (obs == "help") then
     call print_help()
     stop
  end if
  call getarg(2,buf); read(buf,*) A1
  call getarg(3,buf); read(buf,*) A2
  call getarg(4,buf); read(buf,*) beta

  call mpinit(internal_precision)

  select case (obs)
  case ('H')
     allocate(Z1(0:A1),E1(0:A1))
     allocate(Z2(0:A2),E2(0:A2))
     call calc_Z(beta,A1,Z1)
     call calc_E(beta,A1,Z1,E1)
     call calc_Z(beta,A2,Z2)
     call calc_E(beta,A2,Z2,E2)
     val = E1(A1) + E2(A2)
  case ('H_spin')
     allocate(Z1(0:A1),E1(0:A1))
     allocate(Z2(0:A2),E2(0:A2))
     call calc_Z_spin(beta,A1,Z1)
     call calc_E_spin(beta,A1,Z1,E1)
     call calc_Z_spin(beta,A2,Z2)
     call calc_E_spin(beta,A2,Z2,E2)
     val = E1(A1) + E2(A2)
  case ('F')
     allocate(Z1(0:A1))
     allocate(Z2(0:A2))
     call calc_Z(beta,A1,Z1)
     call calc_Z(beta,A2,Z2)
     val = (-log(Z1(A1)) - log(Z2(A2)))/beta
  case ('F_spin')
     allocate(Z1(0:A1))
     allocate(Z2(0:A2))
     call calc_Z_spin(beta,A1,Z1)
     call calc_Z_spin(beta,A2,Z2)
     val = (-log(Z1(A1)) - log(Z2(A2)))/beta
  case ('C')
     allocate(CC1(0:A1))
     allocate(CC2(0:A2))
     call calc_C(beta,A1,CC1)
     call calc_C(beta,A2,CC2)
     val = CC1(A1) + CC2(A2)
  case default
     write (0,'(a)') "Error: invalid observable "//trim(obs)//"."
     write (0,'(a)') "Type ""free_fermi help"" for more info."
     stop
  end select

  write (6,fmt) val

contains


  subroutine print_help()
    implicit none
    write (6,'(a)') 'ffree: Calculate canonical-ensemble thermodynamic observables for two species'
    write (6,'(a)') '       of noninteracting fermions in a harmonic trap.'
    write (6,'(a)') ''
    write (6,'(a)') 'Usage: ffree <obs> <A1> <A2> <beta>'
    write (6,'(a)') ''
    write (6,'(a)') 'where the parameters are:'
    write (6,'(a)') ''
    write (6,'(a)') '  <obs>         Observable (or "help" for this page)'
    write (6,'(a)') '  <A1>, <A2>    Number of particles for species 1 & 2'
    write (6,'(a)') '  <beta>        Inverse temperature [units of 1/(hbar * omega)]'
    write (6,'(a)') ''
    write (6,'(a)') 'Observables: '
    write (6,'(a)') ''
    write (6,'(a)') '  H            Energy, canonical ensemble'
    write (6,'(a)') '  C            Heat capacity, canonical ensemble'
    write (6,'(a)') '  F            Free energy, canonical ensemble'
    write (6,'(a)') '  H_spin       Energy, canonical ensemble, w/ spin'
    write (6,'(a)') '  F_spin       Free energy, canonical ensemble, w/ spin'
    write (6,'(a)') ''
    write (6,'(a)') 'Notes:'
    write (6,'(a)') '  The calculations are done for a system of two species of noninteracting fermions'
    write (6,'(a)') '  moving in a trap of frequency omega.'
    write (6,'(a)') ''
    write (6,'(a)') '  Because the particles are noninteracting, the partition function factorizes,'
    write (6,'(a)') '  and the energy of the system is simply the sum of the energy for each species.'
    write (6,'(a)') ''
    write (6,'(a)') '  By default the fermions do not have a spin degree of freedom. The H_spin and F_spin'
    write (6,'(a)') '  observables, however, do include a spin degree of freedom for each species.'
    write (6,'(a)') ''
    write (6,'(a)') '  All energies and temperatures are measured in units of hbar * omega.'
  end subroutine print_help


  type(mp_real) function Sk(beta,k)
    ! Calculate the 1-particle partition function at inverse 
    ! temperature k*beta,
    !   S(k) = Σ exp(-beta k ϵ(j))
    ! Where the sum is over all s.p. states j.
    ! Input:
    !   beta:   Inverse temperature 1/T, T in units of ℏω.
    !   k:   Integer, > 0.
    ! Output:
    !   S(k), as above.
    implicit none
    type(mp_real), intent(in) :: beta
    integer,  intent(in) :: k

    Sk = (2*sinh((beta * k) / 2))**(-d)
  end function Sk


  type(mp_real) function S1k(beta,k)
    ! Calculate the derivative of S(k) w.r.t. beta,
    !   S'(k) = ∂S(k)/∂beta
    implicit none
    type(mp_real), intent(in) :: beta
    integer,  intent(in) :: k

    S1k = -k*d*(Sk(beta,k)/2) * mpcoth((beta*k)/2)
  end function S1k


  type(mp_real) function mpcoth(x)
    implicit none
    type(mp_real), intent(in) :: x
    
    mpcoth = 1/tanh(x)
  end function mpcoth


  type(mp_real) function S2k(beta,k)
    ! Calculate the second derivative of S(k) w.r.t. beta,
    !   S''(k) = ∂²S(k)/∂β².
    implicit none
    type(mp_real), intent(in) :: beta
    integer,  intent(in) :: k

    S2k = d * k**2 * (Sk(beta,k) / 4) * (d * cosh((beta*k)/2)**2 + 1)/(sinh((beta*k)/2)**2)
  end function S2k


  subroutine calc_Z(beta1,N,Z)
    ! Calculate the canonical partition functions for 0,1,...,N free fermions in
    ! a d-dimensional harmonic trap.
    implicit none
    real(rk), intent(in) :: beta1
    integer,  intent(in) :: N
    type(mp_real), intent(out) :: Z(0:N)

    integer :: k, N1
    type(mp_real) :: beta
    type(mp_real) :: Svals(1:N)

    beta = beta1

    do k=1,N
       Svals(k) = Sk(beta,k)
    end do

    Z(0) = 1
    do N1=1,N
       Z(N1) = 0
       do k=1,N1
          Z(N1) = Z(N1) + (-1)**(k+1) * Svals(k) * Z(N1-k)
       end do
       Z(N1) = Z(N1) / N1
    end do
  end subroutine calc_Z


  subroutine calc_E(beta1,N,Z,E)
    ! Calculate the canonical thermal energy E for 0,1,...,N free fermions in a
    ! d-dimensional harmonic trap.
    implicit none
    real(rk), intent(in)  :: beta1
    integer,  intent(in)  :: N
    type(mp_real), intent(in)  :: Z(0:N)
    type(mp_real), intent(out) :: E(0:N)
    type(mp_real) :: Svals(1:N), S1vals(1:N)
    integer :: k, N1
    type(mp_real) :: beta

    beta = beta1
    do k=1,N
       Svals(k) = Sk(beta,k)
       S1vals(k) = S1k(beta,k)
    end do

    E(0) = 0
    do N1=1,N
       E(N1) = 0
       do k=1,N1
          E(N1) = E(N1) + (-1)**(k+1) * (-S1vals(k) * Z(N1-K) + Svals(k) * E(N1-k))
       end do
       E(N1) = E(N1) / N1
    end do
    do N1=1,N
       E(N1) = E(N1) / Z(N1)
    end do
  end subroutine calc_E


  subroutine calc_C(beta1,N,C)
    ! Calculate the canonical ensemble heat capacity C for 0,1,...,N free
    ! fermions in a d-dimensional harmonic trap.
    implicit none
    real(rk), intent(in)  :: beta1
    integer,  intent(in)  :: N
    type(mp_real), intent(out) :: C(0:N)

    integer :: k, N1
    type(mp_real) :: Z(0:N), Z1(0:N), Z2(0:N), Svals(1:N), S1vals(1:N), S2vals(1:N), beta

    beta = beta1
    do k=1,N
       Svals(k) = Sk(beta,k)
       S1vals(k) = S1k(beta,k)
       S2vals(k) = S2k(beta,k)
    end do

    ! Z  = Z
    ! Z1 = ∂Z/∂β
    ! Z2 = ∂²Z/∂β²
    Z(0)  = 1
    Z1(0) = 0
    Z2(0) = 0
    do N1=1,N
       Z(N1)  = 0
       Z1(N1) = 0
       Z2(N1) = 0
       do k=1,N1
          Z(N1)  = Z(N1)  + (-1)**(k+1) * Svals(k) * Z(N1-k)
          Z1(N1) = Z1(N1) + (-1)**(k+1) * (S1vals(k) * Z(N1-K) + &
               Svals(k) * Z1(N1-k))
          Z2(N1) = Z2(N1) + (-1)**(k+1) * (S2vals(k) * Z(N1-k) + &
               2 * S1vals(k) * Z1(N1-k) + Svals(k) * Z2(N1-k))
       end do
       Z(N1) = Z(N1) / N1
       Z1(N1) = Z1(N1) / N1
       Z2(N1) = Z2(N1) / N1
    end do
    C(0) = 0
    do N1=1,N
       C(N1) = beta**2 * (-1 * Z1(N1)**2/Z(N1)**2 + Z2(N1)/Z(N1))
    end do
  end subroutine calc_C


!  subroutine calc_CN_dT(beta,N,C)
!    ! Calculate the canonical heat capacity C for 0,1,...,N free fermions in a
!    ! d-dimensional harmonic trap.
!    ! This version computes it via a numerical derivative. Not sure how accurate
!    ! it is. The extra precision should help.
!    implicit none
!    real(rk), intent(in)  :: beta
!    integer,  intent(in)  :: N
!    real(rk), intent(out) :: C(0:N)
!
!    real(rk), parameter :: ΔT_over_T = 0.001
!
!    real(rk) :: beta⁺, beta⁻, ΔT
!    real(rk) :: Z⁺(0:N), E⁺(0:N), Z⁻(0:N), E⁻(0:N)
!
!    ΔT = ΔT_over_T / beta
!    beta⁻ = beta/(1.d0 - ΔT * beta)  ! T - ΔT
!    beta⁺ = beta/(1.d0 + ΔT * beta)  ! T + ΔT
!
!    call calc_ZN(beta⁻,N,Z⁻)
!    call calc_EN(beta⁻,N,Z⁻,E⁻)
!    call calc_ZN(beta⁺,N,Z⁺)
!    call calc_EN(beta⁺,N,Z⁺,E⁺)
!    C = (E⁺-E⁻)/(2*ΔT)
!  end subroutine calc_CN_dT


  subroutine calc_Z_spin(beta1,N,Z)
    ! Calculate the canonical partition functions for 0,1,...,N free fermions in
    ! a d-dimensional harmonic trap, INCLUDING a spin degree of freedom.
    implicit none
    real(rk), intent(in) :: beta1
    integer,  intent(in) :: N
    type(mp_real), intent(out) :: Z(0:N)

    integer :: k, N1
    type(mp_real) :: beta
    type(mp_real) :: Svals(1:N)

    beta = beta1
    do k=1,N
       Svals(k) = Sk(beta,k)
    end do

    Z(0) = 1
    do N1=1,N
       Z(N1) = 0
       do k=1,N1
          Z(N1) = Z(N1) + (-1)**(k+1) * 2 * Svals(k) * Z(N1-k)
       end do
       Z(N1) = Z(N1) / N1
    end do
  end subroutine calc_Z_spin


  subroutine calc_E_spin(beta1,N,Z,E)
    ! Calculate the canonical thermal energy E for 0,1,...,N free fermions in a
    ! d-dimensional harmonic trap, INCLUDING a spin degree of freedom.
    implicit none
    real(rk), intent(in)  :: beta1
    integer,  intent(in)  :: N
    type(mp_real), intent(in)  :: Z(0:N)
    type(mp_real), intent(out) :: E(0:N)

    integer :: k, N1
    type(mp_real) :: beta
    type(mp_real) :: Svals(1:N), S1vals(1:N)

    beta = beta1
    do k=1,N
       Svals(k) = Sk(beta,k)
       S1vals(k) = S1k(beta,k)
    end do

    E(0) = 0
    do N1=1,N
       E(N1) = 0
       do k=1,N1
          E(N1) = E(N1) + (-1)**(k+1) * 2 * (-S1vals(k) * Z(N1-K) + Svals(k) * E(N1-k))
       end do
       E(N1) = E(N1) / N
    end do
    do N1=1,N
       E(N1) = E(N1) / Z(N1)
    end do
  end subroutine calc_E_spin


! integer function nsp(nmax)
!    ! Number of single-particles harmonic oscillator states with 2*n + l <=
!    ! Nmax.
!    implicit none
!    integer, intent(in) :: nmax
!
!    integer :: n, l
!
!    nsp = 0;
!    do n=0,nmax/2
!       do l=0,nmax-2*n
!          nsp = nsp + 2*L + 1
!       end do
!    end do
!  end function nsp

end program free_fermi
