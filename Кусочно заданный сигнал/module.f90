module chill
  
  use types, only: wp
  
  implicit none
  
contains
  
  elemental real(wp) function f(x) result(y)
    real(wp), intent(in) :: x
    
    if (x .lt. 0.4_wp) then
      y = 0.6_wp
    elseif (x .lt. 0.6_wp) then
      y = 0.2_wp
    else
      y = 0.6_wp
    end if
    
  end function f
  
  subroutine solver(n, u, c, h, dt, time)
    integer, intent(in) :: n
    real(wp), dimension(0:n+1), intent(inout) :: u
    real(wp), intent(in) :: c, h, dt, time
    real(wp), dimension(:), allocatable :: aa, bb, cc, dd
    real(wp) :: t, C1, C2, C3, C4
    integer :: i
    
    allocate(aa(n), bb(n), cc(n), dd(n))
    
    C1 = - (c + abs(c)) / (2.0_wp * h)
    C3 = (c - abs(c)) / (2.0_wp * h)
    C4 = 1.0_wp / dt
    C2 = C4 - C1 - C2
    
    aa = C1
    bb = C2
    cc = C3
    
    t = dt
    do while (t <= time)
      dd(1) = C4 * u(1) - aa(1) * u(0)
      do i = 2, n - 1
        dd(i) = C4 * u(i)
      end do
      dd(n) = C4 * u(n) - cc(n-2) * u(n+1)
      call sweep_method(aa, bb, cc, dd, u(1:n))
      call bound_condition(n, u)
      t = t + dt
    end do
    
    deallocate(aa, bb, cc, dd)
    
  end subroutine solver
  
  subroutine init_value(n, x, u)
    integer, intent(in) :: n
    real(wp), dimension(0:n+1), intent(in) :: x
    real(wp), dimension(0:n+1), intent(out) :: u
    
    u = f(x)
    
    call bound_condition(n, u)
    
  end subroutine init_value
  
  subroutine bound_condition(n, u)
    integer, intent(in) :: n
    real(wp), dimension(0:n+1), intent(inout) :: u
    
    u(0) = u(n-1)
    u(n+1) = u(2)
    
  end subroutine bound_condition
  
  subroutine sweep_method(a, b, c, d, x)
    real(wp), dimension(:), intent(in) :: a, b, c, d
    real(wp), dimension(:), intent(out) :: x
    real(wp), dimension(:), allocatable :: alpha, beta, gama
    integer :: i, n
    
    n = size(x, dim=1)
    
    allocate(alpha(n-1), beta(n), gama(n))
    
    gama(1) = b(1)
    alpha(1) = - c(1) / gama(1)
    beta(1) = d(1) / gama(1)
    
    do i = 2, n - 1
      gama(i) = b(i) + a(i) * alpha(i-1)
      alpha(i) = - c(i) / gama(i)
      beta(i) = (d(i) - a(i) * beta(i-1)) / gama(i)
    end do
    
    gama(n) = b(n) + a(n) * alpha(n-1)
    beta(n) = (d(n) - a(n) * beta(n-1)) / gama(n)
    
    x(n) = beta(n)
    
    do i = n - 1, 1, -1
      x(i) = beta(i) + alpha(i) * x(i+1)
    end do
    
    deallocate(alpha, beta, gama)
    
  end subroutine sweep_method
  
  subroutine output(x, u)
    real(wp), dimension(:), intent(in) :: x, u
    integer, parameter :: io = 102
    integer :: ios, i, n
    character(len=256) :: str
    
    n = size(x, dim=1)
    
    open(unit=io, file='result.dat', iostat=ios, iomsg=str, status="replace", action="write")
    if (ios /= 0) stop trim(str)
    do i = 1, n
      write(io, *) x(i), u(i), f(x(i))
    end do
    close(io)
    
  end subroutine output
  
  subroutine reading(n, L, c, CFL, time)
    integer, intent(out) :: n
    real(wp), intent(out) :: L, c, CFL, time
    integer, parameter :: io = 101
    integer :: ios
    character(len=256) :: str
    
    open(unit=io, file='input.dat', iostat=ios, iomsg=str, status="old", action="read")
    if (ios /= 0) stop trim(str)
    read(io, *) n
    read(io, *) L
    read(io, *) c
    read(io, *) CFL
    read(io, *) time
    close(io)
    
  end subroutine reading
  
end module chill
