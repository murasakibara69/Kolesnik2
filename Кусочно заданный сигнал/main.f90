program main
  
  use types, only: wp
  use chill, only: reading, init_value, solver, output
  
  implicit none
  
  integer :: n, i
  real(wp) :: L, c, CFL, time, h, dt
  real(wp), dimension(:), allocatable :: u, x
  
  call reading(n, L, c, CFL, time)
  
  h = L / (n - 1)
  dt = h * CFL / c
  
  allocate(x(0:n+1), u(0:n+1))
  
  x = [(h * (i - 2), i = 1, n + 2)]
  
  call init_value(n, x, u)
  
  call solver(n, u, c, h, dt, time)
  
  call output(x(1:n), u(1:n))
  
  deallocate(x, u)
  
  call system("gnuplot -p plot.plt")
  
end program main
