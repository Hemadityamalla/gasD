program sodShock

        implicit none
        !Parameters and variables
        double precision, parameter :: tubeHalfLength = 1.0, Y = 1.4
        double precision, parameter :: Tfinal = 100.0  !This needs to be made clear
        double precision :: ul, ur, pl, pr, rhol, rhor, El, Er
        double precision, allocatable, dimension(:) :: u, rho, E, p, x
        integer :: i, N, iter, nGhost, ncv
        double precision :: xl, xr, dx, dt, t


        dx = 0.1 !This needs to be made clear

        dt = 0.01 !This needs to be made clear

        xl = -1.0*tubeHalfLength
        xr = tubeHalfLength

        N = int((xr - xl)/dx)
        print *, N
        nGhost = 2
        ncv = N + nGhost
        !Initializing the arrays and setting initial conditions

        allocate(x(N), u(ncv), rho(ncv),E(ncv), p(ncv))


        !Initial conditions
        ul = 0.0
        ur = 0.0
        rhol = 1.0
        rhor = 0.125
        pl = 1e5
        pr = 1e4
        El = pl/(Y-1) + 0.5*rhol*ul**2

        Er = pr/(Y-1) + 0.5*rhor*ur**2
        u = 0.0
        rho = 0.0
        E = 0.0
        p = 0.0
        do i=1,N
                x(i) = xl + (dx + (i-1)*dx)
        end do
        do i=1,N
                if (x(i) .le. 0.5*(xl + xr)) then
                        u(i+1) = ul
                        rho(i+1) = rhol
                        p(i+1) = pl
                        E(i+1) = El
                else
                        u(i+1) = ur
                        rho(i+1) = rhor
                        p(i+1) = pr
                        E(i+1) = Er
                end if
        end do
        !Take care of ghost cells


        !Testing for initial conditions
        !Space derivative discretization

        !Time integration


        !Data writing


        !De-allocation of variable arrays
        deallocate(x)
        deallocate(u, p, E, rho)
contains

        subroutine writeData(timeStep, x, u, p, rho, N)
                implicit none
                integer, intent(in) :: timeStep, N
                double precision, dimension(N), intent(in) :: x, u, p, rho
                integer :: i
                character(21) :: filename
                write(filename, '(a,i0.6,a)') 'fortran_op/', timeStep, '.txt'
                print *, 'testStuff'
                open(1, file=filename, status='new')
                write(1, *) 'x ', 'velocity ', 'rho ', 'pressure '
                do i=1,N
                        write(1, *) x(i), u(i), rho(i), p(i)
                end do
                close(1)
        end subroutine writeData
end program sodShock
