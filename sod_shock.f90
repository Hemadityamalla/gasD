program sodShock

        implicit none
        !Parameters and variables
        double precision, parameter :: tubeHalfLength = 5.0, Y = 1.4
        double precision, parameter :: Tfinal = 7e-3  !This needs to be made clear
        double precision :: ul, ur, pl, pr, rhol, rhor, El, Er
        double precision, allocatable, dimension(:) :: u, rho, E, p, x, unew, rhonew, Enew, pnew
        double precision, allocatable, dimension(:) :: mflux, vflux, eflux
        integer :: i, N, iter, nGhost, ncv
        double precision :: xl, xr, dx, dt, t


        dx = 0.025 !This needs to be made clear

        dt = (dx*0.1)/375.0 !This setting is based on the sound speed and a courant number of 0.2

        xl = -1.0*tubeHalfLength
        xr = tubeHalfLength

        N = int((xr - xl)/dx)
        nGhost = 2
        ncv = N + nGhost
        !Initializing the arrays and setting initial conditions

        allocate(x(N), u(ncv), rho(ncv),E(ncv), p(ncv), unew(ncv), rhonew(ncv),Enew(ncv), pnew(ncv))
        allocate(mflux(N),vflux(N),eflux(N))


        !Initial conditions
        ul = 0.0
        ur = 0.0
        rhol = 1.0
        rhor = 0.125
        pl = 1e5
        pr = 1e4
        El = pl/(rhol*(Y-1)) + 0.5*rhol*ul**2 

        Er = pr/(rhor*(Y-1)) + 0.5*rhor*ur**2
        u = 0.0
        rho = 0.0
        E = 0.0
        p = 0.0
        unew = 0.0
        rhonew = 0.0
        Enew = 0.0
        pnew = 0.0

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
        !Take care of ghost cells- we take that the tube is infinite- so we must make sure that the shock/expansion does not get too
        !close to the bondaries
        u(1) = ul
        u(ncv) = ur
        unew(1) = ul
        unew(ncv) = ur

        rho(1) = rhol
        rho(ncv) = rhor
        rhonew(1) = rhol
        rhonew(ncv) = rhor

        p(1) = pl
        p(ncv) = pr
        pnew(1) = pl
        pnew(ncv) = pr


        E(1) = El
        E(ncv) = Er
        Enew(1) = El
        Enew(ncv) = Er

        !Time loop
        t = 0.0
        iter = 1
        print *,'Starting time integration..'
        do while (t < Tfinal)
                
                call computeFluxes(dx, rho, u, p, E, mflux, vflux, eflux, N, nGhost)
                
                rhonew(2:ncv-1) = 0.5*(rho(3:ncv) + rho(1:N)) - 0.5*dt*(mflux)
                unew(2:ncv-1) = (0.5*(rho(3:ncv)*u(3:ncv) + rho(1:N)*u(1:N)) - 0.5*dt*vflux)/rhonew(2:ncv-1) !What if we divide by zero here?
                Enew(2:ncv-1) = 0.5*(E(1:ncv) + E(1:N)) - 0.5*dt*(eflux)
                pnew(2:ncv-1) = computePressure(Enew(2:ncv-1), rhonew(2:ncv-1), unew(2:ncv-1), Y, size(unew(2:ncv-1)))!rhonew(2:ncv-1)*(Y-1)*(Enew(2:ncv-1) - 0.5*rhonew(2:ncv-1)*unew(2:ncv-1)**2)

                !Updating variables
                rho = rhonew
                u = unew
                E = Enew
                p = pnew
               
                !Check for solution divergence
                if (maxval(abs(rho)) .ge. 10) stop 'Solution diverging!'
                
                !Writing data
                if (mod(iter,20) .le. 1e-15) then
                        call writeData(iter, x, u(2:ncv-1), p(2:ncv-1), rho(2:ncv-1), E(2:ncv-1), N)
                end if
                iter = iter+1
                t = t + dt
                print *,'Time= ',t

        end do
        print *,'Finished!'



        !De-allocation of variable arrays
        !deallocate(x)
        !deallocate(u, p, E, rho, unew, pnew, Enew, rhonew)
        !deallocate(mflux, vflux, eflux)
contains

        subroutine writeData(timeStep, x, u, p, rho, E, N)
                implicit none
                integer, intent(in) :: timeStep, N
                double precision, dimension(N), intent(in) :: x, u, p, rho, E
                integer :: i
                character(21) :: filename
                write(filename, '(a,i0.6,a)') 'fortran_op/', timeStep, '.txt'
                open(1, file=filename, status='new')
                write(1, *) 'x ', 'v ', 'rho ', 'p ', 'E '
                do i=1,N
                        write(1, *) x(i), u(i), rho(i), p(i), E(i)
                end do
                close(1)
        end subroutine writeData

        subroutine computeFluxes(dx, rho, u, p, E, mflux, vflux, eflux, N, nghost)
                implicit none
                integer, intent(in) :: N, nghost
                double precision, dimension(N+nghost), intent(in) :: u, p, E, rho
                double precision, intent(in) :: dx
                double precision, dimension(N), intent(out) :: mflux, vflux, eflux
                double precision, dimension(N+1) :: mf, vf, ef
                double precision, parameter :: eps = 1e-10
                integer :: ncv, i
                ncv = N + nghost


                do i=1,N+1

                        thetar = 
                end do
                mflux = (rho(3:ncv)*u(3:ncv) - rho(1:N)*u(1:N))/dx
                vflux = ((rho(3:ncv)*u(3:ncv)**2 + p(3:ncv)) - (rho(1:N)*u(1:N)**2 + p(1:N)))/dx
                eflux = ((E(3:ncv)*u(3:ncv) + p(3:ncv)*u(3:ncv)) - (E(1:N)*u(1:N) + p(1:N)*u(1:N)))/dx


        end subroutine computeFluxes

        function computeEnergy(p, Y, rho, u, arrsize) !Works only for arrays!
                implicit none
                double precision, intent(in) :: Y
                integer, intent(in) :: arrsize
                double precision, intent(in), dimension(arrsize) :: p, rho, u
                double precision, dimension(arrsize) :: computeEnergy

                computeEnergy = p/(rho*(Y-1)) + 0.5*rho*u**2
                
        end function computeEnergy

        function computePressure(e, rho, u, Y, arrsize)
                implicit none
                double precision, intent(in) :: Y
                integer, intent(in) :: arrsize
                double precision, intent(in), dimension(arrsize) :: e, rho, u
                double precision, dimension(arrsize) :: computePressure

                computePressure = rho*(Y - 1)*(e - 0.5*rho*u**2)
        end function computePressure


        function phi(r)
                implicit none
                double precision, intent(in) :: r
                double precision :: phi
                phi = 0.0 !Upwind 
                !phi = 1.0 !Lax-Wendroff
                !phi = r !Beam-Warming
                !phi = 0.5*(1.0 + r) !Fromm
                
                !phi = (r + abs(r))/(1 + abs(r)) !van Leer
                !phi = maxval((/0.0, min(1.0, 2*r), min(2.0, r)/)) !supberbee
                !phi = int(r .gt. 0.0)*min(r, 1.0) + int(r .le. 0.0)*0.0 !minmod
                !phi = max(0.0, min(1.0/3.0 + r/6.0, r)) !koren

        end function phi

end program sodShock
