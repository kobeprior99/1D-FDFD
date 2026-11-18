%% Simplified FDFD Solver function for my own reference
function [reflect, transmit] = FDFD_simplified(frequency, x, eps_r, mu_r, ...
                                    sigma_e, sigma_m, eps0, mu0, ...
                                    left_boundary, right_boundary)

    %% Grid Parameters
    N  = length(x);
    dx = x(2) - x(1);

    %% Angular Frequency and Free-Space Parameters
    omega = 2*pi*frequency;
    k0    = omega*sqrt(eps0*mu0);
    eta0  = sqrt(mu0/eps0);

    %% Incident Fields
    % Unit-amplitude plane wave propagating in the +x direction
    Einc = exp(-1j*k0*x(:));

    % Incident magnetic field (assuming exp(+jwt) convention)
    Hinc = -Einc/eta0;

    %% Effective Complex Permittivity and Permeability
    eps_z = eps0*(eps_r(:) - 1j*sigma_e(:)/(omega*eps0));
    mu_y  = mu0*(mu_r(:)   - 1j*sigma_m(:)/(omega*mu0));

    %% Initialize FDFD Matrix and Excitation Vector
    A = sparse(N, N);
    Y = complex(zeros(N, 1));

    %% Construct FDFD Matrix
    for n = 2:N-1

        % Tridiagonal coefficients
        a = 1/(eps_z(n)*mu_y(n-1)*omega^2*dx^2);
        c = 1/(eps_z(n)*mu_y(n)*omega^2*dx^2);
        b = 1 - a - c;

        % f(i) knowns
        term1 = -((eps_z(n) - eps0)/eps_z(n))*Einc(n);

        term2 = -(mu_y(n) - mu0)/(1j*omega*dx*eps_z(n)*mu_y(n))*Hinc(n);

        term3 = (mu_y(n-1) - mu0)/(1j*omega*dx*eps_z(n)*mu_y(n-1))*Hinc(n-1);

        f = term1 + term2 + term3;
        % put coefficients in their appropriate positions in the matrix
        A(n,n-1) = a;
        A(n,n)   = b;
        A(n,n+1) = c;

        % construct y vector interior
        Y(n) = f;
    end

    %% Boundary Conditions

    % Left boundary
    if strcmpi(left_boundary, 'o')

        % Scattered fields only propagate outward
        A(1,1) = -1/(2*dx);
        A(1,2) = -1j*k0;
        A(1,3) = 1/(2*dx);
        Y(1)   = 0;

    else

        % Perfect electric conductor boundary
        % scattered field is the incident field completely reflected
        A(1,1) = 1;
        Y(1)   = -Einc(1);

    end

    % Right boundary
    if strcmpi(right_boundary, 'o')

        % Scattered fields only propagate outward
        A(N,N-2) = -1/(2*dx);
        A(N,N-1) = 1j*k0;
        A(N,N)   = 1/(2*dx);
        Y(N)     = 0;

    else

        % Perfect electric conductor boundary
        % scattered field is the incident field completely reflected
        A(N,N) = 1;
        Y(N)   = -Einc(N);

    end

    %% Solve Linear System
    Escat = A\Y;
    Etot  = Einc + Escat;

    %% Extract Reflection and Transmission Coefficients

    % Reflected field at the left boundary
    reflect = Escat(1)/Einc(1);

    % Total transmitted field at the right boundary
    transmit = Etot(N)/Einc(N);

end