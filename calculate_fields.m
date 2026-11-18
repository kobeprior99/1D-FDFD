%% Computational efficient FDFD Solver 
% Note that FDFD_simplified is computationally more expensive but is more
% intuative to understand.
function [E_total, E_scat] = calculate_fields(frequency, x, eps_r, mu_r, ...
                                    sigma_e, sigma_m, eps0, mu0, ...
                                    left_boundary, right_boundary)
    % Grid Parameters
    N  = length(x);
    dx = x(2) - x(1);
    % Angular Frequency and Free-Space Parameters
    omega = 2*pi*frequency;
    k0    = omega*sqrt(eps0*mu0);
    eta0  = sqrt(mu0/eps0);
    % Incident Fields
    Einc = exp(-1j*k0*x(:));
    Hinc = -Einc/eta0;
    % Effective Complex Permittivity and Permeability
    eps_z = eps0*(eps_r(:) - 1j*sigma_e(:)/(omega*eps0));
    mu_y  = mu0*(mu_r(:)   - 1j*sigma_m(:)/(omega*mu0));
    
    %Interior n grid points
    n = (2:N-1).';

    % Tridiagonal coefficients (interior)
    a = 1 ./ (eps_z(n) .* mu_y(n-1) * omega^2 * dx^2);
    c = 1 ./ (eps_z(n) .* mu_y(n)   * omega^2 * dx^2);
    b = 1 - a - c;

    % Right-hand side (interior)
    term1 = -((eps_z(n) - eps0) ./ eps_z(n)) .* Einc(n);
    term2 = -(mu_y(n)   - mu0) ./ (1j*omega*dx*eps_z(n).*mu_y(n))   .* Hinc(n);
    term3 =  (mu_y(n-1) - mu0) ./ (1j*omega*dx*eps_z(n).*mu_y(n-1)) .* Hinc(n-1);
    
    %make Y and fill interior points
    Y = complex(zeros(N,1));
    Y(n) = term1 + term2 + term3;

    % Boundary conditions
    if strcmpi(left_boundary, 'o')
        iL = [1; 1; 1];  jL = [1; 2; 3];  vL = [-1/(2*dx); -1j*k0; 1/(2*dx)];
        Y(1) = 0;
    else
        iL = 1;  jL = 1;  vL = 1;
        Y(1) = -Einc(1);
    end

    if strcmpi(right_boundary, 'o')
        iR = [N; N; N];  jR = [N-2; N-1; N];  vR = [-1/(2*dx); 1j*k0; 1/(2*dx)];
        Y(N) = 0;
    else
        iR = N;  jR = N;  vR = 1;
        Y(N) = -Einc(N);
    end

    % Assemble
    A = sparse([n; n; n; iL; iR], ...
               [n-1; n; n+1; jL; jR], ...
               [a; b; c; vL; vR], N, N);
    %% Solve Linear System
    E_scat = A\Y;
    E_total  = Einc + E_scat;
    
end