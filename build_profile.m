function [eps_r, sigma_e] = build_profile(x, edges, eps_list, sig_list)
    dx = x(2) - x(1);
    x  = x(:);
    eps_r   = ones(size(x));
    sigma_e = zeros(size(x));
    for i = 1:numel(eps_list)
        a = edges(i);  b = edges(i+1);
        cov = max(0, min(b, x + dx/2) - max(a, x - dx/2)) / dx;  % 0..1 coverage
        eps_r   = eps_r   + cov*(eps_list(i) - 1);
        sigma_e = sigma_e + cov*sig_list(i);
    end
end