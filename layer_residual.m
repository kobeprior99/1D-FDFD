function res = layer_residual(p, x, f, eps0, mu0, Ltot, domainLength, bc)
    eps_list = p(1:5);
    sig_list = p(6:10);
    w        = p(11:15);

    d     = Ltot * w / sum(w);                       % thicknesses, sum = Ltot
    x0    = (domainLength - Ltot)/2;                 % center the stack
    edges = x0 + [0, cumsum(d(:)).'];

    [eps_r, sigma_e] = build_profile(x, edges, eps_list, sig_list);
    mu_r    = ones(size(x));
    sigma_m = zeros(size(x));

    G = FDFD(f, x, eps_r, mu_r, sigma_e, sigma_m, eps0, mu0, bc, bc);
    res = 1e4 * [real(G); imag(G)];                  % scaled so target ~ 0.3
end