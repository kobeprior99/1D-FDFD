function [eps_r, mu_r, sigma_e, sigma_m] = get_material_profile(x, material_boundaries, eps_r_list, mu_r_list, sigma_e_list, sigma_m_list)
    % Generates vectors for material properties for each cell lin the
    % domain based on user-defined material boundaries, and properties
    numGridPoints = length(x);
    
    % Initialize the entire domain to free space (air) by default
    eps_r   = ones(1, numGridPoints);
    mu_r    = ones(1, numGridPoints);
    sigma_e = zeros(1, numGridPoints);
    sigma_m = zeros(1, numGridPoints);
    
    numMaterials = size(material_boundaries, 1);
    
    % Loop through each layer and apply properties where x falls within [start, end]
    for m = 1:numMaterials
        x_start = material_boundaries(m, 1);
        x_end   = material_boundaries(m, 2);
        
        % Mask for the current layer interval
        material_mask = x >= x_start & x <= x_end;
        
        eps_r(material_mask)   = eps_r_list(m);
        mu_r(material_mask)    = mu_r_list(m);
        sigma_e(material_mask) = sigma_e_list(m);
        sigma_m(material_mask) = sigma_m_list(m);
    end
end