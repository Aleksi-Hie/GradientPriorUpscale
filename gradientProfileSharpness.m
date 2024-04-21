function [profile_sharpness, curve_lengths] = gradientProfileSharpness(I, treshold, length)
    [Gx,Gy] = imgradient(I);
    magnitude = sqrt(Gx.^2 + Gy.^2);
    direction = atan2(Gy, Gx);

    edge_pixels = magnitude > treshold;
    zeta1 = 0.15; 
    zeta2 = 0.08;
    gamma_ = 5;
    neighborhood_size = 5;
    % Trace gradient profile
    profile_length = length;
    profile_sharpness = zeros(size(I));
    curve_lengths = zeros(size(I));

    for i = 1:size(I,1)
        for j = 1:size(I,2)
            if edge_pixels(i,j)
                profile = zeros(1, profile_length);
                x0 = [i,j];
                profile(1) = magnitude(i,j);
                direction_x0 = direction(i,j);
                k = 1;
                while k < profile_length
                    x = round(x0 + k * [cos(direction_x0), sin(direction_x0)]);
                    if x(1) > 0 && x(1) <= size(I,1) && x(2) > 0 && x(2) <= size(I,2)
                        new_magnitude = magnitude(x(1), x(2));
                        if new_magnitude == 0
                            break;
                        end
                        profile(k) = new_magnitude;
                        direction_x0 = direction(x(1), x(2));
                    else
                        break; % Stop if out of bounds
                    end
                    k = k + 1;
                end
                % Compute variance of profile eg the sharpness
                profile_sum = sum(profile);
                if profile_sum == 0 % Check for zero sum
                    profile_sharpness(i,j) = 0;
                else
                    m_prime = profile ./ profile_sum;
                    curve_lengths(i,j) = sqrt(sum(diff([x0; x], 1, 1).^2, 2));
                    profile_sharpness(i,j) = sqrt(sum(m_prime .* curve_lengths(i,j).^2));
                end
            end
        end
    end







% Find edge pixels
[edge_rows, edge_cols] = find(edge_pixels);
num_edge_pixels = numel(edge_rows);

edges = zeros(num_edge_pixels * (2 * neighborhood_size + 1)^2, 2);
weights = zeros(num_edge_pixels * (2 * neighborhood_size + 1)^2, 1);

idx = 0;

for k = 1:num_edge_pixels
    i = edge_rows(k);
    j = edge_cols(k);

    % Get indices of neighboring edge pixels within the neighborhood
    neighborhood_rows = max(1, min(size(I, 1), i + (-neighborhood_size:neighborhood_size)));
    neighborhood_cols = max(1, min(size(I, 2), j + (-neighborhood_size:neighborhood_size)));
    [neighborhood_i, neighborhood_j] = meshgrid(neighborhood_rows, neighborhood_cols);
    neighborhood_i = neighborhood_i(:);
    neighborhood_j = neighborhood_j(:);
    
    % Exclude the current pixel from the neighborhood
    exclude_idx = find(neighborhood_i == i & neighborhood_j == j);
    neighborhood_i(exclude_idx) = [];
    neighborhood_j(exclude_idx) = [];

    % Compute gradient differences and distances
    grad_diff = sqrt((Gx(i,j) - Gx(sub2ind(size(Gx), neighborhood_i, neighborhood_j))).^2 + ...
                    (Gy(i,j) - Gy(sub2ind(size(Gy), neighborhood_i, neighborhood_j))).^2);
    dist = sqrt((i - neighborhood_i).^2 + (j - neighborhood_j).^2);



    % Compute edge weights and store edges
    num_neighbors = numel(neighborhood_i);
    edges(idx+1:idx+num_neighbors, :) = [i * ones(num_neighbors, 1), j * ones(num_neighbors, 1)];
    weights(idx+1:idx+num_neighbors) = exp(-zeta1 * grad_diff.^2 - zeta2 * dist.^2);
    idx = idx + num_neighbors;
end

    edges = edges(1:idx, :);
    weights = weights(1:idx);
    nodes = num_edge_pixels;
    
    % Estimate sharpness
    sigma_hat = profile_sharpness(edge_pixels);
    sigma = sigma_hat;

    % Create graph
    G = graph(edges(:,1), edges(:,2), weights, nodes);
    G = simplify(G);
    A = adjacency(G, 'weighted');
    sigma = sigma(:);
    % Minimize energy
    for i = 1:20
        sigma_prev = sigma;
        sigma = (sigma_hat + gamma_ * (A * sigma)) ./ (1 + gamma_ * sum(A, 2));
        if norm(sigma - sigma_prev) < 1e-6
            break;
        end
    end
    profile_sharpness = zeros(size(I));
    profile_sharpness(edge_pixels) = sigma;
end

