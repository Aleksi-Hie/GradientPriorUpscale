function gradient_E = GradientEnergy(Ih, Il, G, beta)
    Ih_smooth_down = imfilter(Ih, G, 'same');
    Ih_smooth_down = imresize(Ih_smooth_down, size(Il), 'bicubic');
    
    diff = Ih_smooth_down - Il;
    
    % Upsample the difference and convolve with G
    diff_up = imresize(diff, size(Ih), 'bicubic');
    diff_up_smooth = imfilter(diff_up, G, 'same');
    
    % Compute the Laplacian of Ih and ITh
    lap_Ih = del2(Ih);
    lap_ITh = imresize(del2(HRGradientTransform(Ih, Il)), size(Ih), 'bicubic');
    lap_diff = lap_Ih - lap_ITh;
    
    % Compute the gradient of the energy function
    gradient_E = diff_up_smooth - beta * lap_diff;
end
