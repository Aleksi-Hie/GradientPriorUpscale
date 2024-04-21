% Take low resolution images and upscale them using gradient profile prior
origimg = imread("swan500.png");
origimg = im2gray(origimg);
origimg = double(origimg);
Il = imresize(origimg, 0.5);
rescaleFactor = 2;

% Define spatial filter G for downsampling and smoothing
sigma = [0.8, 1.2, 1.6]; % Kernel standard deviations for different upsampling factors
s = sigma(rescaleFactor-1);
filter_size = 2 * ceil(2 * s) + 1;
G = fspecial('gaussian', filter_size, s);

% Define parameters
beta = 0.5; % Weight parameter for gradient constraint
tau = 0.2; % Step size for gradient descent
num_iterations = 3; % Number of iterations for gradient descent

% Initialize HR image Ih
Ih = imresize(Il, rescaleFactor, 'bicubic'); % Initialize with upsampling
% Gradient descent optimization
regularupscale = Ih;
for iter = 1:num_iterations
    % Compute gradient of the energy function
    gradient_E = GradientEnergy(Ih, Il, G, beta);
    % Update HR image using gradient descent
    Ih = Ih - tau * gradient_E;
end

% Display the original and the new image in side-by-side subplots
figure;
subplot(1, 3, 1);
imshow(uint8(origimg));
title('Original Image');

subplot(1, 3, 2);
imshow(uint8(Ih));
title('Gradient prior');

subplot(1,3, 3);
imshow(uint8(regularupscale));
title('Bicubic interpolation')


rmseGrad = sqrt(mean((origimg(:)-Ih(:)).^2))
rmseBicubic = sqrt(mean((origimg(:)-regularupscale(:)).^2))