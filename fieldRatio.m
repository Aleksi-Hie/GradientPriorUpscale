function ratio = fieldRatio(HR, LR)
    treshold = 0.1;
    length = 5;
    lambda_h = 1.6;
    lambda_l = 1.6; % Assuming the same lambda for both HR and LR
    
    HR_downsampled = imresize(HR, size(LR), 'bicubic');
    
    [hr_s, hr_curves] = gradientProfileSharpness(HR_downsampled, treshold, length);
    [lr_s, lr_curves] = gradientProfileSharpness(LR, treshold, length);
    
    hr_s = hr_s + 1e-10;
    lr_s = lr_s + 1e-10;
    
    c = (lambda_h * a(lambda_h) * lr_s * gamma(1/lambda_l)) ./ (lambda_l * a(lambda_l) * hr_s * gamma(1/lambda_h));
    x = (-((a(lambda_h) * abs(hr_curves))./hr_s).^lambda_h) + ((a(lambda_l) * abs(lr_curves)./lr_s).^lambda_l);
    exp_term = exp(x);
    ratio = c .* exp_term;
    
end

function a_val = a(lambda)
    a_val = sqrt(gamma(3/lambda) / gamma(1/lambda));
end
