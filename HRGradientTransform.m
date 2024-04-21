function transformed = HRGradientTransform(HR, LR)
    ratio = fieldRatio(HR,LR); 
    lrgrad = imgradient(LR);
    transformed = ratio.*lrgrad;
end