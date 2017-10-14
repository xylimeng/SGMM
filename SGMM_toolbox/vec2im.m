% image from an univariate vector + image mask 
% mask - matrix - logical 
function y = vec2im(x, mask) 
y = NaN + mask; 
y(mask) = x;

