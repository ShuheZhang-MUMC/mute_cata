function out = Gaussian_filter(in,xigma,theta,type)
% image: raw image
% xigma: radius of Gaussian filter
% theta: weight of Gaussian filter



a =  1/(2*pi*xigma^2);

[m,n,c] = size(in);
pix_x = 1:n;
pix_y = 1:m;
[pix_x,pix_y] = meshgrid(pix_x,pix_y);
r2 = (pix_x-n/2).^2 + (pix_y-m/2).^2;
GF = theta * a * exp(-r2/(2*xigma^2));
out = zeros(m,n,c);
for c_ = 1:c
if type == 'fft'
    r = in(:,:,c_);r(2*m-1,2*n-1) = 0;
    GF(2*m-1,2*n-1) = 0;
    r  = ifft2(fft2(r).*fft2(GF));
    out(:,:,c_) = real(r(fix(m/2)+1:fix(m/2)+m,fix(n/2)+1:fix(n/2)+n));
elseif type == 'cov'
    out(:,:,c_) = conv2(in(:,:,c_),GF,'same');
end
end
end