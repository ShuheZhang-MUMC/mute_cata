function out = L0smoothing(img,lambda)
[m,n,c] = size(img);
out = zeros(m,n,c);

beta0 = 2*lambda;
beta_max = 10^5;

fil_dx = [1,-1];
fil_dy = [1;-1];

    dx = [-1,1];
    dy = [-1;1];

tic
for con = 1:c
    blur_B_w = img(:,:,con);
    O = blur_B_w;
    beta_ = beta0;
    ft_S = fft2(O);

    oft_dx = psf2otf(dx,size(O));
    oft_dy = psf2otf(dy,size(O));
    denGrad = abs(oft_dx).^2 + abs(oft_dy).^2;
    
    while beta_ < beta_max

       %% solving G subproblem
        for nnn = 1:3
        Gx = imfilter(O,fil_dx,'replicate');
        Gy = imfilter(O,fil_dy,'replicate');
        t = (Gx.^2 + Gy.^2) < lambda/beta_;
        Gx(t) = 0;
        Gy(t) = 0;
        
        %% solving O subproblem
        
        Fs = (ft_S + beta_*(fft2(Gx).*conj(oft_dx)...
                           +fft2(Gy).*conj(oft_dy)))./(1 + beta_*denGrad);
        
        O = real(ifft2(Fs)); 
        
        end
        beta_ = beta_ * 2;
    end
    out(:,:,con) =  O;
end

end