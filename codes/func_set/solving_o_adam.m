function o = solving_o_adam(s,t,lambda,iter)


[m,n] = size(s);


Dx = psf2otf([1,-2,1],[m,n]);
Dy = psf2otf([1;-2;1],[m,n]);
Dxy = 1*psf2otf([1,-1;-1,1],[m,n]);


rho1 = 0.9;
rho2 = 0.999;
step = 0.05;

o0.para = (s - 1)./(t.^2 + eps) + 1;
o0.mom1 = 0;
o0.mom2 = 0;


clc
disp('Hessian denoising....')

err_bef = Inf;
err_now = 0;

for con = 1:iter
    
    den = t.^2.*(t.^2.*o0.para+1-t.^2-s);
   
    [L1_reg,ss] = prox_L1(o0.para);
    Grad_o = 2 * den + lambda * L1_reg;  
    
    err_now = mean(mean((t.^2.*o0.para+1-t.^2-s).^2 + ss));
    ratio = abs(err_bef - err_now)/err_bef;
    disp(['iteration:',num2str(con),', res = ',num2str(ratio)]);
    o0 = optimizer_nadam(o0,Grad_o,rho1,rho2,step,con);
    if ratio < 1e-4
        break;
    end
    err_bef = err_now;
end


o = o0.para;




end

function [out,ss] = prox_L1(o)
    
    gxx = conv2(o,[1,-2,1],'same');
    gyy = conv2(o,[1;-2;1],'same');
    gxy = conv2(o,[-1,1;1,-1],'same');
    ss = sqrt(gxx.^2 + gyy.^2 + 2*gxy.^2);   

    v1 = gxx./(ss + eps);
    v2 = gyy./(ss + eps);
    v3 = gxy./(ss + eps);

    gxx = conv2(v1,[0,1,-2,1],'same');
    gyy = conv2(v2,[0;1;-2;1],'same');
    gxy = conv2(v3,[0,0,0;0,-1,1;0,1,-1],'same');

    out = (gxx + gyy + 2 * gxy);
end

function para_new = optimizer_nadam(para_old,grad,rho1,rho2,step,iter)

para_new = para_old;
para_new.mom1 = rho1 * para_old.mom1 + (1 - rho1) * grad;
para_new.mom2 = rho2 * para_old.mom2 + (1 - rho2) * grad.^2;

m1 = para_new.mom1/(1 - rho1^iter);
m2 = para_new.mom2/(1 - rho2^iter);

para_new.para = para_old.para - step * 1./(sqrt(m2) + eps).*...
                                 (rho1 * m1 + (1 - rho1) * grad);
end