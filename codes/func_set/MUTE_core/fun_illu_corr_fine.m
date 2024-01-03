function [illu_corr,illumination_est] = fun_illu_corr_fine(input_img,size_window,mask)
[m,n] = size(input_img);
illu_dark = get_dark_channel( 1 - input_img, size_window);

omega_a = 0.9;
omega_b = 1;
itter_n = 1; % iteration time
e = 0.0001;
while true
    F = get_F(itter_n);
    Fn = F(itter_n+2);
    if Fn < ((omega_b-omega_a)/e)
        itter_n = itter_n+1;
    else
        break
    end
end
this_itter = itter_n;
F = get_F(itter_n);

x1 = omega_a + F(itter_n-2)/F(itter_n) * (omega_b - omega_a);
x2 = omega_a + F(itter_n-1)/F(itter_n) * (omega_b - omega_a);


while (this_itter > 3)
    illu_est = 1 - x1 * illu_dark;
    x = guided_filter(1-input_img, illu_est, 15, 0.001);
    illumination_est = max(reshape(x, m, n),0.1);
    illu_corr = 1 - (((1 - input_img - 1) ./ illumination_est) + 1);
    fx1 = (sum(sum((illu_corr(mask)>1))) / sum(sum(mask)) - 0.02)^2;

    illu_est = 1 - x2 * illu_dark;
    x = guided_filter(1-input_img, illu_est, 15, 0.001);
    illumination_est = max(reshape(x, m, n),0.1);
    illu_corr = 1 - (((1 - input_img - 1) ./ illumination_est) + 1);
    fx2 = (sum(sum((illu_corr(mask)>1))) / sum(sum(mask)) - 0.02)^2;
    
    if fx2 < fx1
        omega_a = x1;
        x1 = x2;
        x2 = omega_a + F(this_itter-1)/F(this_itter)*(omega_b - omega_a);
    else
        omega_b = x2;
        x2 = x1;
        x1 = omega_a + F(this_itter-2)/F(this_itter)*(omega_b - omega_a);
    end
     this_itter = this_itter - 1;
end


illu_est = 1 - (omega_a + omega_b) / 2 * illu_dark;
x = guided_filter(1-input_img, illu_est, 15, 0.001);
illumination_est = max(reshape(x, m, n),0.1);
illu_corr = 1 - (((1 - input_img - 1) ./ illumination_est) + 1);

disp(['ill_corr = ',num2str( (omega_a + omega_b) / 2 ),...
    '; mean = ',num2str(mean(mean(illu_corr(mask)))),...
    '; avg = ',num2str(sum(sum((illu_corr>1))) / sum(sum(mask)))])


end