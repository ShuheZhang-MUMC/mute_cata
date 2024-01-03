function [dehazed,tran_est] = fun_haze_corr_multiscale_TV_NAdam(illu_corr,I0_aim,TV_reg,layer_all,mask,name_)

[m,n] = size(illu_corr);
imwrite(illu_corr.*mask,'Y_raw.png');
Smooth_layer = zeros(m,n,layer_all);

ret = 1/2;
L = (sqrt(numel(illu_corr))/20)*ret.^(layer_all:-1:1);

for num = 1:layer_all
    out = Gaussian_filter(illu_corr,L(num),1,'fft');
%     out = RF(illu_corr,3,L(num)/2,L(num)/2,L(num)/2);
%     out = deconvtv(illu_corr,(1/(L(num))/2)^2);
    Smooth_layer(:,:,num) = out;
end

lambda_inte = 1;

begin_layer = 2;

para = 1/(layer_all-(begin_layer-1));

I0 = I0_aim*mean(illu_corr(mask));



wa = zeros(m,n,layer_all);
ap = 0.1;

cost_function = [];

step_size = 0.01;
err_bef = Inf;

momentum_wa_order1 = zeros(m,n,layer_all);beta1_wa = 0.9;
momentum_wa_order2 = zeros(m,n,layer_all);beta2_wa = 0.99;

momentum_aa_order1 = 0;beta1_aa = 0.9;
momentum_aa_order2 = 0;beta2_aa = 0.99;


for con = 1:10000
    if step_size == 0 || err_bef < 0.015
        break;
    end
    
    den = Smooth_layer(:,:,1);
    for layer = begin_layer:layer_all
        w = para./(1+exp(-wa(:,:,layer)));
        den = den + w.*Smooth_layer(:,:,layer);
    end
    T = 1 - ap*den;
    
    O = (illu_corr - 1)./(T.^2) + 1;
    
    
    
%     ant = ['iterations = ',num2str(con)];
%     
%     figure1 = figure(3);
%     imshow([imresize(O.*mask,0.6),imresize(T.*mask,0.6)],[0,1]);title(ant);
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.547915167095112 0.873724489795902 0.341544987146529 0.0778061224489772],...
%     'Color',[1 1 1],...
%     'String','Latent T_{sc}',...
%     'FontSize',24,...
%     'FitBoxToText','off');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%     [0.0696690766168355 0.874138361638351 0.341544987146529 0.0778061224489773],...
%     'Color',[1 1 1],...
%     'String','Latent Y channel',...
%     'FontSize',24,...
%     'FitBoxToText','off');
    
    
    temp_do = mean(mean(conv2(O,[-1,1],'same').^2 + conv2(O,[-1;1],'same').^2));
    cost_function(con) = mean(mean(O(mask)-I0)) - mean(mean(temp_do));
    
    err_now = (mean(mean(O(mask))) - I0)^2;
    
    temp_all = lambda_inte*2*(O(mask)-I0);
    ddO_dwa = 2*(conv2(conv2(O,[-1,1],'same'),[0,1,-1],'same') + conv2(conv2(O,[-1;1],'same'),[0;1;-1],'same'));

    % update wa
    if TV_reg ~= 0
        DTxx = conv2(T,[1,-2,1],'same');
        DTyy = conv2(T,[1;-2;1],'same');
        DTxy = conv2(T,[-1,1;1,-1],'same');

        DT2 = sqrt(DTxx.^2 + DTyy.^2 + 2 * DTxy.^2);
        DTxx = DTxx ./ (DT2 + eps);
        DTyy = DTyy ./ (DT2 + eps);
        DTxy = DTxy ./ (DT2 + eps);
    end
    for layer = begin_layer:layer_all
        temp_wa = wa(:,:,layer);
        temp_Smooth = Smooth_layer(:,:,layer);
        
        dO_dwa_fenzi = 2 * ap * para * temp_Smooth(mask) .* (illu_corr(mask) - 1);
        dO_dwa_fenmu = (1+exp(-temp_wa(mask))).^2 .* T(mask).^3;
        
        if TV_reg~=0

            temp_TV = ap * temp_Smooth * para .* exp(-temp_wa)./(1+exp(-temp_wa)).^2;
            dGdW_xx_temp = DTxx .* temp_TV;
            dGdW_yy_temp = DTyy .* temp_TV;    
            dGdW_xy_temp = DTxy .* temp_TV;  

            dGdW_xx = conv2(dGdW_xx_temp,[0,1,-2,1],'same');
            dGdW_yy = conv2(dGdW_yy_temp,[0;1;-2;1],'same');
            dGdW_xy = conv2(dGdW_xy_temp,[0,0,0;0,-1,1;0,1,-1],'same');
            
            TV_term = TV_reg * (dGdW_xx + dGdW_yy + 2 * dGdW_xy);
            
            dLdW = (temp_all - ddO_dwa(mask)) .* dO_dwa_fenzi ./ dO_dwa_fenmu + TV_term(mask);
        else
            dLdW = (temp_all - ddO_dwa(mask)) .* dO_dwa_fenzi ./ dO_dwa_fenmu;
        end

        %% NAdam
        temp = momentum_wa_order1(:,:,layer);
        temp(mask) = beta1_wa * temp(mask) + (1 - beta1_wa) * dLdW;
        momentum_wa_order1(:,:,layer) =  temp;           
        
        
        temp = momentum_wa_order2(:,:,layer);
        temp(mask) = beta2_wa * temp(mask) + (1 - beta2_wa) * dLdW.^2;
        momentum_wa_order2(:,:,layer) =  temp; 
        
        m1 = momentum_wa_order1(:,:,layer) / (1 - beta1_wa^con);
        m2 = momentum_wa_order2(:,:,layer) / (1 - beta2_wa^con);
                        
%                         
        temp_wa(mask) = temp_wa(mask) - step_size * 1./(sqrt(m2(mask)) + eps).*...
                                 (beta1_wa * m1(mask) + (1 - beta1_wa) * dLdW);
                             
        wa(:,:,layer) = temp_wa;
    end
    
    % update ap
    
    if TV_reg~=0
        dGdA_xx_temp = -den .* DTxx;
        dGdA_yy_temp = -den .* DTyy;
        dGdA_xy_temp = -den .* DTxy;
        
        dGdA_xx = conv2(dGdA_xx_temp,[0,1,-2,1],'same');
        dGdA_yy = conv2(dGdA_yy_temp,[0;1;-2;1],'same');
        dGdA_xy = conv2(dGdA_xy_temp,[0,0,0;0,-1,1;0,1,-1],'same');
        
        DLDap = mean(mean(2*(temp_all - ddO_dwa(mask)).*(illu_corr(mask)-1).*den(mask)./(T(mask).^3))) + ...
                                        TV_reg * mean(mean(dGdA_xx(mask)+dGdA_yy(mask)+dGdA_xy(mask)));
    else
        DLDap = mean(mean(2*(temp_all - ddO_dwa(mask)).*(illu_corr(mask)-1).*den(mask)./(T(mask).^3)));
    end
    
    momentum_aa_order1 = beta1_aa * momentum_aa_order1 + (1 - beta1_aa) * DLDap;
    momentum_aa_order2 = beta2_aa * momentum_aa_order2 + (1 - beta2_aa) * DLDap.^2;
    
    m1 = momentum_aa_order1 / (1 - beta1_aa^con);
    m2 = momentum_aa_order2 / (1 - beta2_aa^con);
    
    ap = ap - step_size*1./(sqrt(m2) + eps).*...
                                 (beta1_wa * m1 + (1 - beta1_wa) * DLDap);


    err_bef = err_now;
    disp(['at the ',num2str(con),'-th itter, ','the error = ',...
                    num2str(err_now),', step size = ',num2str(step_size)]); 
                

end


dehazed = O;
tran_est = T;


end