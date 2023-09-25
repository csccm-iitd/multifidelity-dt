function ysimu = MF_predict(xsimu,dmodel1,dmodel2)

cd('./H-PCFE_p')

y_LF = predictor1(xsimu,dmodel1);

xsimu_aug = [xsimu,y_LF(:)];

ysimu = predictor1(xsimu_aug,dmodel2);

cd('../')
end



% function ysimu = MF_predict(xsimu,dmodel1,dmodel2)
% 
% cd('./H-PCFE_p')
% 
% npred = 100;
% 
% [y_LF,mse_LF] = predictor1(xsimu,dmodel1);
% 
% for ipred = 1:npred
%     ypred         = y_LF + mse_LF.*randn(size(y_LF));
%     
%     xsimu_aug = [xsimu,ypred(:)];
% 
%     ysimup(:,ipred) = predictor1(xsimu_aug,dmodel2);
% end
% ysimu       = mean(ysimup,2);
% 
% cd('../')
% end