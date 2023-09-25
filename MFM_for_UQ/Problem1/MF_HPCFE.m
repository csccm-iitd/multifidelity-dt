% main code

function [dmodel_LF, dmodel_HF] = MF_HPCFE(xc,yc,xf,yf)

nvar                    = size(xc,2);
lob                     = 0.1*ones(nvar,1);
upb                     = 20*ones(nvar,1);
theta0                  = 10*ones(nvar,1);

cd('./H-PCFE_p/')

dmodel_LF               = PK_fit(xc, yc, @PCFE, @corrgauss, theta0, lob, upb);

yf_pred                 = predictor1(xf,dmodel_LF);

xf_aug                  = [xf, yf_pred(:)];

lob                     = 0.1*ones(nvar+1,1);
upb                     = 20*ones(nvar+1,1);
theta0                  = 10*ones(nvar+1,1);

dmodel_HF               = PK_fit(xf_aug, yf, @PCFE, @corrgauss, theta0, lob, upb);

cd('../')

end