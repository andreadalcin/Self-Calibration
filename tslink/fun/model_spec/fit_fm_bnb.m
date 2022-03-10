function f = fit_fm_bnb(mss)
%FIT_FM 
 %f = fund(mss(1:3,:),mss(4:6,:));
 % f = fund_lin(mss(1:2,:),mss(4:5,:));
  f = FundamentalMatrix_Sampson_Bilinear_SeDuMi_ScaleUP(mss(1:3,:),mss(4:6,:));
  f = f(:);

end

