function Pw=water_fill(Plinear,Le,ee)

if (Le==1)
  Pw=Plinear;
else
  %
  % Coefficients
  %
  wf_coef              = zeros(Le);
  wf_coef(:,1)         = ones(Le,1);
  wf_coef(Le,:)        = ones(1,Le);
  wf_coef(1:Le-1,2:Le) = (-1)*eye(Le-1);
  inv_wf_coef          = inv(wf_coef);
  Pw                   = zeros(Le,1);
  %
  % Independent term
  %
  wf_indep         = zeros(Le,1);
  wf_indep(1:Le-1) = (ones(Le-1,1)./ee(2:Le))-(ones(Le-1,1)./ee(1));
  wf_indep(Le)     = Plinear;
  %
  % Result
  %
  Pw = inv_wf_coef*wf_indep;
  %
  % Consistency check 1 (P_i <> NaN)
  %
  if ((~(Pw(Le,1)<=0))&(~(Pw(Le,1)>=0)))
    Pt = water_fill(Plinear,Le-1,ee(1:Le-1));
    Pw = [Pt;0];
  end;
  %
  % Consistency check 2 (P_i > 0 for all i)
  %
  if (Pw(Le,1)<0)
    Pt = water_fill(Plinear,Le-1,ee(1:Le-1));
    Pw = [Pt;0];
  end;
end;