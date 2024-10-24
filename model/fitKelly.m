function [err,prob] = fitKelly(p,sf,tf,c,resp)
% [err,prob] = fitKelly(p,sf,tf,c,resp)
%
% calculates the -log likelihood of data in 'resp' of the TCSF model from
% DH Kelly (1979) and the weibull psychometric fuction.

% Calculate the TCSF surface
G = calculateKelly(p,sf,tf);

% Convert from log10 sensitivity to contrast thresholds
p.t = .1.^G(:);

% Calculate P(correct) for each trial
prob = weibull(p,c(:));
prob = min(prob,.99); % standard hack to avoid log(0) NaNs

% -log likelihood calculation
err = -sum(resp.*log(prob) + (1-resp).*log(1-prob));

