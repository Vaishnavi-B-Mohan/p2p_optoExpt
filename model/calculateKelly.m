function sens = calculateKelly(p, tf, sf)
% sens = calculateKelly(p, tf, sf)
%
% 4/5/2024  gmb  removed adding NaN's to low senstivity values

% calculate Kelly values
alpha = 2*pi*sf;	%used in the formula, related to spatial frequency
v = tf ./ sf;	
k = p.param1 + p.param2 *(abs(log10( v / 3 ))) .^ 3; %p.param1 = 6.1; p.param2 = 7.3;
amax = p.param3 ./ (v + 2); % p.param3 = 45.9;
sens = k .* v .* (alpha .* alpha) .* exp(-2*alpha ./ amax);

%	Pull out the ones above 2;
% sens(sens<2) = NaN;% ones(size(sens(l)))*(0/0);  %	This sets the bad ones to NaN
sens = log10(sens/2);

% Student in Australia pointed out that to make this workf
% for standing, rather than traveling, waves we need to plot
% G/2 (see Kelly, 1979, p. 1341 and related discussion
