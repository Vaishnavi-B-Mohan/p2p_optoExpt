function h=etawaitbar(h,i,n)
% h=etawaitbar([h],[i],[n])
%
% Shows a waitbar including time remaining and eta for finishing
%
% Inputs:
%  h        handle for waitbar.  Obtained by initializing with h=waitbar;
%  i        index of current value in loop (1:n)
%  n        number of interations in the loop
%
% Example:
%
% n=100;
% h = etawaitbar;
% for i=1:n
%   etawaitbar(h,i,n);
%   pause(.1);
% end

% 6/22/24  gmb  wrote it.


if nargin==0  % initialize the waitbar
    tic
    h = waitbar(0);
    pos = get(h,'Position');
    set(h,'Position',[pos(1:3),70])  % increase height for 3 rows of text
else  % update the waitbar
    timeLeft = (toc/i)*(n-i);  % in secondsd
    str{1} = sprintf('%d of %d',i,n);
    str{2} = sprintf('%d mins, %d secs',floor(timeLeft/60),floor(mod(timeLeft,60)));
    str{3} = sprintf('eta: %s',datestr(now+timeLeft/(24*60*60),'HH:MM PM'));
    waitbar(i/n,h,str)
    if i==n
        close(h)  % we're done, close the window
    end
end



