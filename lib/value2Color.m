function color = value2Color(x, cm, xmin, xmax)
% function color = value2Color(x, cm, xmin, xmax)
if nargin<4
    xmin = min(x(:));
    xmax = max(x(:));
end

x(x<xmin) = xmin;
x(x>xmax) = xmax;

n_color = size(cm,1);

x_norm = (x-xmin)./(xmax - xmin);
x_color = round(x_norm*(n_color-1))+1;

color = cm(x_color,:);
