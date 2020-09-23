
function x = invxlnx3_fast(y,f)
%
% Fast solving y = x*log(x-2) for x, y can be a vector
%
i_large = y>1000;
i_small = y<=1000;

iyL = floor(y(i_small)/0.01)+1;
iyR = iyL + 1;
iyR(iyR>100001) = 100001;

x = zeros(size(y));
x(i_small) = f(iyL) + (y(i_small)-(iyL-1)*0.01).*(f(iyR)-f(iyL));

z = y(i_large)./log(y(i_large)-2);
for j = 1 : 20
    z = y(i_large)./log(z-2);
end
x(i_large) = z;

% z = y/log(y-2);
% for j = 1 : 10
%     z = y/log(z-2);
% end
% [z invxlnx(y)]