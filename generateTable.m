function [ Table ] = generateTable( num )
%GENERATETABLE Summary of this function goes here
seg = floor( num/4 ); %%% make sure the num is the multiple of 4
x=zeros(1,4*seg);
pmin = 1e-7;
p1st = 1e-6;
p2nd = 0.001;
p3rd = 0.01;
p4th = 0.1;
pmax = 1/3;
x1 = linspace(p1st, p2nd - 1e-6, seg);
x2 = linspace(p2nd, p3rd - 1e-6, seg);
x3 = linspace(p3rd, p4th, seg);
x4 = linspace(p4th, pmax, seg);
x = [pmin, x1, x2, x3, x4(2:end)];
y = -(1-2*x).*log2(1-2*x)-2*x.*log2(x); 
Table=[x;y];  

end

