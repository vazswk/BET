function biao = generate_Table
pmin = 1e-6;
p1st = 1e-4;
p2nd = 0.001;
p3rd = 0.01;
p4th = 0.1;
pmax = 1/3;
x1 = linspace(pmin, p1st, 100);
x2 = linspace(p1st, p2nd, 100);
x3 = linspace(p2nd, p3rd, 300);
x4 = linspace(p3rd, p4th, 500);
x5 = linspace(p4th, pmax, 1000);
x = [x1,x2,x3,x4,x5];
y = -2*x.*log2(x) - (1-2*x).*log2(1-2*x);
biao=[x;y];   

        
