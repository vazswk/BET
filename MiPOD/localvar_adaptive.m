function v = localvar_adaptive(X,BS,KernelGroupSel)

% This function calculates the local standard deviation for matrix X.
% The standard deviation (std) is evaluated for a square region KxK
% pixels surrounding each pixel. At the boundary, the matrix is NOT
% padded. Instead, the std is calculated from available pixels only.
% If K is not an odd integer, it is floored to the closest odd integer.
%
% Input:  X   MxN matrix
%         K   size of the square region for calculating std
% Output: s   local std calculated for KxK regions
% Typical use: s=localstd(X,3);
X = double(X);

switch KernelGroupSel
    case 1
        v = Var_Kernel_N (X,BS);
    case 2
        v = Var_Kernel_I (X,BS);
    case 3
        v = Var_Kernel_C (X,BS);
    case 4
        v = Var_Kernel_D (X,BS);
    case 5
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_I (X,BS);
        v = min(v,[],3);
    case 6
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_C (X,BS);
        v = min(v,[],3);
    case 7
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_D (X,BS);
        v = min(v,[],3);
    case 8
        v(:,:,1) = Var_Kernel_I (X,BS);
        v(:,:,2) = Var_Kernel_C (X,BS);
        v = min(v,[],3);
    case 9
        v(:,:,1) = Var_Kernel_I (X,BS);
        v(:,:,2) = Var_Kernel_D (X,BS);
        v = min(v,[],3);
    case 10
        v(:,:,1) = Var_Kernel_C (X,BS);
        v(:,:,2) = Var_Kernel_D (X,BS);
        v = min(v,[],3);
    case 11
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_I (X,BS);
        v(:,:,3) = Var_Kernel_C (X,BS);
        v = min(v,[],3);
    case 12
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_I (X,BS);
        v(:,:,3) = Var_Kernel_D (X,BS);
        v = min(v,[],3);
    case 13
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_D (X,BS);
        v(:,:,3) = Var_Kernel_C (X,BS);
        v = min(v,[],3);
    case 14
        v(:,:,1) = Var_Kernel_I (X,BS);
        v(:,:,2) = Var_Kernel_D (X,BS);
        v(:,:,3) = Var_Kernel_C (X,BS);
        v = min(v,[],3);
    case 15
        v(:,:,1) = Var_Kernel_N (X,BS);
        v(:,:,2) = Var_Kernel_I (X,BS);
        v(:,:,3) = Var_Kernel_C (X,BS);
        v(:,:,4) = Var_Kernel_D (X,BS);
        v = min(v,[],3);
    case 16
        v = var_kernel(X,BS,'W');
    case 17
        v = Var_Kernel_NH (X,BS);
    case 18
        v = Var_Kernel_IH (X,BS);
    case 19
        v(:,:,1) = Var_Kernel_NH (X,BS);
        v(:,:,2) = Var_Kernel_IH (X,BS);
        v = min(v,[],3);
end

function v = var_kernel(X,KernBS,KernShape)
load('Vkernels.mat');
kern = eval(['Vkernels.BK' num2str(KernBS*11) '.' KernShape]);
x1 = conv2(X,kern,'same');	       % Local sums
x2 = conv2(X.*X,kern,'same');	   % Local quadratic sums
R = conv2(ones(size(X)),kern,'same');  % Number of matrix elements in each square region
v = x2./R-(x1./R).^2;	           % Local variance

function VN = Var_Kernel_N (X,BS)
v(:,:,1) = var_kernel(X,BS,'U');
v(:,:,2) = var_kernel(X,BS,'R');
v(:,:,3) = var_kernel(X,BS,'L');
v(:,:,4) = var_kernel(X,BS,'D');
VN = min(v,[],3);

function VNH = Var_Kernel_NH (X,BS)
v(:,:,1) = var_kernel(X,BS,'URH');
v(:,:,2) = var_kernel(X,BS,'ULH');
v(:,:,3) = var_kernel(X,BS,'RUH');
v(:,:,4) = var_kernel(X,BS,'RDH');
v(:,:,5) = var_kernel(X,BS,'LUH');
v(:,:,6) = var_kernel(X,BS,'LDH');
v(:,:,7) = var_kernel(X,BS,'DRH');
v(:,:,8) = var_kernel(X,BS,'DLH');
VNH = min(v,[],3);

function VI = Var_Kernel_I (X,BS)
v(:,:,1) = var_kernel(X,BS,'UI');
v(:,:,2) = var_kernel(X,BS,'RI');
v(:,:,3) = var_kernel(X,BS,'LI');
v(:,:,4) = var_kernel(X,BS,'DI');
VI = min(v,[],3);

function VIH = Var_Kernel_IH (X,BS)
v(:,:,1) = var_kernel(X,BS,'UIRH');
v(:,:,2) = var_kernel(X,BS,'UILH');
v(:,:,3) = var_kernel(X,BS,'RIUH');
v(:,:,4) = var_kernel(X,BS,'RIDH');
v(:,:,5) = var_kernel(X,BS,'LIUH');
v(:,:,6) = var_kernel(X,BS,'LIDH');
v(:,:,7) = var_kernel(X,BS,'DIRH');
v(:,:,8) = var_kernel(X,BS,'DILH');
VIH = min(v,[],3);

function VC = Var_Kernel_C (X,BS)
v(:,:,1) = var_kernel(X,BS,'UR');
v(:,:,2) = var_kernel(X,BS,'DR');
v(:,:,3) = var_kernel(X,BS,'DL');
v(:,:,4) = var_kernel(X,BS,'UL');
VC = min(v,[],3);

function VD = Var_Kernel_D (X,BS)
v(:,:,1) = var_kernel(X,BS,'URD');
v(:,:,2) = var_kernel(X,BS,'DRD');
v(:,:,3) = var_kernel(X,BS,'DLD');
v(:,:,4) = var_kernel(X,BS,'ULD');
VD = min(v,[],3);
