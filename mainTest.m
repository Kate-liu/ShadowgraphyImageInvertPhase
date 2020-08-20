%% main test

clear; 
clc;  
clf;
close all;


%% Laser Source（激光光源）

% Field size and sampling
% Set 10 * 10 mm field
% Sampling 4096+ 1 pixel
L0 = 5e-3;
Nx = 64;
Ny = 64;

x = L0 * linspace(-1, 1, Nx);  
y = L0 * linspace(-1, 1, Ny);
[X, Y] = meshgrid(x, y); 

% Laser Standard deviation 
% Set 1 mm
A = 1;
sigma_r = 1e-3;

% Wave length
% Red
% Variable lambda
lambda = 632e-9;

% Gaussian function with a = A, b = x-scale, c = y-scale, d = laser standard deviation
f_gauss2D = @(a, x, y, d) (a .* exp(- ((x.^2 + y.^2) / (2 * ((d)^2)))));
% E0 = f_gauss2D( A, X, Y, sigma_r );
E0 = ones(Nx, Ny);

% Figure
figure(1);
mesh(X, Y, E0);
colormap(gray);

figure(2);
mesh(X, Y, E0);

if (mod(Ny, 2)==0)
	halfNy = Ny / 2;
else
	halfNy = (Ny + 1) / 2;
end
figure(3);
plot(x, E0(halfNy, :), 'c');
axis([-5e-3 5e-3 -0.2 1.4]);



%% Test Object of the Cylinder

Delta_Phi = zeros(Nx, Ny);

r = 0.8e-3;  % Cylinder radius 0.8mm
n1 = 1;  % air Refractive
n2 = 1 - 4 * 10^(-3); % Cylinder Refractive

kAir = 2 * pi * n1 / lambda;
kPlasma = 2 * pi * n2 / lambda;

% Cylinder object influence（圆柱形被测物对光场的相位影响）
for i = 1 : Nx
    if abs(X(1, i)) >= r
        Delta_Phi(:, i) = kAir * 2 * r;

    else
        DeltaZ = sqrt(r^2 - X(1, i)^2);
        Delta_Phi(:, i) = kAir * 2 * ( r - DeltaZ) + kPlasma * 2 * DeltaZ;

    end
end

figure(11);
mesh(X, Y, Delta_Phi);
colormap(gray);
% title('被测对象对于光场的相位影响');

figure(12);
mesh(X, Y, Delta_Phi);


%% use the main_forward function
sourceMap = E0;
Phi = Delta_Phi;
verbose = 1;

% get the shadowgraph image with the corresponding deflection potential
if (verbose) disp('Obtaining the shadowgram ...'); end;

forwardTic = tic;
targetMap = main_forward(sourceMap, Phi); % usually it takes ~2s in my computer
forwardTime = toc(forwardTic);
if (verbose) disp(sprintf('Finish in %fs', forwardTime)); end


% displaying the source map, deflection potential, and the shadowgraphy image
% close all;
figure(21);
subplot(2,2,1);
imagesc(Phi); colormap default; colorbar;
title('Deflection potential');
subplot(2,2,2);
imagesc(targetMap); colormap gray;
title('Shadowgram image');
subplot(2,2,4);
plot(targetMap(ceil(end/2),:));
title('Central horizontal slice');


%% use the main_inverse_extended function

% set the inversion algorithm parameters
num_sites = floor(numel(targetMap) * 0.8); % number of sites to be tried (I recommend 0.8 * number of pixels)
% algorithm = 'quasi-newton'; % 'lbfgs' uses much less memory, 'quasi-newton' gives slightly better performance (for small input size, use quasi-newton)
algorithm = 'lbfgs';

% now invert the image
if (verbose) disp('Retrieving the deflection potential ...'); end
inverseTic = tic;
[PhiI, sites, w] = main_inverse_extended(sourceMap, targetMap, num_sites, algorithm, verbose); % it takes ~7 minutes in my computer
inverseTime = toc(inverseTic);
if (verbose) disp(sprintf('Finish in %fs', inverseTime)); end

% normalise the retrieved potential
PhiI = PhiI - min(PhiI(:));

% displaying the retrieved deflection potential
figure(31);
subplot(2,2,1);
imagesc(PhiI); colormap default; colorbar;
title('Retrieved deflection potential');
subplot(2,2,3);
plot(Phi(ceil(end/2),:), 'b-'); hold on;
plot(PhiI(ceil(end/2),:), 'g--'); hold off;
title('Horizontal slice of the deflection potentials');
subplot(2,2,2);
imagesc(targetMap); colormap gray;
title('Shadowgram image');
subplot(2,2,4);
plot(targetMap(ceil(end/2),:));
title('Central horizontal slice');



