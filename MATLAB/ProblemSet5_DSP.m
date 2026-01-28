% Taha Akhlaq - DSP Problem Set 5
clear; close all; clc;

%% Problem 2

% Part A
% Laplacian
h_Lap = (1/6) * [1 4 1; 4 -20 4; 1 4 1];

% Sobel Filters
hx = (1/4) * [1 0 -1; 2 0 -2; 1 0 -1];
hy = (1/4) * [-1 -2 -1; 0 0 0; 1 2 1];

% h_LapSob
h_LapSob = conv2(hx, hx) + conv2(hy, hy);

fprintf('2a: h_LapSob Coefficients\n');
disp(h_LapSob);
fprintf('size(h_LapSob) = %dx%d\n', size(h_LapSob,1), size(h_LapSob,2));
fprintf('even symmetry about center (zero-phase) = %d\n', isequal(h_LapSob, rot90(h_LapSob,2)));

% Part B

% normalized freq response
[H_Lap, fx, fy] = freqz2(h_Lap);
[H_LapSob, ~, ~] = freqz2(h_LapSob);

fprintf('\n2b\n');

imag_max_Lap = max(abs(imag(H_Lap(:))));
imag_max_LapSob = max(abs(imag(H_LapSob(:))));
fprintf('max |imag(H_Lap)| = %e\n', imag_max_Lap);
fprintf('max |imag(H_LapSob)| = %e\n', imag_max_LapSob);

H_Lap_r = real(H_Lap);
H_LapSob_r = real(H_LapSob);

max_Lap = max(H_Lap_r(:));
min_Lap = min(H_Lap_r(:));
max_LapSob = max(H_LapSob_r(:));
min_LapSob = min(H_LapSob_r(:));

fprintf('H_Lap: max = %e, min = %e\n', max_Lap, min_Lap);
fprintf('H_LapSob: max = %e, min = %e\n', max_LapSob, min_LapSob);

tol = 1e-12;
fprintf('H_Lap non-positive within tol: %d\n', max_Lap <= tol);
fprintf('H_LapSob non-positive within tol: %d\n', max_LapSob <= tol);

H_Lap = H_Lap_r;
H_LapSob = H_LapSob_r;

% Part C
% surface & contour
w_x = fx * pi;
w_y = fy * pi;

figure;
colormap(jet);

subplot(2,2,1);
surf(w_x, w_y, H_Lap); shading interp; 
title('Surface: H_{Lap}'); xlabel('wx'); ylabel('wy');

subplot(2,2,3);
contour(w_x, w_y, H_Lap, 20); grid on;
title('Contour: H_{Lap}'); xlabel('wx'); ylabel('wy');

subplot(2,2,2);
surf(w_x, w_y, H_LapSob); shading interp; 
title('Surface: H_{LapSob}'); xlabel('wx'); ylabel('wy');

subplot(2,2,4);
contour(w_x, w_y, H_LapSob, 20); grid on;
title('Contour: H_{LapSob}'); xlabel('wx'); ylabel('wy');

% Part D
% filter images & shared display scaling
data_lily = load('LilyImg.mat');
data_rodan = load('RodanImg.mat');

field_names_lily = fieldnames(data_lily);
lily = data_lily.(field_names_lily{1});

field_names_rodan = fieldnames(data_rodan);
rodan = data_rodan.(field_names_rodan{1});

lily_Lap = filter2(h_Lap, double(lily));
rodan_Lap = filter2(h_Lap, double(rodan));
lily_LapSob = filter2(h_LapSob, double(lily));
rodan_LapSob = filter2(h_LapSob, double(rodan));

mxLap = max(abs([lily_Lap(:); rodan_Lap(:)]));
mxLapSob = max(abs([lily_LapSob(:); rodan_LapSob(:)]));

figure;
colormap(gray(256));

ax = subplot(2,3,1); image(ax, lily); axis(ax,'image'); axis(ax,'off'); title(ax,'Original Lily');

ax = subplot(2,3,2);
image(ax, lily_Lap, 'CDataMapping','scaled');
clim(ax, [-mxLap mxLap]);
axis(ax,'image'); axis(ax,'off'); title(ax,'Lily Lap');

ax = subplot(2,3,3);
image(ax, lily_LapSob, 'CDataMapping','scaled');
clim(ax, [-mxLapSob mxLapSob]);
axis(ax,'image'); axis(ax,'off'); title(ax,'Lily LapSob');

ax = subplot(2,3,4); image(ax, rodan); axis(ax,'image'); axis(ax,'off'); title(ax,'Original Rodan');

ax = subplot(2,3,5);
image(ax, rodan_Lap, 'CDataMapping','scaled');
clim(ax, [-mxLap mxLap]);
axis(ax,'image'); axis(ax,'off'); title(ax,'Rodan Lap');

ax = subplot(2,3,6);
image(ax, rodan_LapSob, 'CDataMapping','scaled');
clim(ax, [-mxLapSob mxLapSob]);
axis(ax,'image'); axis(ax,'off'); title(ax,'Rodan LapSob');

%% Problem 3

% Part A (function below)

% Part B
% upsample both images by 2 (zero insertion)
lily_up = upsample_2d(lily);
rodan_up = upsample_2d(rodan);

figure;
colormap(gray(256));
ax = subplot(1,2,1); image(ax, lily_up);  axis(ax,'image'); axis(ax,'off'); title(ax,'Upsampled Lily');
ax = subplot(1,2,2); image(ax, rodan_up); axis(ax,'image'); axis(ax,'off'); title(ax,'Upsampled Rodan');

fprintf('\n3b\n');
fprintf('Lily Original (10x10):\n');  disp(lily(1:10,1:10));
fprintf('Lily Upsampled (10x10):\n'); disp(lily_up(1:10,1:10));

fprintf('Rodan Original (10x10):\n');   disp(rodan(1:10,1:10));
fprintf('Rodan Upsampled (10x10):\n');  disp(rodan_up(1:10,1:10));

% Part C
% magnitude 2-D DFT before/after upsampling (fftshift DC to center)
F_lily = fftshift(fft2(double(lily)));
F_lily_up = fftshift(fft2(double(lily_up)));
F_rodan = fftshift(fft2(double(rodan)));
F_rodan_up = fftshift(fft2(double(rodan_up)));

Mag_lily = abs(F_lily);
Mag_lily_up = abs(F_lily_up);
Mag_rodan = abs(F_rodan);
Mag_rodan_up = abs(F_rodan_up);

Disp_lily = log(Mag_lily     + 1);
Disp_lily_up = log(Mag_lily_up  + 1);
Disp_rodan = log(Mag_rodan    + 1);
Disp_rodan_up = log(Mag_rodan_up + 1);

lim_lily = [min([Disp_lily(:);  Disp_lily_up(:)]),  max([Disp_lily(:);  Disp_lily_up(:)])];
lim_rodan = [min([Disp_rodan(:); Disp_rodan_up(:)]), max([Disp_rodan(:); Disp_rodan_up(:)])];

figure;
colormap(gray(256));

ax = subplot(2,2,1);
image(ax, Disp_lily, 'CDataMapping','scaled');  clim(ax, lim_lily);
axis(ax,'image'); axis(ax,'off'); title(ax,'Lily: Spectrum (Original)');

ax = subplot(2,2,2);
image(ax, Disp_lily_up, 'CDataMapping','scaled'); clim(ax, lim_lily);
axis(ax,'image'); axis(ax,'off'); title(ax,'Lily: Spectrum (Upsampled)');

ax = subplot(2,2,3);
image(ax, Disp_rodan, 'CDataMapping','scaled'); clim(ax, lim_rodan);
axis(ax,'image'); axis(ax,'off'); title(ax,'Rodan: Spectrum (Original)');

ax = subplot(2,2,4);
image(ax, Disp_rodan_up, 'CDataMapping','scaled'); clim(ax, lim_rodan);
axis(ax,'image'); axis(ax,'off'); title(ax,'Rodan: Spectrum (Upsampled)');

%% Function
function out = upsample_2d(in)
    [R, C] = size(in);
    out = zeros(2*R, 2*C, 'like', in);
    out(1:2:end, 1:2:end) = in;
end