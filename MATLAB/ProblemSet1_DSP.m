% Taha Akhlaq - DSP Problem Set 1
clear; close all; clc;

%% Problem 2:

% paramaters
num_w_pts = 1000;
omega_grid = linspace(0, pi, num_w_pts);
db_floor_val = 1e-12;

% Rectangular
len_rect = 30;
rect_win_ts = ones(len_rect,1);
rect_win_ts = rect_win_ts / sum(rect_win_ts); % sum = 1 -> W(0) = 1 (0 dB)
if abs(sum(rect_win_ts)-1) > 1e-12
    rect_win_ts = rect_win_ts / sum(rect_win_ts);
end

[H_rect_z, ~] = freqz(rect_win_ts, 1, omega_grid);
mag_rect_db_vec = 20*log10(abs(H_rect_z)+db_floor_val);

% first null after dc: detect first local minimum of |H|
mag_rect_abs_vec = abs(H_rect_z);
idx_min1 = find(mag_rect_abs_vec(2:end-1) <= mag_rect_abs_vec(1:end-2) & ...
                mag_rect_abs_vec(2:end-1) <= mag_rect_abs_vec(3:end), 1, 'first');
if isempty(idx_min1)
    omega_null_rect = NaN;
    mlw_rect_val = NaN;
else
    idx_null_rect1 = idx_min1 + 1;
    omega_null_rect = omega_grid(idx_null_rect1);
    mlw_rect_val = 2 * omega_null_rect; % null-to-null mainlobe width
end

% Chebyshev
psl_db_target = 30;
best_err_cheb = inf;
len_cheb_sel = NaN;
cheb_win_ts = [];
H_cheb_z = [];
mag_cheb_db_vec = [];
mlw_cheb_val = NaN;

for nn = 26:40 % trial N near 30 to match rectangular mainlobe width
    tmp_win = chebwin(nn, psl_db_target);
    tmp_win = tmp_win / sum(tmp_win); % sum=1
    if abs(sum(tmp_win)-1) > 1e-12
        tmp_win = tmp_win / sum(tmp_win);
    end

    H_tmp = freqz(tmp_win, 1, omega_grid);
    mag_tmp_abs_vec = abs(H_tmp);
    idx_min2 = find(mag_tmp_abs_vec(2:end-1) <= mag_tmp_abs_vec(1:end-2) & ...
                    mag_tmp_abs_vec(2:end-1) <= mag_tmp_abs_vec(3:end), 1, 'first');
    if isempty(idx_min2)
        continue;
    end

    omega_null_tmp = omega_grid(idx_min2+1);
    mlw_tmp_val = 2 * omega_null_tmp;
    err_tmp = abs(mlw_tmp_val - mlw_rect_val);

    if err_tmp < best_err_cheb
        best_err_cheb = err_tmp;
        len_cheb_sel = nn;
        cheb_win_ts = tmp_win;
        H_cheb_z = H_tmp;
        mag_cheb_db_vec = 20*log10(abs(H_tmp)+db_floor_val);
        mlw_cheb_val = mlw_tmp_val;
    end
end

% Kaiser
beta_grid_vec = 0:0.5:12; % trial beta values
best_err_kais = inf;
beta_kais_sel = NaN;
kais_win_ts = [];
H_kais_z = [];
mag_kais_db_vec = [];
mlw_kais_val = NaN;

for beta_try = beta_grid_vec
    tmp_win = kaiser(len_cheb_sel, beta_try);
    tmp_win = tmp_win / sum(tmp_win); % sum = 1
    if abs(sum(tmp_win)-1) > 1e-12
        tmp_win = tmp_win / sum(tmp_win);
    end

    H_tmp = freqz(tmp_win, 1, omega_grid);
    mag_tmp_abs_vec = abs(H_tmp);
    idx_min3 = find(mag_tmp_abs_vec(2:end-1) <= mag_tmp_abs_vec(1:end-2) & ...
                    mag_tmp_abs_vec(2:end-1) <= mag_tmp_abs_vec(3:end), 1, 'first');
    if isempty(idx_min3)
        continue;
    end

    omega_null_tmp = omega_grid(idx_min3+1);
    mlw_tmp_val = 2 * omega_null_tmp;
    err_tmp = abs(mlw_tmp_val - mlw_cheb_val); % match chebyshev (and rectangular) width

    if err_tmp < best_err_kais
        best_err_kais = err_tmp;
        beta_kais_sel = beta_try;
        kais_win_ts = tmp_win;
        H_kais_z = H_tmp;
        mag_kais_db_vec = 20*log10(abs(H_tmp)+db_floor_val);
        mlw_kais_val = mlw_tmp_val;
    end
end

% plot
figure('Color',[1 1 1]); hold on; grid on; box on;
plot(omega_grid, mag_rect_db_vec, 'LineWidth', 1.25);
plot(omega_grid, mag_cheb_db_vec, 'LineWidth', 1.25);
plot(omega_grid, mag_kais_db_vec, 'LineWidth', 1.25);
xlim([0 pi]); ylim([-50 0]);
xlabel('\omega (rad/sample)', 'Interpreter','tex');
ylabel('Magnitude (dB)', 'Interpreter','tex');
title('Problem 2: |W(e^{j\omega})| - Rect vs Chebwin vs Kaiser (sum(w)=1)', 'Interpreter','tex');
legend('Rect, N=30', ...
       sprintf('Chebwin, N=%d, PSL=%d dB', len_cheb_sel, psl_db_target), ...
       sprintf('Kaiser, N=%d, \\beta=%.1f', len_cheb_sel, beta_kais_sel), ...
       'Location','southwest');

%% Problem 3:

% parameters
fs_hz = 10e3;
M_samp = 1000;
N_dft = 1024;
a = [1.0 0.8 0.9];
phi = [pi/2 pi/3 -2*pi/5];
f_hz_tones = [1e3 2e3 3.5e3];

% signal
t = (0:M_samp-1)/fs_hz;
x = zeros(1,M_samp);
for m = 1:3
    x = x + a(m)*cos(2*pi*f_hz_tones(m)*t + phi(m)); % sum of cosines
end

% window + dft
w_ham = hamming(M_samp);
x_win = x(:) .* w_ham;
X_dft = fft(x_win, N_dft);
X_dft_shift = fftshift(X_dft); % centered spectrum

% bin-aligned frequency vectors
f_hz = (0:N_dft-1) * (fs_hz/N_dft); % unshifted
f_hz_shift = (-N_dft/2:N_dft/2-1) * (fs_hz/N_dft); % shifted

% magnitude
mag_db = 20*log10(abs(X_dft) + db_floor_val); % unshifted
mag_db_shift = 20*log10(abs(X_dft_shift) + db_floor_val); % shifted

% plots
% full spectrum
figure('Color',[1 1 1]); hold on; grid on; box on;
plot(f_hz_shift, mag_db_shift, 'LineWidth', 1.6);
xlim([-fs_hz/2 fs_hz/2]);
ylim([-110 50]);
xlabel('f (Hz)');
ylabel('Magnitude (dB)');
title('Problem 3: N=1024 DFT of Hamming-Windowed x[n] (Full, fftshift)');

% nonnegative
figure('Color',[1 1 1]); hold on; grid on; box on;
idx_pos = 1:(N_dft/2+1);
plot(f_hz(idx_pos), mag_db(idx_pos), 'LineWidth', 1.6);
xlim([0 fs_hz/2]);
ylim([-110 50]);
xlabel('f (Hz)');
ylabel('Magnitude (dB)');
title('Problem 3: N=1024 DFT of Hamming-Windowed x[n] (Nonnegative)');


% peak indices
k_list = zeros(1, 2*numel(f_hz_tones));
for m = 1:numel(f_hz_tones)
    k_plus = round(f_hz_tones(m) * N_dft / fs_hz); % +freq peak
    k_minus = mod(N_dft - k_plus, N_dft); % -freq peak
    k_list(2*m-1) = k_plus;
    k_list(2*m) = k_minus;
end
k_list = sort(k_list);
fprintf('Problem 3: peak bin indices k (0..N-1): ');
fprintf('%d ', k_list);
fprintf('\n');

%% Problem 4:

% parameters
fs_hz = 44.1e3; % sampling rate
tempo_bpm = 100; % tempo
sec_per_quarter = 60/tempo_bpm;
tones_hz = [392 440 587.33];
snr_db = 40;

% welch (2 hz bins, 50% overlap, 100 blocks)
df_target_hz = 2;
N0_block = round(fs_hz/df_target_hz);
N_fft = 2^nextpow2(N0_block);
noverlap_welch = floor(N0_block/2);
K_blocks = 100;
hop_welch = N0_block - noverlap_welch;
L_needed = (K_blocks-1)*hop_welch + N0_block;

% integer number of quarter notes long enough
samp_per_quarter = round(sec_per_quarter*fs_hz);
M_notes = ceil(L_needed/samp_per_quarter);
L_total = M_notes*samp_per_quarter;

% synthesize
t = (0:L_total-1)/fs_hz;
sig_clean = zeros(1,L_total);

for q = 1:M_notes
    n0 = (q-1)*samp_per_quarter + 1;
    n1 = q*samp_per_quarter;
    seg_t = t(n0:n1);
    omit_idx = randi(3);
    keep_idx = setdiff(1:3, omit_idx);
    xq = cos(2*pi*tones_hz(keep_idx(1))*seg_t) + cos(2*pi*tones_hz(keep_idx(2))*seg_t);
    sig_clean(n0:n1) = xq;
end

% add white gaussian noise at target snr
p_sig = mean(sig_clean.^2);
p_noise = p_sig/10^(snr_db/10);
sig_noisy = sig_clean + sqrt(p_noise)*randn(size(sig_clean));


% periodogram

win_welch = hamming(N0_block);
[Pxx_welch, F_welch] = pwelch(sig_noisy, win_welch, noverlap_welch, N_fft, fs_hz);
pwr_db = 10*log10(Pxx_welch + eps);

% full 0..fs/2
figure('Color',[1 1 1]); hold on; grid on; box on;
plot(F_welch, pwr_db, 'LineWidth', 1.25);
xlim([0 fs_hz/2]);
xlabel('f (Hz)');
ylabel('Power (dB)');
title('Problem 4: Periodogram via pwelch (Hamming, ~2 Hz bins)');

% zoom 300-600 hz
figure('Color',[1 1 1]); hold on; grid on; box on;
plot(F_welch, pwr_db, 'LineWidth', 1.25);
xlim([300 600]);
xlabel('f (Hz)');
ylabel('Power (dB)');
title('Problem 4: Periodogram (Zoom 300-600 Hz)');


% spectrogram

% window/nfft same as welch; hop = eighth note
M_offset = round((sec_per_quarter/2)*fs_hz);
noverlap_spec = N0_block - M_offset;
if noverlap_spec < 0
    noverlap_spec = 0;
end

[S_stft, F_spec, T_spec] = spectrogram(sig_noisy, win_welch, noverlap_spec, N_fft, fs_hz);
S_db = 10*log10(abs(S_stft).^2 + eps);

% full band
figure('Color',[1 1 1]);
imagesc(T_spec, F_spec, S_db); axis xy; colormap jet; colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Problem 4: Spectrogram (window = N0, nfft = periodogram, offset = eighth note)');
ylim([0 fs_hz/2]);


%% Problem 5

% params
N = 200; % length
n = 0:N-1; % index grid

% (a) windows
w = ones(1,N); % rectangular on 0..N-1
w0 = zeros(1,N); % subsampling mask: 1 at n = 4m+3
w0(4:4:end) = 1; % (MATLAB index shift)
% relation (dtft): w0 = w .* p where p selects every 4th with offset 3

% (b) N-point dfts and stem plots
W_dft = fft(w,N);
W0_dft = fft(w0,N);
k = 0:N-1; % bin indices

figure('Color',[1 1 1]); hold on; grid on; box on;
stem(k, abs(W_dft), 'filled'); % |W[k]|
stem(k, abs(W0_dft), 'filled'); % |W0[k]|
xlabel('k'); ylabel('|W(k)|, |W_0(k)|'); 
title('Problem 5(b): N-point DFT magnitudes of w[n] and w_0[n]');
legend('|W(k)|','|W_0(k)|','Location','best');

% (c) high-resolution (N0-point) spectra, dc normalized
N0 = 2^nextpow2(16*N); % smallest power of 2 >= 16N
W_N0 = fft(w,N0); % high-res dft of w
W0_N0 = fft(w0,N0); % high-res dft of w0
W_N0 = W_N0 / W_N0(1); % normalize dc to 1
W0_N0 = W0_N0 / W0_N0(1); % normalize dc to 1

W_N0_s = fftshift(W_N0);
W0_N0_s = fftshift(W0_N0);
omega = (-N0/2:N0/2-1) * (2*pi/N0);

figure('Color',[1 1 1]); hold on; grid on; box on;
plot(omega, abs(W_N0_s), 'LineWidth', 1.2); % |W(e^{j\omega})|
plot(omega, abs(W0_N0_s), 'LineWidth', 1.2); % |W0(e^{j\omega})|
xlim([-pi pi]);
xlabel('\omega (rad/sample)');
ylabel('magnitude (dc-normalized)');
title('Problem 5(c): |W(e^{j\omega})| and |W_0(e^{j\omega})| (N_0-point)');
legend('|W|','|W_0|','Location','best');

% (d) cosine, window, zero-pad to N0, compare 0-pi
x = cos(0.4*n);
x_w = x .* w; % rectangular
x_w0 = x .* w0; % subsampled mask

X_w = fft([x_w zeros(1,N0-N)] , N0);
X_w0 = fft([x_w0 zeros(1,N0-N)], N0);
k_pos = 0:N0/2; % nonnegative bins
omega_pos = k_pos * (2*pi/N0);

figure('Color',[1 1 1]); hold on; grid on; box on;
plot(omega_pos, abs(X_w(k_pos+1)), 'LineWidth', 1.2); % |X_w|
plot(omega_pos, abs(X_w0(k_pos+1)), 'LineWidth', 1.2); % |X_{w0}|
xlim([0 pi]);
xlabel('\omega (rad/sample)');
ylabel('magnitude');
title('Problem 5(d): Cosine Spectrum with w[n] vs w_0[n] (N_0-point)');
legend('|X\_w|','|X\_{w0}|','Location','best');

% (e)
fprintf('\nProblem 5(e):\n');
fprintf('w0[n] keeps every fourth sample with an offset of three, which is periodic subsampling \n');
fprintf('this introduces spectral imaging: the magnitude spectrum shows extra image copies of the rectangular-window spectrum shifted by multiples of pi/2 \n');
fprintf('because w0[n] has only N/4 nonzero samples, its effective aperture is four times shorter, so the mainlobe is wider and the spectral resolution is lower than with w[n] \n');
fprintf('compared with w[n], spectra using w0[n] show additional images near +- pi/two and near pi, and the tone peaks are less sharp; w0[n] therefore causes imaging distortion and poorer resolution \n\n');