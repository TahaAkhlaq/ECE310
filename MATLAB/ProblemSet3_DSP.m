% Taha Akhlaq - DSP Problem Set 3
clear; close all; clc;

%% Problem 4:

% parameters
Rp = 1.5; % passband ripple
Rs = 30; % stopband attenuation
Wp = [0.3 0.6]; % passband edges
Ord = 8; % bandpass order

% 4a
Nbp = Ord/2;
[z,p,k] = ellip(Nbp, Rp, Rs, Wp, 'bandpass');

if length(p) ~= Ord
    error('ellip returned order %d, expected %d', length(p), Ord);
end

[b,a] = zp2tf(z,p,k);

nfft = 4096;
[H,w] = freqz(b,a,nfft);

figure;
plot(w/pi, 20*log10(abs(H)),'LineWidth',1.1);
grid on;
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
title('4a: 8th-order Elliptic BPF');

% 4b
[sos_up, G_up] = zp2sos(z,p,k,'up','inf');
[sos_dn, G_dn] = zp2sos(z,p,k,'down','inf');

L = size(sos_up,1);

[b_up,a_up] = sos2tf(sos_up,G_up);
[b_dn,a_dn] = sos2tf(sos_dn,G_dn);

fprintf('4b: SOS vs original H(z)\n');
fprintf('  up:   max|b-b| = %.5e, max|a-a| = %.5e\n', max(abs(b - b_up)), max(abs(a - a_up)));
fprintf('  down: max|b-b| = %.5e, max|a-a| = %.5e\n', max(abs(b - b_dn)), max(abs(a - a_dn)));

% 4c
pz_up = struct('poles',cell(L,1),'zeros',cell(L,1));
pz_dn = struct('poles',cell(L,1),'zeros',cell(L,1));

for i = 1:L
    % up
    b_i = sos_up(i,1:3);
    a_i = sos_up(i,4:6);
    z_i = roots(b_i);
    p_i = roots(a_i);
    pz_up(i).poles = sort(abs(p_i)).';
    pz_up(i).zeros = sort(abs(z_i)).';

    % down
    b_i = sos_dn(i,1:3);
    a_i = sos_dn(i,4:6);
    z_i = roots(b_i);
    p_i = roots(a_i);
    pz_dn(i).poles = sort(abs(p_i)).';
    pz_dn(i).zeros = sort(abs(z_i)).';
end

fprintf('\n4c: pole magnitudes (up):\n');
for i = 1:L
    fprintf('  %d: [%g, %g]\n', i, pz_up(i).poles);
end

fprintf('\n4c: pole magnitudes (down):\n');
for i = 1:L
    fprintf('  %d: [%g, %g]\n', i, pz_dn(i).poles);
end

% 4d
w = linspace(0,pi,nfft);
z_w = exp(1j*w);

[F_up,  ~, ~] = compute_Fk_sos(sos_up, G_up, z_w);
[F_dn,  ~, ~] = compute_Fk_sos(sos_dn, G_dn, z_w);

maxF_up = max(abs(F_up),[],2).';
maxF_dn = max(abs(F_dn),[],2).';

fprintf('\n4d: max |F_k| (up):\n');
disp(maxF_up);
fprintf('4d: max |F_k| (down):\n');
disp(maxF_dn);

figure;
hold on;
for k = 1:L
    plot(w/pi, 20*log10(abs(F_up(k,:))),'LineWidth',1.0);
end
hold off;
grid on;
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
title('4d: Cumulative F_k(z) for up ordering');
legend(arrayfun(@(k) sprintf('F_%d(z)',k), 1:L, 'UniformOutput',false), ...
       'Location','Best');

figure;
hold on;
for k = 1:L
    plot(w/pi, 20*log10(abs(F_dn(k,:))),'LineWidth',1.0);
end
hold off;
grid on;
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
title('4d: Cumulative F_k(z) for down ordering');
legend(arrayfun(@(k) sprintf('F_%d(z)',k), 1:L, 'UniformOutput',false), ...
       'Location','Best');

% 4e
tol = 1e-6;
fprintf('4e: up vs down SOS (reverse-order)\n');

for i = 1:L
    j  = L + 1 - i;

    Au = sos_up(i,4:6);
    Ad = sos_dn(j,4:6);
    Bu = sos_up(i,1:3);
    Bd = sos_dn(j,1:3);

    den_err = max(abs(Au - Ad));

    c = (Bd * Bu.') / (Bu * Bu.');
    num_err = max(abs(Bd - c*Bu));

    fprintf('  up %d vs down %d: ', i, j);
    fprintf('denominator_error = %.3e, ', den_err);
    fprintf('c = %.6g, ', c);
    fprintf('numerator_error = %.3e\n', num_err);
end

%% Problem 5

% original H(z)
bH = [0.1336 0.0563 0.0563 0.1336];
aH = [1 -1.5055 1.2630 -0.3778];

% allpass form
bA1 = [-0.4954 1];
aA1 = [1 -0.4954];
bA2 = [0.7626 -1.0101 1];
aA2 = [1 -1.0101 0.7626];

[bHa,aHa] = parallel_allpass(bA1,aA1,bA2,aA2);

% 5a
Nw = 4096;
w  = linspace(0,pi,Nw);

H  = freqz(bH, aH, w);
Ha = freqz(bHa,aHa,w);

mag_err = max(abs(abs(H) - abs(Ha)));

pH  = sort(roots(aH));   zH  = sort(roots(bH));
pHa = sort(roots(aHa));  zHa = sort(roots(bHa));

pole_err = max(abs(pH - pHa));
zero_err = max(abs(zH - zHa));

fprintf('\n5a:\n');
fprintf('  max |H|-|Ha|   = %.5e\n', mag_err);
fprintf('  max pole error = %.5e\n', pole_err);
fprintf('  max zero error = %.5e\n', zero_err);

% 5b
b = 4;

[qH,f_H]   = quantize_coef([bH aH],  b);
[bHq,aHq]  = split_coef(qH,length(bH));

[qA1,f_A1] = quantize_coef([bA1 aA1],b);
[bA1q,aA1q]= split_coef(qA1,length(bA1));

[qA2,f_A2] = quantize_coef([bA2 aA2],b);
[bA2q,aA2q]= split_coef(qA2,length(bA2));

[bHAQ,aHAQ] = parallel_allpass(bA1q,aA1q,bA2q,aA2q);

fprintf('\n5b:\n');
fprintf('  H(z) frac bits = %d\n', f_H);
fprintf('  allpass1 frac bits = %d\n', f_A1);
fprintf('  allpass2 frac bits = %d\n', f_A2);

% 5c
p      = roots(aH);
z0     = roots(bH);
p_qH   = roots(aHq);
z_qH   = roots(bHq);
p_qAQ  = roots(aHAQ);
z_qAQ  = roots(bHAQ);

pe_qH  = max(abs(sort(p)  - sort(p_qH)));
pe_qAQ = max(abs(sort(p)  - sort(p_qAQ)));
ze_qH  = max(abs(sort(z0) - sort(z_qH)));
ze_qAQ = max(abs(sort(z0) - sort(z_qAQ)));

figure; hold on; grid on; axis equal;
th = linspace(0,2*pi,400);
plot(cos(th), sin(th), 'k--');

plot(real(z0),    imag(z0),    'bo', 'MarkerSize',6);
plot(real(p),     imag(p),     'bx', 'MarkerSize',6);

plot(real(z_qH),  imag(z_qH),  'ro', 'MarkerSize',6);
plot(real(p_qH),  imag(p_qH),  'rx', 'MarkerSize',6);

plot(real(z_qAQ), imag(z_qAQ), 'go', 'MarkerSize',6);
plot(real(p_qAQ), imag(p_qAQ), 'gx', 'MarkerSize',6);

title('5c: Poles and zeros');
xlabel('Real Part');
ylabel('Imaginary Part');
legend('Unit circle', ...
       'H zeros','H poles', ...
       'H_Q zeros','H_Q poles', ...
       'H_{AQ} zeros','H_{AQ} poles', ...
       'Location','Best');
hold off;

fprintf('\n5c:\n');
fprintf('  max pole error (H quantized)       = %.3e\n', pe_qH);
fprintf('  max pole error (Ha quantized)      = %.3e\n', pe_qAQ);
fprintf('  max zero error (H quantized)       = %.3e\n', ze_qH);
fprintf('  max zero error (Ha quantized)      = %.3e\n', ze_qAQ);

% 5d
Href0  = eval_H(bH,  aH,  1);
Hrefpi = eval_H(bH,  aH, -1);

HQ0    = eval_H(bHq, aHq, 1);
HQpi   = eval_H(bHq, aHq,-1);
HAQ0   = eval_H(bHAQ,aHAQ,1);
HAQpi  = eval_H(bHAQ,aHAQ,-1);

eHQ0_dB   = gain_err_db(Href0,HQ0);
eHQpi_dB  = gain_err_db(Hrefpi,HQpi);
eHAQ0_dB  = gain_err_db(Href0,HAQ0);
eHAQpi_dB = gain_err_db(Hrefpi,HAQpi);

fprintf('\n5d\n');
fprintf('  H quantized  : w=0: %+6.3f dB, w=pi: %+6.3f dB\n',  eHQ0_dB,  eHQpi_dB);
fprintf('  Ha quantized : w=0: %+6.3f dB, w=pi: %+6.3f dB\n', eHAQ0_dB, eHAQpi_dB);

% 5e
Nw = 1e4;
w  = linspace(0,pi,Nw);

H    = freqz(bH,   aH,   w);
HQ   = freqz(bHq,  aHq,  w);
Ha1q = freqz(bA1q, aA1q, w);
Ha2q = freqz(bA2q, aA2q, w);
HAQ  = 0.5*(Ha1q + Ha2q);

e_HQ  = max(abs(H - HQ));
e_HAQ = max(abs(H - HAQ));

fprintf('\n5e:\n');
fprintf('  max |H - H quantized|  = %.3e\n',  e_HQ);
fprintf('  max |H - Ha quantized| = %.3e\n', e_HAQ);

% 5f
H_dB   = 20*log10(abs(H));
HQ_dB  = 20*log10(abs(HQ));
HAQ_dB = 20*log10(abs(HAQ));

H_ph   = unwrap(angle(H))*180/pi;
HQ_ph  = unwrap(angle(HQ))*180/pi;
HAQ_ph = unwrap(angle(HAQ))*180/pi;

figure;
plot(w/pi,H_dB,'LineWidth',1.1); hold on;
plot(w/pi,HQ_dB,'LineWidth',1.1);
plot(w/pi,HAQ_dB,'LineWidth',1.1);
hold off; grid on;
ylim([-40, max(H_dB)+2]);
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Magnitude (dB)');
title('5f: Magnitude responses');
legend('H','H_Q','H_AQ','Location','Best');

figure;
plot(w/pi,H_ph,'LineWidth',1.1); hold on;
plot(w/pi,HQ_ph,'LineWidth',1.1);
plot(w/pi,HAQ_ph,'LineWidth',1.1);
hold off; grid on;
xlabel('Normalized Frequency (x pi rad/sample)');
ylabel('Phase (degrees)');
title('5f: Phase responses');
legend('H','H_Q','H_AQ','Location','Best');

%% local functions

function [F_all,B_all,A_all] = compute_Fk_sos(sos,G,z_w)
L = size(sos,1);
nfft = length(z_w);
z_inv  = 1./z_w;
z_inv2 = z_inv.^2;
B_all = zeros(L,nfft);
A_all = zeros(L,nfft);
for i = 1:L
    b = sos(i,1:3);
    a = sos(i,4:6);
    B_all(i,:) = b(1) + b(2).*z_inv + b(3).*z_inv2;
    A_all(i,:) = a(1) + a(2).*z_inv + a(3).*z_inv2;
end
F_all = zeros(L,nfft);
F_all(1,:) = G ./ A_all(1,:);
for k = 2:L
    F_all(k,:) = F_all(k-1,:) .* B_all(k-1,:) ./ A_all(k,:);
end
end

function [bHa,aHa] = parallel_allpass(b1,a1,b2,a2)
bHa = 0.5*(conv(b1,a2) + conv(b2,a1));
aHa = conv(a1,a2);
end

function [q,f] = quantize_coef(coef,bits)
maxc = max(abs(coef));
for f = bits-1:-1:0
    max_val = (2^(bits-1)-1)/2^f;
    if maxc <= max_val
        break;
    end
end
q_int = round(coef*2^f);
q_int = max(min(q_int, 2^(bits-1)-1), -2^(bits-1));
q = q_int / 2^f;
end

function [bq,aq] = split_coef(q,n_b)
bq = q(1:n_b);
aq = q(n_b+1:end);
end

function Hval = eval_H(b,a,z)
k   = 0:length(b)-1;
num = b * (z.^(-k)).';
k   = 0:length(a)-1;
den = a * (z.^(-k)).';
Hval = num/den;
end

function err = gain_err_db(Href,Htest)
Href  = max(abs(Href),eps);
Htest = max(abs(Htest),eps);
err = 20*log10(Htest/Href);
end