% Taha Akhlaq - DSP Problem Set 2
clear; close all; clc;

%% Problem 1:

% parameters
Fs = 40e6; % sampling rate
Rp = 1.5; % passband ripple
Rs = 30; % stopband attenuation
fp = [10e6 11e6]; % passband edges
fsb = [9e6  12e6]; % stopband edges
Fmax = 20e6; % plot limit
MAG_YLIM = [-80 5]; % magnitude axis
N = 4096; % frequency grid size
edgeHz = [9e6 10e6 11e6 12e6]; % edge freqs

% normalized (digital) and rad/s (analog) edges
Wp_d = fp/(Fs/2);
Ws_d = fsb/(Fs/2);
Wp_a = 2*pi*fp;
Ws_a = 2*pi*fsb;

% frequency grids for responses
f_lin = linspace(0, Fmax, N);
w_analog = 2*pi*f_lin;
w_dig = linspace(0, pi, N);
f_from_w = (w_dig/(2*pi))*Fs;

% 8 filters (analog/digital Butterworth, Cheb I, Cheb II, Elliptic)
F = {
  struct('ord',@() buttord(Wp_a,Ws_a,Rp,Rs,'s'),'des',@(n,W) butter(n,W,'bandpass','s'),'isA',true,'name',"Analog Butterworth")
  struct('ord',@() cheb1ord(Wp_a,Ws_a,Rp,Rs,'s'),'des',@(n,W) cheby1(n,Rp,W,'bandpass','s'),'isA',true,'name',"Analog Chebyshev I")
  struct('ord',@() cheb2ord(Wp_a,Ws_a,Rp,Rs,'s'),'des',@(n,W) cheby2(n,Rs,W,'bandpass','s'),'isA',true,'name',"Analog Chebyshev II")
  struct('ord',@() ellipord(Wp_a,Ws_a,Rp,Rs,'s'),'des',@(n,W) ellip(n,Rp,Rs,W,'bandpass','s'),'isA',true,'name',"Analog Elliptic")
  struct('ord',@() buttord(Wp_d,Ws_d,Rp,Rs),'des',@(n,W) butter(n,W,'bandpass'),'isA',false,'name',"Digital Butterworth")
  struct('ord',@() cheb1ord(Wp_d,Ws_d,Rp,Rs),'des',@(n,W) cheby1(n,Rp,W,'bandpass'),'isA',false,'name',"Digital Chebyshev I")
  struct('ord',@() cheb2ord(Wp_d,Ws_d,Rp,Rs),'des',@(n,W) cheby2(n,Rs,W,'bandpass'),'isA',false,'name',"Digital Chebyshev II")
  struct('ord',@() ellipord(Wp_d,Ws_d,Rp,Rs),'des',@(n,W) ellip(n,Rp,Rs,W,'bandpass'),'isA',false,'name',"Digital Elliptic")
};

% build responses first to set a common phase scale
for k = 1:numel(F)
    [n,Wn] = F{k}.ord();
    [b,a]  = F{k}.des(n,Wn);

    if F{k}.isA
        H = freqs(b,a,w_analog);   f = f_lin;
    else
        H = freqz(b,a,w_dig);      f = f_from_w;
    end

    D(k).name = F{k}.name;
    D(k).isA = F{k}.isA;
    D(k).b = b;
    D(k).a = a;
    D(k).f = f;
    D(k).mag = 20*log10(abs(H));      
    D(k).ph = unwrap(angle(H))*180/pi;  
    D(k).ordDen = numel(a)-1; % part (a) denominator degree
    D(k).z = roots(b);   
    D(k).p = roots(a);
end

phMin = min(cellfun(@min,{D.ph}));
phMax = max(cellfun(@max,{D.ph}));
pad = 0.05*(phMax - phMin + eps);
PHASE_YLIM = [phMin - pad, phMax + pad];

fprintf('\nProblem 1 (a–d) \n');
fprintf('Passband 10–11 MHz, Stopbands <9 & >12 MHz, Rp=%.2f dB, Rs=%.2f dB, Fs=%.1f MHz\n\n', Rp, Rs, Fs/1e6);

for k = 1:numel(D)
    % (c) magnitude, unwrapped phase, 0–20 MHz, uniform scales
      figure('Name',sprintf('%s — Response',D(k).name),'Color','w'); 
      tiledlayout(2,1);

      nexttile;
      plot(D(k).f*1e-6, D(k).mag,'LineWidth',1.2); grid on; xlim([0 20]); ylim(MAG_YLIM);
      xlabel('Frequency (MHz)'); ylabel('|H| (dB)');
      title(sprintf('%s — |H| (dB)  [Den. Order = %d]', D(k).name, D(k).ordDen));

      nexttile;
      plot(D(k).f*1e-6, D(k).ph,'LineWidth',1.2); grid on; xlim([0 20]); ylim(PHASE_YLIM);
      xlabel('Frequency (MHz)'); ylabel('Phase (degree)');
      title(sprintf('%s — Unwrapped Phase', D(k).name));

      % (b) pole–zero plots
      if D(k).isA
        figure('Name',sprintf('%s — s-plane PZ',D(k).name),'Color','w');
        hold on; grid on; axis equal;
        plot(real(D(k).z)/(2*pi*1e6), imag(D(k).z)/(2*pi*1e6), 'o','LineWidth',1.2);
        plot(real(D(k).p)/(2*pi*1e6), imag(D(k).p)/(2*pi*1e6), 'x','LineWidth',1.2);
        xlabel('Re\\{s\\}/(2\\pi·MHz)'); ylabel('Im\\{s\\}/(2\\pi·MHz)');
        title(sprintf('%s — Zeros (o) & Poles (x)', D(k).name));
        legend({'zeros','poles'},'Location','best');
      else
        figure('Name',sprintf('%s — z-plane PZ',D(k).name),'Color','w');
        zplane(D(k).z, D(k).p);
        title(sprintf('%s — z-plane Zeros & Poles', D(k).name));
      end

      % (d) Edge checks at 9, 10, 11, 12 MHz and passband peak in 10–11 MHz
      EdB = edge_dB(D(k).b, D(k).a, D(k).isA, edgeHz, Fs);
      passIdx = (D(k).f >= fp(1)) & (D(k).f <= fp(2));
      peakPB = max(D(k).mag(passIdx));

      v9  = verdict_stop(EdB(1), Rs);
      v12 = verdict_stop(EdB(4), Rs);
      v10 = verdict_pass(EdB(2), Rp);
      v11 = verdict_pass(EdB(3), Rp);
      vpk = verdict_peak(peakPB, Rp);

      fprintf('%-22s | DenOrder=%2d | Peak(10–11)= %+7.3f dB  -> %s\n', D(k).name, D(k).ordDen, peakPB, vpk);
      fprintf('  @ 9 MHz: %+7.3f dB -> %s   @ 10 MHz: %+7.3f dB -> %s\n', EdB(1), v9,  EdB(2), v10);
      fprintf('  @ 11 MHz: %+7.3f dB -> %s  @ 12 MHz: %+7.3f dB -> %s\n\n', EdB(3), v11, EdB(4), v12);
end

% helpers

% magnitude (dB) at edge frequencies
function dB = edge_dB(b,a,isAnalog,edgeHz,Fs)
if isAnalog
    H = freqs(b,a,2*pi*edgeHz);
else
    H = freqz(b,a,2*pi*(edgeHz/Fs));
end
dB = 20*log10(abs(H)).';
end

% stopband meets if amplitude <= -Rs
function s = verdict_stop(val,Rs)
epsdB = 1e-6;
if val <= -Rs - epsdB
    s = "overdesign";
elseif val <= -Rs + epsdB
    s = "meets";
else
    s = "violates";
end
end

% passband edge meets if -Rp <= value <= 0 dB
function s = verdict_pass(val,Rp)
epsdB = 1e-3;
if val >= -Rp - epsdB && val <= 0 + epsdB
    s = "meets";
elseif val < -Rp - epsdB
    s = "under";
else
    s = "over";
end
end

% passband peak target -Rp dB <= peak <= 0 dB
function s = verdict_peak(pk,Rp)
epsdB = 1e-3;
if pk >= -Rp - epsdB && pk <= 0 + epsdB
    s = "≈0 dB (ok)";
elseif pk < -Rp - epsdB
    s = "< -Rp (low)";
else
    s = "> 0 dB";
end
end

%{
- Peak in 10-11 MHz ~ 0 dB for all filters -> consistent with unity-gain
- Passband edges at 10 and 11 MHz are within +/-1.5 dB for all
- Stopband edges (9 and 12 MHz, target <= -30 dB):
  - Meets both edges exactly: Digital Chebyshev II (order 6)
  - Meets at 9, slight overdesign at 12: Digital Butterworth (order 8)
  - Overdesign at one or both edges: Analog Butterworth, Analog Chebyshev I,
    Analog Chebyshev II, Analog Elliptic, Digital Chebyshev I, Digital Elliptic
- extra attenuation is expected from integer order choices and
  the sharper transitions of Chebyshev/Elliptic families. Butterworth needs a
  higher order for the same specs
%}


%% Problem 2

% linear tolerances (passband, stopband)
delta_p = (10^(Rp/20)-1)/(10^(Rp/20)+1);
delta_s = 10^(-Rs/20);

% reuse grids
wz = w_dig;
fz = f_from_w;

% kaiser-window fir
[filtOrderK, WnK, betaK] = kaiserord([fsb(1) fp(1) fp(2) fsb(2)], [0 1 0], [delta_s delta_p delta_s], Fs);

% odd length for symmetric bandpass
if rem(filtOrderK,2) ~=0
    filtOrderK = filtOrderK + 1;
end

bK = fir1(filtOrderK, [fp(1) fp(2)]/(Fs/2), 'bandpass', kaiser(filtOrderK+1, betaK));

% coefficient stem
figure('Name','Kaiser FIR: coefficients','Color','w');
stem(0:numel(bK)-1, bK, 'filled'); grid on;
xlabel('n'); ylabel('h[n]');
title(sprintf('Kaiser FIR: length = %d', numel(bK)));

% z-plane
figure('Name','Kaiser FIR: z-plane','Color','w');
zplane(bK, 1);
title('Kaiser FIR: z-plane zeros and poles');

HK = freqz(bK, 1, wz);
magK = 20*log10(abs(HK));
phK = unwrap(angle(HK)) * 180/pi;

% equiripple (parks–mcclellan)
[Npm, foPM, aoPM, wPM] = firpmord([fsb(1) fp(1) fp(2) fsb(2)], [0 1 0], [delta_s delta_p delta_s], Fs);

% odd length
if rem(Npm, 2) ~= 0
    Npm = Npm + 1;
end

bPM = firpm(Npm, foPM, aoPM, wPM);  

% coefficient stem
figure('Name','Equiripple FIR: coefficients','Color','w');
stem(0:numel(bPM)-1, bPM, 'filled'); grid on;
xlabel('n'); ylabel('h[n]');
title(sprintf('Equiripple FIR: length = %d', numel(bPM)));

% z-plane
figure('Name','Equiripple FIR: z-plane','Color','w');
zplane(bPM, 1);
title('Equiripple FIR: z-plane zeros and poles');

HPM = freqz(bPM, 1, wz);
magPM = 20*log10(abs(HPM));
phPM = unwrap(angle(HPM)) * 180/pi;

% plots
phMin2 = min([min(phK) min(phPM)]);
phMax2 = max([max(phK) max(phPM)]);
pad2 = 0.05*(phMax2 - phMin2 + eps);
PHASE2_YLIM = [phMin2 - pad2, phMax2 + pad2];

FIRs = {
  struct('name',"Kaiser FIR", 'mag',magK, 'ph',phK)
  struct('name',"Equiripple FIR", 'mag',magPM, 'ph',phPM)
};

for i = 1:numel(FIRs)
    figure('Name',sprintf('%s — Response',FIRs{i}.name),'Color','w');
    tiledlayout(2,1);

    nexttile;
    plot(fz*1e-6, FIRs{i}.mag, 'LineWidth', 1.2);
    grid on; xlim([0 20]); ylim(MAG_YLIM);
    xlabel('Frequency (MHz)'); ylabel('|H| (dB)');
    title(sprintf('%s — |H| (dB)', FIRs{i}.name));

    nexttile;
    plot(fz*1e-6, FIRs{i}.ph, 'LineWidth', 1.2);
    grid on; xlim([0 20]); ylim(PHASE2_YLIM);
    xlabel('Frequency (MHz)'); ylabel('Phase (degree)');
    title(sprintf('%s — Unwrapped Phase', FIRs{i}.name));
end

% 2(b) equiripple weights vs tolerances
W_expected = [1/delta_s, 1/delta_p, 1/delta_s];
wRatios = wPM / wPM(2);
fprintf('\nProblem 2(b)\n');
fprintf('delta_p = %.6g, delta_s = %.6g\n', delta_p, delta_s);
fprintf('expected weight ratios = [%.3f 1.000 %.3f]\n', W_expected(1)/W_expected(2), W_expected(3)/W_expected(2));
fprintf('firpmord weight ratios = [%.3f 1.000 %.3f]\n\n', wRatios(1), wRatios(3));

% 2(c) fir passband can not peak at 0 dB
RK  = check_fir_specs(fz,  magK,  fp, fsb, Rp, Rs);
RPM = check_fir_specs(fz,  magPM, fp, fsb, Rp, Rs);

fprintf('Problem 2(c)\n');
% Kaiser
fprintf('Kaiser:\n');
fprintf('  passband max/min  = %+6.3f / %+6.3f dB\n', RK.pbMax, RK.pbMin);
fprintf('  passband spread   = %6.3f dB  (target 1.5)  -> %s\n', RK.pbDiff, tern(RK.pb_ok,'meets','fails'));
fprintf('  stop peaks  <=9 MHz:  %+6.3f dB   >=12 MHz:  %+6.3f dB   (limit -30)\n', RK.s1Peak, RK.s2Peak);
fprintf('  stop verdicts      <=9: %s   >=12: %s\n\n', RK.s1_v, RK.s2_v);

% Equiripple
fprintf('Equiripple:\n');
fprintf('  passband max/min  = %+6.3f / %+6.3f dB\n', RPM.pbMax, RPM.pbMin);
fprintf('  passband spread   = %6.3f dB  (target 1.5)  -> %s\n', RPM.pbDiff, tern(RPM.pb_ok,'meets','fails'));
fprintf('  stop peaks  <=9 MHz:  %+6.3f dB   >=12 MHz:  %+6.3f dB   (limit -30)\n', RPM.s1Peak, RPM.s2Peak);
fprintf('  stop verdicts      <=9: %s   >=12: %s\n\n', RPM.s1_v, RPM.s2_v);

% helpers
function R = check_fir_specs(fz, magdB, fp, fsb, Rp, Rs)
    inPB = (fz >= fp(1)) & (fz <= fp(2));
    inS1 = (fz <= fsb(1));
    inS2 = (fz >= fsb(2));

    R.pbMax  = max(magdB(inPB));
    R.pbMin  = min(magdB(inPB));
    R.pbDiff = R.pbMax - R.pbMin;

    R.s1Peak = max(magdB(inS1));
    R.s2Peak = max(magdB(inS2));

    R.pb_ok  = (R.pbDiff <= Rp + 1e-3);
    R.s1_v   = stop_verdict(R.s1Peak, Rs);
    R.s2_v   = stop_verdict(R.s2Peak, Rs);
end

function s = stop_verdict(val, Rs)
    epsdB = 1e-6;
    if val <= -Rs - epsdB
        s = "overdesign";
    elseif val <= -Rs + epsdB
        s = "meets";
    else
        s = "fails";
    end
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end

%{
Kaiser
- passband max/min = +0.000 dB / -5.327 dB
- passband spread  = 5.327 dB vs 1.5 dB => FAILS
- stop peaks: <=9 MHz = -33.603 dB, >=12 MHz = -35.234 dB vs -30 dB => OVERDESIGN
- fix: increase length and/or widen transitions

Equiripple
- passband max/min = +0.912 dB / -1.008 dB
- passband spread  = 1.920 dB vs 1.5 dB => FAILS
- stop peaks: <=9 MHz = -27.765 dB, >=12 MHz = -27.794 dB vs -30 dB => FAILS
- fix: raise order and increase stopband weights
%}