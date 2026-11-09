% CW1_LPC_VOWELS  EEEM030 Assignment 1: linear predictive vowel synthesiser.
%
% This script does the whole coursework in Matlab: it reads the vowel
% recordings, estimates the LPC filters, finds the formants and F0,
% compares model orders and segment lengths, synthesises the vowels and
% writes the figures, tables and wav files.
%
% How to run: open this file in Matlab and press Run. It looks for the folder
% data/speech (next to this file, or one to three folders above it) and
% writes its results to a folder called output beside that data folder:
%   output/figures/   the two figure types used in the report
%   output/tables/    csv tables, results.json and summary.txt (console log)
%   output/audio/     the synthesised wav files
% To use another data folder, set the environment variable CW1_ROOT.
%
% What it does, in order:
%   1. Six self-tests of the functions at the end of this file.
%   2. For whod_m, whod_f (the two main vowels) and hod_f (an extra, no audio):
%      picks a 100 ms segment, estimates the order 26 LPC filter with the
%      autocorrelation method (own Levinson-Durbin recursion), reads off the
%      formants from the poles, estimates F0, checks against the covariance
%      method, shifted windows and a cepstral envelope, sweeps the model order
%      and the segment length, tests pre-emphasis and 10 kHz resampling,
%      synthesises one second of the vowel and analyses it again.
%   3. Saves everything the report uses.
%
% Only base Matlab is needed. If the Signal Processing Toolbox is installed,
% the script also compares its own coefficients with the built-in lpc() and
% arcov() and prints the differences. Local functions at the end of the file
% need Matlab R2016b or later.

clear; close all; clc;

rootDir = getenv('CW1_ROOT');
if isempty(rootDir)
    here = fileparts(mfilename('fullpath'));
    rootDir = '';
    cand = here;
    for up = 0:3                                   % look for data/speech here and in the folders above
        if exist(fullfile(cand, 'data', 'speech'), 'dir'), rootDir = cand; break; end
        cand = fileparts(cand);
    end
end
if isempty(rootDir) || ~exist(fullfile(rootDir, 'data', 'speech'), 'dir')
    rootDir = uigetdir(pwd, 'Select the folder that contains data/speech');
end
dataDir = fullfile(rootDir, 'data', 'speech');
outDir = fullfile(rootDir, 'output');
figDir = fullfile(outDir, 'figures');
tabDir = fullfile(outDir, 'tables');
audDir = fullfile(outDir, 'audio');
for d = {outDir, figDir, tabDir, audDir}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

% ---- Settings --------------------------------------------------------
vowels   = {'whod_m', 'whod_f', 'hod_f'};   % file names without .wav
speaker  = {'male', 'female', 'female'};
isMain   = [true true false];               % hod_f is an extra: no audio, no figures
segMs    = 100;                             % segment length asked for in the brief
finalP   = 26;                              % fs/1000 + 2 for fs = 24 kHz
orders   = [4 8 12 16 20 26 32 40];
lengths  = [20 30 50 100 150 200];
shiftsMs = -20:5:20;
preCoef  = 0.97;
synthSec = 1.0;
nfft     = 8192;
cepK     = 40;                              % cepstrum coefficients kept for the envelope

fid = fopen(fullfile(tabDir, 'summary.txt'), 'w');
say(fid, 'EEEM030 CW1 Matlab results (order %d, %d ms segments)\n', finalP, segMs);
if exist('OCTAVE_VERSION', 'builtin') == 0
    say(fid, 'Matlab version: %s\n', version);
else
    say(fid, 'Run in Octave %s (not Matlab)\n', version);
end
say(fid, '\n');

% ---- 1. Self-tests -----------------------------------------------------
nPass = run_self_tests(fid);
say(fid, '%d of 6 self-tests passed\n\n', nPass);

orderRows = {}; lengthRows = {}; variantRows = {}; shiftRows = {};
res = struct();

for v = 1:numel(vowels)
    name = vowels{v};
    [x, fs] = audioread(fullfile(dataDir, [name '.wav']));
    x = x(:).';                                    % row vector
    n0 = round(segMs / 1000 * fs);

    % ---- Segment: most periodic and steady window near the middle ------
    s0 = choose_segment(x, fs, n0);
    seg = x(s0 + 1 : s0 + n0);
    centre = s0 + floor(n0 / 2);
    say(fid, '== %s (%s): segment %.3f to %.3f s, fs %d Hz\n', name, speaker{v}, s0/fs, (s0+n0)/fs, fs);

    % ---- LPC estimation (autocorrelation method) ----------------------
    [a, G, pg] = lpc_autocorr(seg, finalP, true);
    say(fid, 'Order %d: prediction gain %.1f dB, G = %.5f, max pole radius %.4f\n', finalP, pg, G, max(abs(roots(a))));
    if exist('lpc', 'file') == 2 || exist('lpc', 'file') == 6
        try
            aRef = lpc((seg .* hamming_win(n0)).', finalP);
            say(fid, 'Built-in lpc(): max |a - a_lpc| = %.2e\n', max(abs(a - aRef(:).')));
        catch
            say(fid, 'lpc() could not be used, built-in check skipped\n');
        end
    end

    % ---- Covariance method (own least squares) -------------------------
    [aCov, gCov] = lpc_covariance(seg, finalP);
    fmCov = formants_from_lpc(aCov, fs);
    if exist('arcov', 'file') == 2
        try
            aArc = arcov(seg.', finalP);
            say(fid, 'Built-in arcov(): max |a_cov - a_arcov| = %.2e\n', max(abs(aCov - aArc(:).')));
        catch
            say(fid, 'arcov() could not be used, built-in check skipped\n');
        end
    end

    % ---- Formants ------------------------------------------------------
    [fm, bw] = formants_from_lpc(a, fs);
    [fGrid, envA] = lpc_response_db(a, G, fs, nfft);
    [~, envC] = lpc_response_db(aCov, gCov, fs, nfft);
    m5 = fGrid <= 5000;
    dEnv = envA(m5) - envC(m5);
    offsetDb = mean(dEnv);
    rmsAfter = sqrt(mean((dEnv - offsetDb).^2));
    pkLpc = envelope_peaks(fGrid, envA, 150, 5000);
    pkCep = cepstral_envelope_peaks(seg, fs, cepK, nfft);
    say(fid, 'Formants F1-F3 = %s Hz; bandwidths %s Hz\n', num2str(round(fm)), num2str(round(bw)));
    say(fid, 'Covariance formants %s Hz; envelope offset %.2f dB, RMS after offset %.2f dB\n', ...
        num2str(round(fmCov)), offsetDb, rmsAfter);
    say(fid, 'LPC curve peaks <5 kHz: %s Hz\n', num2str(round(pkLpc(1:min(5,end)))));
    say(fid, 'Cepstral envelope (%d coefficients) peaks <5 kHz: %s Hz\n', cepK, num2str(round(pkCep(1:min(5,end)))));

    % Shifted windows
    sh = nan(numel(shiftsMs), 3);
    for i = 1:numel(shiftsMs)
        sg = segment_at(x, centre + round(shiftsMs(i)/1000*fs), n0);
        sh(i, :) = formants_from_lpc(lpc_autocorr(sg, finalP, true), fs);
        shiftRows{end+1} = sprintf('%s,%d,%s,%s,%s', name, shiftsMs(i), numstr(sh(i,1)), numstr(sh(i,2)), numstr(sh(i,3))); %#ok<SAGROW>
    end
    say(fid, 'Window shifts: min %s, max %s Hz\n', num2str(round(min(sh, [], 1))), num2str(round(max(sh, [], 1))));

    % ---- F0 ------------------------------------------------------------
    [~, f0frames, pk] = f0_track(seg, fs);
    voiced = pk >= 0.5;
    f0 = mean(f0frames(voiced));
    f0cep = f0_cepstrum(seg, fs);
    say(fid, 'F0: autocorrelation %.1f Hz (%d frames, sd %.1f), cepstrum %.1f Hz\n', f0, sum(voiced), std(f0frames(voiced)), f0cep);

    % ---- Order sweep ---------------------------------------------------
    say(fid, 'Order sweep: order, gain dB, F1, F2, F3\n');
    pgs = zeros(size(orders));
    envOrders = cell(size(orders));
    for i = 1:numel(orders)
        [ap, gp, pgs(i)] = lpc_autocorr(seg, orders(i), true);
        fp = formants_from_lpc(ap, fs);
        [~, envOrders{i}] = lpc_response_db(ap, gp, fs, nfft);
        orderRows{end+1} = sprintf('%s,%d,%.1f,%s,%s,%s,%.4f', name, orders(i), pgs(i), numstr(fp(1)), numstr(fp(2)), numstr(fp(3)), max(abs(roots(ap)))); %#ok<SAGROW>
        say(fid, '  %2d %5.1f  %s\n', orders(i), pgs(i), num2str(round(fp)));
    end

    % ---- Length sweep --------------------------------------------------
    say(fid, 'Length sweep (order %d): ms, gain dB, F1, F2, F3\n', finalP);
    envLengths = cell(size(lengths));
    for i = 1:numel(lengths)
        nL = round(lengths(i) / 1000 * fs);
        [sg, sgStart] = segment_at(x, centre, nL);
        [al, gl, pgl] = lpc_autocorr(sg, finalP, true);
        fl = formants_from_lpc(al, fs);
        [~, envLengths{i}] = lpc_response_db(al, gl, fs, nfft);
        lengthRows{end+1} = sprintf('%s,%d,%d,%.3f,%.1f,%s,%s,%s', name, lengths(i), nL, sgStart/fs, pgl, numstr(fl(1)), numstr(fl(2)), numstr(fl(3))); %#ok<SAGROW>
        say(fid, '  %3d %5.1f  %s\n', lengths(i), pgl, num2str(round(fl)));
    end

    % ---- Variants: pre-emphasis and 10 kHz -----------------------------
    segPre = pre_emphasis(seg, preCoef);
    seg10 = resample_24k_to_10k(seg);
    fs10 = fs * 5 / 12;
    variants = { ...
        '24 kHz, order 26 (final)',      formants_from_lpc(a, fs); ...
        '24 kHz, order 26, pre-emph.',   formants_from_lpc(lpc_autocorr(segPre, 26, true), fs); ...
        '24 kHz, order 12',              formants_from_lpc(lpc_autocorr(seg, 12, true), fs); ...
        '24 kHz, order 12, pre-emph.',   formants_from_lpc(lpc_autocorr(segPre, 12, true), fs); ...
        '10 kHz, order 12',              formants_from_lpc(lpc_autocorr(seg10, 12, true), fs10); ...
        '10 kHz, order 26',              formants_from_lpc(lpc_autocorr(seg10, 26, true), fs10)};
    for i = 1:size(variants, 1)
        fv = variants{i, 2};
        variantRows{end+1} = sprintf('%s,"%s",%s,%s,%s', name, variants{i,1}, numstr(fv(1)), numstr(fv(2)), numstr(fv(3))); %#ok<SAGROW>
        say(fid, '%s: %s\n', variants{i,1}, num2str(round(fv)));
    end

    % ---- Synthesis and whole-chain check (main vowels) -------------------
    exc = impulse_train(f0, fs, synthSec);
    y = normalise_peak(filter(G, a, exc), 0.9);
    mid = floor(numel(y) / 2);
    ySeg = y(mid - floor(n0/2) + 1 : mid - floor(n0/2) + n0);
    aChk = lpc_autocorr(ySeg, finalP, true);
    fmChk = formants_from_lpc(aChk, fs);
    [~, f0y, pky] = f0_track(ySeg, fs);
    f0Chk = mean(f0y(pky >= 0.5));
    say(fid, 'Synthetic vowel re-analysed: F1-F3 = %s Hz, F0 = %.1f Hz, %d impulses\n\n', num2str(round(fmChk)), f0Chk, sum(exc));

    if isMain(v)
        write_wav(fullfile(audDir, [name '_synth_final.wav']), y, fs);
        write_wav(fullfile(audDir, [name '_original_segment.wav']), normalise_peak(seg, 0.9), fs);
        write_wav(fullfile(audDir, [name '_impulse_train.wav']), 0.5 * exc, fs);
        for p = [4 8 12 26 40]
            [ap, gp] = lpc_autocorr(seg, p, true);
            write_wav(fullfile(audDir, sprintf('%s_synth_order%02d.wav', name, p)), normalise_peak(filter(gp, ap, exc), 0.9), fs);
        end
        for L = [20 50 100 200]
            sg = segment_at(x, centre, round(L / 1000 * fs));
            [al, gl] = lpc_autocorr(sg, finalP, true);
            write_wav(fullfile(audDir, sprintf('%s_synth_len%03dms.wav', name, L)), normalise_peak(filter(gl, al, exc), 0.9), fs);
        end
        [aPre, gPre] = lpc_autocorr(segPre, finalP, true);
        yPre = filter(1, [1, -preCoef], filter(gPre, aPre, exc));   % de-emphasis
        write_wav(fullfile(audDir, [name '_synth_preemph.wav']), normalise_peak(yPre, 0.9), fs);
        [a10, g10] = lpc_autocorr(seg10, 12, true);
        write_wav(fullfile(audDir, [name '_synth_10kHz_order12.wav']), normalise_peak(filter(g10, a10, impulse_train(f0, fs10, synthSec)), 0.9), fs10);

        make_figures(figDir, name, seg, fs, s0, segMs, finalP, fm, fSeg_spec(seg, fs, nfft), fGrid, envA, ...
            orders, pgs, envOrders);
        make_synthesis_figure(figDir, name, seg, ySeg, exc(mid - floor(n0/2) + 1 : mid - floor(n0/2) + n0), ...
            fs, f0, fSeg_spec(seg, fs, nfft), fGrid, envA, nfft);
    end

    % ---- Collect results -------------------------------------------------
    r = struct();
    r.speaker = speaker{v};
    r.fs = fs;
    r.segment_start_s = s0 / fs;
    r.segment_end_s = (s0 + n0) / fs;
    r.order = finalP;
    r.gain = G;
    r.prediction_gain_db = pg;
    r.max_pole_radius = max(abs(roots(a)));
    r.a_coefficients = a;
    r.formants_hz = fm;
    r.formant_bandwidths_hz = bw;
    r.formants_covariance_hz = fmCov;
    r.covariance_envelope_offset_db_0_5k = offsetDb;
    r.covariance_envelope_rms_after_offset_db_0_5k = rmsAfter;
    r.envelope_peaks_hz_0_5k = pkLpc(1:min(5, numel(pkLpc)));
    r.cepstral_envelope_k40_peaks_hz = pkCep(1:min(5, numel(pkCep)));
    r.formants_shift_min_hz = min(sh, [], 1);
    r.formants_shift_max_hz = max(sh, [], 1);
    r.f0_autocorr_mean_hz = f0;
    r.f0_autocorr_sd_hz = std(f0frames(voiced));
    r.f0_frames = sum(voiced);
    r.f0_cepstrum_hz = f0cep;
    r.synthetic_formants_reestimated_hz = fmChk;
    r.synthetic_f0_hz = f0Chk;
    r.impulse_count = sum(exc);
    res.(name) = r;
end

% ---- Save tables and results -------------------------------------------
write_lines(fullfile(tabDir, 'order_sweep.csv'), 'vowel,order,prediction_gain_db,F1,F2,F3,max_pole_radius', orderRows);
write_lines(fullfile(tabDir, 'length_sweep.csv'), 'vowel,length_ms,samples,start_s,prediction_gain_db,F1,F2,F3', lengthRows);
write_lines(fullfile(tabDir, 'variants.csv'), 'vowel,variant,F1,F2,F3', variantRows);
write_lines(fullfile(tabDir, 'window_shift_formants.csv'), 'vowel,shift_ms,F1,F2,F3', shiftRows);
try
    txt = jsonencode(res, 'PrettyPrint', true);
catch
    txt = jsonencode(res);
end
fj = fopen(fullfile(tabDir, 'results.json'), 'w'); fwrite(fj, txt); fclose(fj);
fclose(fid);
fprintf('Done. Outputs are in %s\n', outDir);


% =======================================================================
% Local functions
% =======================================================================
function say(fid, fmt, varargin)
    % Print to the command window and to the summary file.
    fprintf(fmt, varargin{:});
    fprintf(fid, fmt, varargin{:});
end

function s = numstr(v)
    % Formant value for the csv files: rounded to whole hertz, NaN if none.
    if isnan(v), s = 'NaN'; else, s = sprintf('%d', round(v)); end
end

function write_lines(path, header, rows)
    f = fopen(path, 'w');
    fprintf(f, '%s\n', header);
    for i = 1:numel(rows), fprintf(f, '%s\n', rows{i}); end
    fclose(f);
end

function write_wav(path, y, fs)
    % audiowrite wants one column per channel, so the row vector is transposed.
    audiowrite(path, y(:), fs, 'BitsPerSample', 16);
end

function y = normalise_peak(y, peak)
    y = y * (peak / max(abs(y)));
end

function w = hamming_win(n)
    % Symmetric Hamming window (the same as hamming(n) in the toolbox).
    w = 0.54 - 0.46 * cos(2 * pi * (0:n-1) / (n - 1));
end

function r = acorr_biased(x, maxLag)
    % Biased autocorrelation r(k+1) = (1/N) sum_n x(n) x(n+k), k = 0..maxLag.
    n = numel(x);
    r = zeros(1, maxLag + 1);
    for k = 0:maxLag
        r(k + 1) = sum(x(1:n-k) .* x(1+k:n)) / n;
    end
end

function [a, err] = levinson_durbin_cw(r, p)
    % Levinson-Durbin recursion. a(1) = 1, so A(z) = 1 + a(2) z^-1 + ...
    a = zeros(1, p + 1);
    a(1) = 1;
    err = r(1);
    for m = 1:p
        acc = r(m + 1) + sum(a(2:m) .* r(m:-1:2));
        k = -acc / err;
        aNew = a;
        aNew(2:m) = a(2:m) + k * a(m:-1:2);
        aNew(m + 1) = k;
        a = aNew;
        err = err * (1 - k * k);
    end
end

function [a, G, predGainDb] = lpc_autocorr(x, p, useWindow)
    % Autocorrelation-method LPC on one segment (Hamming window if useWindow).
    x = x(:).';
    if useWindow, x = x .* hamming_win(numel(x)); end
    r = acorr_biased(x, p);
    [a, err] = levinson_durbin_cw(r, p);
    G = sqrt(err);
    predGainDb = 10 * log10(r(1) / err);
end

function [a, G] = lpc_covariance(x, p)
    % Covariance-method LPC: least squares on the interior samples, no window.
    x = x(:).';
    n = numel(x);
    M = zeros(n - p, p);
    for k = 1:p
        M(:, k) = x(p - k + 1 : n - k).';
    end
    target = x(p + 1 : n).';
    coef = M \ (-target);
    a = [1, coef.'];
    resid = target + M * coef;
    G = sqrt(sum(resid .^ 2) / (n - p));
end

function [f, db] = lpc_response_db(a, G, fs, nfft)
    % Frequency response G / A(e^jw) in dB from 0 to fs/2.
    A = fft(a, nfft);
    H = G ./ A(1:nfft/2 + 1);
    f = (0:nfft/2) * fs / nfft;
    db = 20 * log10(abs(H) + 1e-12);
end

function [f, db] = segment_spectrum_db(x, fs, nfft)
    % Amplitude spectrum of the windowed segment, scaled to compare with G/|A|.
    xw = x .* hamming_win(numel(x));
    S = fft(xw, nfft) / sqrt(numel(x));
    f = (0:nfft/2) * fs / nfft;
    db = 20 * log10(abs(S(1:nfft/2 + 1)) + 1e-12);
end

function s = fSeg_spec(x, fs, nfft)
    % Bundle the segment spectrum for the plotting function.
    [f, db] = segment_spectrum_db(x, fs, nfft);
    s.f = f; s.db = db;
end

function [freqs, bws] = formants_from_lpc(a, fs)
    % First three formants from the roots of A(z): F > 150 Hz and B < 500 Hz.
    rts = roots(a);
    rts = rts(imag(rts) > 0);
    fr = angle(rts) * fs / (2 * pi);
    bw = -log(abs(rts)) * fs / pi;
    [fr, idx] = sort(fr);
    bw = bw(idx);
    keep = fr > 150 & bw < 500;
    fr = fr(keep).';
    bw = bw(keep).';
    freqs = nan(1, 3);
    bws = nan(1, 3);
    m = min(3, numel(fr));
    freqs(1:m) = fr(1:m);
    bws(1:m) = bw(1:m);
end

function pk = envelope_peaks(f, d, lo, hi)
    % Frequencies of the local maxima of a dB curve between lo and hi Hz.
    sel = f >= lo & f <= hi;
    ff = f(sel); dd = d(sel);
    idx = find(dd(2:end-1) > dd(1:end-2) & dd(2:end-1) >= dd(3:end)) + 1;
    pk = ff(idx);
end

function pk = cepstral_envelope_peaks(x, fs, K, nfft)
    % Peaks of the log spectrum smoothed by keeping the first K cepstrum coefficients.
    spec = abs(fft((x - mean(x)) .* hamming_win(numel(x)), nfft)) + 1e-12;
    c = real(ifft(log(spec)));
    lifter = zeros(1, nfft);
    lifter(1 : K + 1) = 1;
    lifter(nfft - K + 1 : nfft) = 1;
    env = real(fft(c .* lifter));
    env = env(1 : nfft/2 + 1);
    f = (0:nfft/2) * fs / nfft;
    pk = envelope_peaks(f, env, 150, 5000);
end

function s0 = choose_segment(x, fs, n)
    % Start sample (0-based) of the best window whose centre lies in the middle
    % 30% to 70% of the word. Score = periodicity (highest normalised
    % autocorrelation peak for 70-400 Hz) minus |ln(RMS first half / RMS second half)|.
    N = numel(x);
    lo = max(floor(0.3 * N - n / 2), 0);
    hi = min(floor(0.7 * N - n / 2), N - n);
    best = -Inf;
    s0 = lo;
    for s = lo : floor(0.005 * fs) : (hi - 1)
        w = x(s+1 : s+n) - mean(x(s+1 : s+n));
        r = acorr_biased(w, floor(fs / 70) + 1);
        r = r / r(1);
        period = max(r(floor(fs / 400) + 1 : end));
        steady = abs(0.5 * log(mean(w(1:floor(n/2)).^2) / mean(w(floor(n/2)+1:end).^2)));
        if period - steady > best
            best = period - steady;
            s0 = s;
        end
    end
end

function [seg, s] = segment_at(x, centre, n)
    % Segment of n samples centred on `centre`, shifted to stay inside x.
    s = min(max(centre - floor(n / 2), 0), numel(x) - n);
    seg = x(s + 1 : s + n);
end

function y = pre_emphasis(x, coef)
    y = [x(1), x(2:end) - coef * x(1:end-1)];
end

function kRef = parabolic_peak(r, k)
    % Refine an integer lag k (lag k is r(k+1) in Matlab indexing).
    y0 = r(k); y1 = r(k + 1); y2 = r(k + 2);
    d = y0 - 2 * y1 + y2;
    if d == 0, kRef = k; else, kRef = k + 0.5 * (y0 - y2) / d; end
end

function [f0, peak] = f0_autocorr_frame(frame, fs)
    % F0 from the highest autocorrelation peak between 70 and 400 Hz.
    frame = frame - mean(frame);
    n = numel(frame);
    lo = floor(fs / 400);
    hi = min(floor(fs / 70) + 1, n - 2);
    r = acorr_biased(frame, hi + 1);
    r = r / r(1);
    [peak, i] = max(r(lo + 1 : hi + 1));
    k = lo + i - 1;
    f0 = fs / parabolic_peak(r, k);
end

function [t, f0, pk] = f0_track(x, fs)
    % Frame-wise F0 with 40 ms frames and a 10 ms hop.
    n = floor(fs * 0.040);
    hop = floor(fs * 0.010);
    starts = 0:hop:(numel(x) - n);
    t = zeros(size(starts)); f0 = t; pk = t;
    for i = 1:numel(starts)
        s = starts(i);
        [f0(i), pk(i)] = f0_autocorr_frame(x(s+1 : s+n), fs);
        t(i) = (s + n / 2) / fs;
    end
end

function f0 = f0_cepstrum(x, fs)
    % Second F0 estimate from the real cepstrum of the whole segment.
    spec = abs(fft((x - mean(x)) .* hamming_win(numel(x)), 16384)) + 1e-12;
    c = real(ifft(log(spec)));
    lo = floor(fs / 400);
    hi = floor(fs / 70);
    [~, i] = max(c(lo + 1 : hi + 1));
    k = lo + i - 1;
    f0 = fs / parabolic_peak(c, k);
end

function x = impulse_train(f0, fs, dur)
    % Unit, one-sample impulses at round(k fs / f0), so the mean rate is f0.
    n = round(dur * fs);
    k = 0:ceil(n * f0 / fs + 1) - 1;
    pos = round(k * fs / f0);
    pos = pos(pos < n);
    x = zeros(1, n);
    x(pos + 1) = 1;
end

function y = resample_24k_to_10k(x)
    % Polyphase resampling by 5/12 with a Kaiser-windowed (beta 5) low-pass
    % FIR, written out so that no toolbox is needed.
    up = 5; down = 12; mx = max(up, down);
    L = 2 * 10 * mx + 1;
    n = (0:L-1) - (L - 1) / 2;
    h = ones(1, L);
    nz = n ~= 0;
    h(nz) = sin(pi * n(nz) / mx) ./ (pi * n(nz) / mx);
    beta = 5;
    t = 2 * (0:L-1) / (L - 1) - 1;
    kw = besseli(0, beta * sqrt(1 - t.^2)) / besseli(0, beta);
    h = h .* kw;
    h = h / sum(h) * up;
    xu = zeros(1, numel(x) * up);
    xu(1:up:end) = x;
    yy = conv(xu, h);
    delay = (L - 1) / 2;
    yy = yy(delay + 1 : delay + numel(x) * up);
    y = yy(1:down:end);
end

function nPass = run_self_tests(fid)
    % Six checks of the functions above against cases with a known answer.
    try, randn('state', 0); catch, end %#ok<RAND>
    fsT = 24000;
    ok = false(1, 6);
    names = {'Levinson-Durbin matches a direct solve', 'recovers a known AR(2) filter', ...
             'covariance and autocorrelation agree on a long signal', 'formants of cascaded resonators', ...
             'F0 of a pulse-driven vowel', 'impulse train rate and shape'};

    x = filter(1, [1 -0.9 0.4], randn(1, 4000));
    r = acorr_biased(x, 6);
    [a, err] = levinson_durbin_cw(r, 6);
    direct = -(toeplitz(r(1:6)) \ r(2:7).');
    ok(1) = max(abs(a(2:end) - direct.')) < 1e-9 && abs(err - (r(1) + sum(direct.' .* r(2:7)))) < 1e-9;

    x = filter(1, [1 -0.9 0.4], randn(1, 60000));
    [a, g] = lpc_autocorr(x, 2, false);
    ok(2) = max(abs(a - [1 -0.9 0.4])) < 0.02 && abs(g - 1) < 0.03;

    x = filter(1, [1 -1.2 0.8], randn(1, 30000));
    a1 = lpc_autocorr(x, 2, false);
    a2 = lpc_covariance(x, 2);
    ok(3) = max(abs(a1 - a2)) < 0.02;

    aR = conv(conv(resonator(500, 60, fsT), resonator(1500, 90, fsT)), resonator(2500, 120, fsT));
    [fT, bT] = formants_from_lpc(aR, fsT);
    ok(4) = max(abs(fT - [500 1500 2500])) < 1 && max(abs(bT - [60 90 120])) < 1;

    aR = conv(resonator(600, 80, fsT), resonator(1800, 100, fsT));
    y = filter(1, aR, impulse_train(130, fsT, 0.4));
    [~, f0T, pkT] = f0_track(y, fsT);
    ok(5) = abs(mean(f0T(pkT >= 0.5)) - 130) < 2 && abs(f0_cepstrum(y, fsT) - 130) < 3;

    x = impulse_train(133.3, fsT, 1.0);
    ok(6) = numel(x) == fsT && all(x == 0 | x == 1) && abs(sum(x) - 133.3) <= 1.5;

    for i = 1:6
        if ok(i), tag = 'PASS'; else, tag = 'FAIL'; end
        say(fid, '  test %d %s: %s\n', i, tag, names{i});
    end
    nPass = sum(ok);
end

function a = resonator(f, bw, fs)
    r = exp(-pi * bw / fs);
    th = 2 * pi * f / fs;
    a = [1, -2 * r * cos(th), r * r];
end

function light_figure(fig)
    % Newer Matlab releases can use a dark figure theme; the report needs black on white.
    try, theme(fig, 'light'); catch, end
    set(fig, 'Color', 'w');
    set(findall(fig, 'Type', 'axes'), 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', [0.15 0.15 0.15]);
    set(findall(fig, 'Type', 'legend'), 'Color', 'w', 'TextColor', 'k', 'EdgeColor', 'k');
    for t = findall(fig, 'Type', 'text').', try, set(t, 'Color', 'k'); catch, end, end
end

function make_synthesis_figure(figDir, name, seg, ySeg, excSeg, fs, f0, spec, fRes, env, nfft)
    % Impulse train, original and synthetic waveforms, and their spectra.
    try
        ncyc = round(0.03 * fs);
        t = (0:ncyc-1) / fs * 1000;
        fig = figure('Visible', 'off', 'Position', [100 100 1300 380]);
        subplot(1, 3, 1); stem(t, excSeg(1:ncyc), 'Marker', 'none'); ylim([0 1.1]);
        xlabel('Time (ms)'); ylabel('Amplitude'); grid on; title('Impulse train (30 ms)');
        subplot(1, 3, 2); plot(t, seg(1:ncyc) / max(abs(seg(1:ncyc)))); hold on;
        plot(t, ySeg(1:ncyc) / max(abs(ySeg(1:ncyc))));
        xlabel('Time (ms)'); ylabel('Amplitude (unit peak)'); grid on; title('Waveforms (30 ms)');
        legend('Original', 'Synthesised', 'Location', 'northeast');
        [fy, dby] = segment_spectrum_db(ySeg, fs, nfft);
        subplot(1, 3, 3); plot(spec.f, spec.db, 'Color', [0.6 0.6 0.6]); hold on;
        plot(fy, dby, 'Color', [0.2 0.6 0.2]); plot(fRes, env, 'r', 'LineWidth', 1.5);
        xlim([0 6000]); ylim([-100 20]); xlabel('Frequency (Hz)'); ylabel('Level (dB)'); grid on; title('Spectra');
        legend('Original', 'Synthesised', 'LPC response', 'Location', 'northeast');
        try, sgtitle(sprintf('%s.wav: synthesis with order 26 and F0 %.1f Hz', strrep(name, '_', '\_'), f0)); catch, end
        light_figure(fig); print(fig, '-dpng', '-r150', fullfile(figDir, ['fig_synthesis_' name '.png']));
        close(fig);
    catch err
        fprintf('Synthesis figure for %s was skipped: %s\n', name, err.message);
    end
end

function make_figures(figDir, name, seg, fs, s0, segMs, finalP, fm, spec, fRes, env, orders, pgs, envOrders)
    % Figure 1/2 style spectrum plot and the order comparison plot.
    try
        fig = figure('Visible', 'off', 'Position', [100 100 1100 400]);
        for p = 1:2
            subplot(1, 2, p);
            plot(spec.f, spec.db, 'Color', [0.55 0.55 0.55]); hold on;
            plot(fRes, env, 'r', 'LineWidth', 1.6);
            yl = ylim;
            for k = 1:numel(fm)
                if ~isnan(fm(k)), plot([fm(k) fm(k)], yl, 'b--'); end
            end
            if p == 1, xlim([0 fs/2]); else, xlim([0 5000]); end
            xlabel('Frequency (Hz)'); ylabel('Level (dB)'); grid on;
            if p == 1
                title('Full band');
                legend('Speech segment spectrum', sprintf('LPC response (order %d)', finalP), 'Formant estimates', 'Location', 'northeast');
            else
                title('0 to 5 kHz');
            end
        end
        try, sgtitle(sprintf('%s.wav, %d ms segment from %.3f s', strrep(name, '_', '\_'), segMs, s0/fs)); catch, end
        light_figure(fig); print(fig, '-dpng', '-r150', fullfile(figDir, ['fig_spectrum_' name '.png']));
        close(fig);

        fig = figure('Visible', 'off', 'Position', [100 100 1100 400]);
        subplot(1, 2, 1); plot(spec.f, spec.db, 'Color', [0.75 0.75 0.75]); hold on;
        leg = {'Segment'};
        cols = [0 0.447 0.741; 0.850 0.325 0.098; 0.466 0.674 0.188; 0.494 0.184 0.556];
        j = 0;
        for i = 1:numel(orders)
            if any(orders(i) == [4 12 26 40])
                j = j + 1;
                plot(fRes, envOrders{i}, 'Color', cols(j, :));
                leg{end+1} = sprintf('order %d', orders(i)); %#ok<AGROW>
            end
        end
        xlim([0 6000]); xlabel('Frequency (Hz)'); ylabel('Level (dB)'); grid on;
        title('LPC response against model order'); legend(leg, 'Location', 'northeast');
        subplot(1, 2, 2); plot(orders, pgs, 'o-'); grid on;
        xlabel('LPC order'); ylabel('Prediction gain (dB)'); title('Prediction gain against order');
        try, sgtitle(sprintf('%s.wav, %d ms segment', strrep(name, '_', '\_'), segMs)); catch, end
        light_figure(fig); print(fig, '-dpng', '-r150', fullfile(figDir, ['fig_orders_' name '.png']));
        close(fig);
    catch err
        fprintf('Figures for %s were skipped: %s\n', name, err.message);
    end
end
