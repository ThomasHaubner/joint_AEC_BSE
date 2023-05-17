%
% Author: Thomas Haubner, LMS, 2022
%
% Exemplary processing of a time-domain microphone and loudspeaker signal 
% by the joint acoustic echo cancellation (AEC) and blind source extraction (BSE) algorithm described in [1].
%
% [1]: T. Haubner, W. Kellermann and Z. Koldovsky, 
% "Joint Acoustic Echo Cancellation and Blind Source Extraction Based on Independent Vector Extraction," 
% 2022 International Workshop on Acoustic Signal Enhancement (IWAENC), Bamberg, Germany, 2022.
%
% BSD 4-Clause License
% 
% Copyright (c) 2022, Thomas Haubner, Friedrich-Alexander-Universität Erlangen-Nürnberg.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. All advertising materials mentioning features or use of this software must
%    display the following acknowledgement:
%        This product includes software which is described in the paper "Joint Acoustic Echo Cancellation and Blind Source Extraction Based on Independent Vector Extraction," International Workshop on Acoustic Signal Enhancement (IWAENC), Bamberg, Germany, 2022 by T. Haubner, W. Kellermann and Z. Koldovsky.
% 
% 4. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
% EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
% OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


clear all; close all; clc;
addpath('lib');

if verLessThan('matlab', '9.12')
    error('pageinv() requires at least Matlab R2022a');
end

save_wav = true;                                                              % save results as wav files in "results" folder


%% settings
% aec
cfg.aec.init        = 'zeros';                 % initialization of AEC filters

% bse
cfg.bse.source_mod = 'bb';                     % probabilistic source model ('bb': broadband, 'nb': narrowband)
cfg.bse.init       = 'refMic';                 % bse filter initialization
cfg.bse.normalize  = true;                     % normalize bse filter

% optimization
cfg.opt.num_iter = 50;                         % number of optimization iterations
cfg.opt.type = 'newton';                       % 'grad' or 'newton'
cfg.opt.step_size_aec = 1;                     % step-size for optimizing AEC filter
cfg.opt.step_size_bse = 1;                     % step-size for optimizing BSE filter
cfg.opt.joint_update = true;                   % joint update of AEC and BSE filters

% data transformation
cfg.stft.frame_len = 2048;
cfg.stft.window = hamming(cfg.stft.frame_len, 'periodic');
cfg.stft.overlap_len = round(cfg.stft.frame_len / 2);


%% load data
disp('1) Loading data')
sig = struct();
[u, cfg.fs] = audioread(fullfile('data', 'loudspeaker_signal.wav'));    % time-domain loudspeaker signal
[x, cfg.fs] = audioread(fullfile('data', 'microphone_signal.wav'));     % time-domain multichannel microphone signal
cfg.stft.fs = cfg.fs;


%% process data
disp('2) Processing data')
[s_hat, e] = proc_joint_aec_ive(x, u, cfg);


%% saving results
disp('3) Saving results')
if ~exist('results', 'dir')
    mkdir('results');
else
    delete(fullfile('results', '*.wav'));
end

if save_wav
    scal_coef = max(abs([x(1:size(e,1), 1); e(:,1); s_hat])) * 1.1;
    audiowrite(fullfile('results', 'microphone_signal_refMic.wav'), x(1:size(e,1), 1) / scal_coef, cfg.fs);
    audiowrite(fullfile('results', 'aec_error_signal_refMic.wav'), e(:, 1) / scal_coef, cfg.fs);
    audiowrite(fullfile('results','target_speech_estimate.wav'), s_hat / scal_coef, cfg.fs);

end
