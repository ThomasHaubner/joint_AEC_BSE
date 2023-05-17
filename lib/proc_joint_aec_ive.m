function [s_hat, e] = proc_joint_aec_ive(x, u, cfg)
%
% Author: Thomas Haubner, LMS, 2022
%
% Processing of time-domain data by the joint acoustic echo cancellation (AEC) and 
% blind source extraction (BSE) algorithm described in [1]..
%
% [1]: T. Haubner, W. Kellermann and Z. Koldovsky, 
% "Joint Acoustic Echo Cancellation and Blind Source Extraction Based on Independent Vector Extraction," 
% 2022 International Workshop on Acoustic Signal Enhancement (IWAENC), Bamberg, Germany, 2022.
%
% Input arguments:
%       x: multichannel microphone signal (time_dim x channel_dim)
%       u: single-channel loudspeaker signal (time_dim x 1)
%       cfg: config struct
%
% Output arguments:
%       s_hat: time-domain single-channel target speech estimate (time_dim x 1)
%       e: time-domain multichannel AEC error signal (time_dim x channel_dim)
%
%
% BSD 4-Clause License
% 
% Copyright (c) 2022, Thomas Haubner, Friedrich-Alexander-Universität Erlangen-Nürnberg.
% All rights reserved.
% 


    %% transforming data to STFT domain
    Uf = my_stft(u, cfg.stft);
    Xf = my_stft(x, cfg.stft);


    %% estimating aec+bse filters
    [Sf_hat, Ef] = proc_freq_joint_aec_ive(Xf, Uf, cfg);

    
    %% inverse time-domain transformation
    e = my_istft(Ef, cfg.stft);                % time-domain AEC error signal
    s_hat = my_istft(Sf_hat, cfg.stft);        % time-domain target speech estimate from BSE


end