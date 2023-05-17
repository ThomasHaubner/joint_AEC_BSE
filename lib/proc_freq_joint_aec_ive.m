function [Sf_hat, Ef] = proc_freq_joint_aec_ive(Xf, Uf, cfg)
%
% Author: Thomas Haubner, LMS, 2022
%
% Processing of stft-domain data by the joint acoustic echo cancellation (AEC) and 
% blind source extraction (BSE) algorithm described in [1].
%
% [1]: T. Haubner, W. Kellermann and Z. Koldovsky, 
% "Joint Acoustic Echo Cancellation and Blind Source Extraction Based on Independent Vector Extraction," 
% 2022 International Workshop on Acoustic Signal Enhancement (IWAENC), Bamberg, Germany, 2022.
%
% Input arguments:
%       Xf: multichannel STFT-domain microphone signal (channel_dim x frame_dim x freq_dim)
%       Uf: single-channel STFT-domain loudspeaker signal (1 x frame_dim x freq_dim)
%       cfg: config struct
%
% Output arguments:
%       Sf_hat: STFT-domain single-channel target speech estimate (1 x frame_dim x freq_dim)
%       Ef: STFT-domain multichannel AEC error signal (channel_dim x frame_dim x freq_dim)
%
%
% BSD 4-Clause License
% 
% Copyright (c) 2022, Thomas Haubner, Friedrich-Alexander-Universität Erlangen-Nürnberg.
% All rights reserved.
% 

    %% preparation
    num_mics = size(Xf, 1);             % number of microphones
    num_frames = size(Xf, 2);           % number of frames
    nnrb = size(Xf, 3);                 % number of non-redundant frequency bins

    % non-linearity for source model
    switch cfg.bse.source_mod
        case 'bb'
            phi_fct = @(x)(conj(x) ./ sqrt(sum(abs(x).^2, 3)));

        case 'nb'
            phi_fct = @(x)(conj(x) ./ abs(x));

        otherwise
        error('source model is not known');

    end

    switch cfg.bse.source_mod
        case 'bb'
            phiphi_fct = @(x)((1 - .5 * abs(x).^2 ./ sum(abs(x).^2, 3)) ./ sqrt(sum(abs(x).^2, 3)));

        case 'nb'
            phiphi_fct = @(x)((1 - .5 * abs(x).^2 ./ (abs(x).^2)) ./ (abs(x)));

        otherwise
        error('source model is not known');

    end


    %% initialization
    cfg.aec.filt_len = 1;

    switch cfg.aec.init 
        case 'zeros'
            Hf_hat = zeros(num_mics, cfg.aec.filt_len, nnrb)+1j*zeros(num_mics, cfg.aec.filt_len, nnrb);
            
        otherwise
            error('Not known');

    end

    switch cfg.bse.init
        case 'refMic'
            wf_bse = zeros(num_mics, 1, nnrb)+1j*zeros(num_mics, 1, nnrb);
            wf_bse(1, 1, :) = 1;

        otherwise
            error('Not known');

    end

    
    %% optimization
    eye_tens = repmat(eye(num_mics-1), 1, 1, nnrb);

    for iter_ind = 1:(cfg.opt.num_iter+1)

        % compute error signal and error signal covariance matrix        
        Df_hat = pagemtimes(Hf_hat, Uf);
        Ef = Xf - Df_hat; 

        C_EfEf_hat = 1/num_frames*pagemtimes(Ef, 'none', Ef, 'ctranspose');

        % compute SOI steering vector via orthognal constraint
        tmp = pagemtimes(C_EfEf_hat, wf_bse);
        af_soi_hat = tmp ./ pagemtimes(wf_bse, 'ctranspose', tmp, 'none');

        % compute AEC demixing filter
        % wf_aec = -pagemtimes(pagetranspose(conj(Hf_hat)), wf_bse);

        % blocking matrix
        gammaf = af_soi_hat(1, :, :);
        gf = af_soi_hat(2:end, :, :);
        Bf = zeros(num_mics-1, num_mics, nnrb);
        Bf(:, 1, :) = gf;
        Bf(:, 2:end, :) = -gammaf.*eye_tens;

        % compute AEC BG filter
        % Hf_bg = -pagemtimes(Bf, Hf_hat);

        % compute SOI and BG signal estimates
        Sf_hat = pagemtimes(wf_bse, 'ctranspose', Ef, 'none');                      % estimated near-end speech signal
        Zf_hat = pagemtimes(Bf, Ef);                                                % estimated background signal

        % compute auxiliary variables for parameter update 
        phi = phi_fct(Sf_hat);
        nu_f = mean(Sf_hat.*phi, 2);
        if strcmp(cfg.opt.type, 'newton')
            phiphi = phiphi_fct(Sf_hat);
            rho_f =  mean(phiphi, 2);
        end
        phi = phi ./ nu_f;                                                                  % rescaling
        C_ZfZf_hat = 1/num_frames*pagemtimes(Zf_hat, 'none', Zf_hat, 'ctranspose');

        if iter_ind ~= (cfg.opt.num_iter+1)     % omit optimization step in the last iteration, i.e., num_iter+1
            if cfg.opt.joint_update
                tmp_psi_wf = pagemtimes(wf_bse, conj(phi));
                inv_C_ZfZf_hat = pageinv(C_ZfZf_hat);                               % inv_C_ZfZf_hat = fastTensorInv(C_ZfZf_hat); 
                tmp_CzzBf = pagemtimes(inv_C_ZfZf_hat, Bf);
                tmp_BfinvCzzBf = pagemtimes(Bf, 'ctranspose', tmp_CzzBf, 'none');
                tmp_BfinvCzzBfEf = pagemtimes(tmp_BfinvCzzBf, Ef);
                dHf_hat = -1/num_frames*pagemtimes((tmp_psi_wf+tmp_BfinvCzzBfEf), 'none', Uf, 'ctranspose');

                if strcmp(cfg.opt.type, 'newton')
                    wf_wfH = pagemtimes(wf_bse, 'none', wf_bse, 'ctranspose');
                    var_u_f = mean(abs(Uf).^2, 2);
                    aec_hessian = (tmp_BfinvCzzBf + rho_f .* wf_wfH ./ nu_f) .* var_u_f;
                    inv_aec_hessian = pageinv(aec_hessian);
                    dHf_hat = pagemtimes(inv_aec_hessian, dHf_hat);
    
                end

            else
                dHf_hat = -1/num_frames*pagemtimes(Ef, 'none', Uf, 'ctranspose');                   % LMS update

                if strcmp(cfg.opt.type, 'newton')
                    var_u_f = mean(abs(Uf).^2, 2);
                    dHf_hat = dHf_hat ./ var_u_f;
                end

            end

            Hf_hat = Hf_hat - cfg.opt.step_size_aec * dHf_hat;                                      % gradient descent

            % update wf_bse
            dwf_bse = 1/num_frames * pagemtimes(Ef, 'none', phi, 'transpose') - af_soi_hat;

            if strcmp(cfg.opt.type, 'newton')
                bse_hessian = -C_EfEf_hat .* conj((nu_f - rho_f) ./ nu_f);
                inv_bse_hessian = pageinv(bse_hessian);
                dwf_bse = pagemtimes(inv_bse_hessian, dwf_bse);

            end

            wf_bse = wf_bse - cfg.opt.step_size_bse * dwf_bse;                                      % gradient descent
            
            if cfg.bse.normalize
                tmp = pagemtimes(C_EfEf_hat, wf_bse);
                wf_bse = wf_bse ./ sqrt(pagemtimes(wf_bse, 'ctranspose', tmp, 'none'));             % normalization

            end

        end

    end

    % fix scale ambiguity by projecting onto first error signal
    scal_coef = mean(conj(Sf_hat) .* Ef(1, :, :)) ./ mean(abs(Sf_hat).^2, 2); 
    wf_bse = conj(scal_coef) .* wf_bse;                                         % multiply bse filter by scaling coefficient to ensure later on correct computation of signal images 

    Sf_hat = scal_coef .* Sf_hat;


end

