function Xf = my_stft(x, cfg)
%
% Author: Thomas Haubner, LMS, 2022
%
    
    if verLessThan('matlab', '9.6')
        error('stft() requires at least Matlab R2019a');
    end
    assert(mod(cfg.frame_len, 2)==0);
    Xf = stft(x, cfg.fs, 'Window', cfg.window, 'OverlapLength', cfg.overlap_len,'FFTLength', cfg.frame_len);
    Xf = Xf(round(cfg.frame_len / 2):end, :, :);
    Xf = permute(Xf, [3, 2, 1]);

end