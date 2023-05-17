function x = my_istft(Xf, cfg)
%
% Author: Thomas Haubner, LMS, 2022
%

    if verLessThan('matlab', '9.6')
        error('stft() requires at least Matlab R2019a');
    end
    assert(mod(cfg.frame_len, 2)==0);
    Xf = permute(Xf, [3, 2, 1]);
    Xf = cat(1, conj(Xf(end-1:-1:2, :, :)), Xf);
    x = istft(Xf, cfg.fs, 'Window', cfg.window, 'OverlapLength',cfg.overlap_len,'FFTLength', cfg.frame_len);

end