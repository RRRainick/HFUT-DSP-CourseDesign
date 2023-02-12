function [Xk] = mydft(x, N)
    WN = exp(-1j*2*pi/N);
    x_m = x(1: N);
    k = (0: N-1); n = (0: N-1);
    %% ak = 1/N*WN.^(k.*n')*x_m';
    Xk = x_m * (WN.^(k'.*n));
end
