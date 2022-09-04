function eseem_fft_raw(x,y, zeroPadding)
% Do ESEEM fft for the raw data: signal ~ 2tau

arguments
    x double
    y double
    zeroPadding logical = false
end

x = x/2;
N = numel(x);
% zero padding
if zeroPadding
    N = 2^nextpow2(N);
end
yFFT = abs(fft(y,N))/N;
freqAxis = (0:N/2-1)/(N*(x(2)-x(1)));
yFFT = yFFT(1:N/2);
yFFT = yFFT/max(yFFT);
freqAxis = freqAxis*1e3; % Convert to MHz

plot(freqAxis, yFFT);
xlim([0, max(freqAxis)]);
ylim([0, 1.1]);
xlabel('Frequency (MHz)');
ylabel('Intensity (arb. u.)');

end