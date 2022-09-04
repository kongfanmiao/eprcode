function eseem_fft(path, keywords, zeroPadding)

% Do Fourier transform for the Tm decay, based on path and keywords
%
% eseem_fft(path, keywords, zeroPadding)
%   zeroPadding (logical): do zero-padding or not


arguments
    path
    keywords
    zeroPadding logical = false
end

[x,y] = load_file(path, keywords);
y = real(y);
[k,c] = exponfit(x,y,1);
func = @(c1, c2, k, x) c1 + c2*exp(-x/k);
bounds = [c(1)/2,c(2)/100, 0.5/k;
          c(1)*2,c(2)*100, 2/k];
startPoint = [c(1), c(2), 1/k];
fo = fitoptions(func, ...
    "StartPoint",startPoint, ...
    "Lower",min(bounds,[],1), ...
    "Upper",max(bounds,[],1));
ft = fittype(func,"options",fo);
curve = fit(x,y,ft);
yfit = curve(x);
y = y-yfit;

eseem_fft_raw(x,y, zeroPadding);

% Create the figure title based on keywords
titleStr = ['Fourier Transform of Tm Decay: ', join(keywords, ' ')];
title(titleStr);

end
