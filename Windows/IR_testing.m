[IRR] = textread('IRR-cave.txt', '%f','delimiter', ',');
[IRL] = textread('IRR-cave.txt', '%f','delimiter', ',');

idx_peaks = find(IRR(:,1)~=0);
impulses = IRR(idx_peaks, 1)
npeaks = sum(IRL(:,1)~= 0)

%IRR(find(IRR(:,1)> 1)) = 1;
%IRL(find(IRL(:,1)>1)) = 1;

% Plot
stem(IRL)
figure

% Sound source
    [x, fs] = audioread('Piano.wav');
    yR = conv(x, IRR);
    yL = conv(x, IRL);
    
% Fast convolution
%     Ly=length(x)+length(h)-1;  % 
%     Ly2=pow2(nextpow2(Ly));    % Find smallest power of 2 that is > Ly
%     X=fft(x, Ly2);		   % Fast Fourier transform
%     H=fft(h, Ly2);	           % Fast Fourier transform
%     Y=X.*H;        	           % 
%     y=real(ifft(Y, Ly2));      % Inverse fast Fourier transform
%     y=y(1:1:Ly);               % Take just the first N elements
%     y=y/max(abs(y)); 

minL = min([length(yR) length(yL)]);

plot(x); hold on;
plot(yL,'r')

soundsc([yL(1:minL) yR(1:minL) ], fs)