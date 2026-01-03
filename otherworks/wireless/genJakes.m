% %-------------------- genJakes.m --------------------%
% Generating a single fading coefficient with the 
% `sum of sinusoids'  Jakes Model. 
% I = frame length
% v = terminal speed in km/h

function s_out=genJakes(s_in, fc, fs, v, seed)


%fc = 5.2e9; 			% Carrier frequency in Herz for HiperLan
c  = 3e8; 		   	% speed of light in meters/second

fdmax=(v*fc)/(3.6*c);   	% Maximum Doppler frequency
lambda=c/fc;			% The wavelength corresponding to fc

N=100;			    	% Number of incident waves
t=1/fs*[1:length(s_in)];      	% The time variable
% The symbol duration in HIPERLAN is 4 us,
% therefore sampling rate is 1/4us 
len=length(t);
theta=rand(1,N)*2*pi;       	% Generating the uniform phases
fd=cos(2*pi*((1:N)/N))*fdmax;   % Generating uniformly spaced 
				% frequencies from -fdmax to +fdmax
E=exp(j.*(2*pi*fd(:)*t(:)'+repmat(theta(:),1,len)));
fadingcoeff=sum(E)/sqrt(N);
s_out = abs(fadingcoeff).*s_in;
%plot(t,abs(fadingcoeff))
%xlabel('time (seconds)');ylabel('Envelope of the fading coefficient');
