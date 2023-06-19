%% ************************************************************************ 
%
% Compute a Butterworth band-pass filter in frequency
%
% [stfinal]=computebandfftfilter(signal,dt,fcU,p,fcL)
%
%         signal - signal to filter
%             dt - time step
%            fcU - upper corner
%            fcL - lower corner (if low-pass choose -fcU)
%              p - number of poles
%  
%    It will add the same amount of 0s as the size of the signal before
%    filtering.
%    The lower part of the filter will have 2*p
%
%% ************************************************************************

function [stfinal]=computebandfftfilter(signal,dt,fcU,p,fcL);

signal=[signal' ]';%signal'*0]';
timeSteps = size(signal);
if(timeSteps(1)==1)
    numsteps=timeSteps(2);
end
if(timeSteps(2)==1)
    numsteps=timeSteps(1);
end
NFFT = 2^nextpow2(length(signal));
fourier=fft(signal,NFFT);

if(fcL <=0)
    fcentral=.5*(fcU+fcL);
    fc=fcU-fcentral;
    f   = (1/dt)*( 0:(NFFT-1) )/NFFT;
    fil =  1./(1.+(abs(f-fcentral)./fc).^(2*p));
else
    fcentral=.5*(fcU+fcL);
    fc=fcU-fcentral;
    f   = (1/dt)*( 0:(NFFT-1) )/NFFT;
    k=find(f >fcentral,1);
    
    fil2 =  1./(1.+(abs(f-fcentral)./fc).^(2*p));   % for upper half of the filter    
    fil1 =  1./(1.+(abs(f-fcentral)./fc).^(2*4*p)); % for lower half of the filter
    kl=find(fil1 > .99,1);
    ones=fil1*0+1;
    ramp=[linspace(0,1,kl)  ones(kl+1:length(ones))];
    fil=[fil1(1:k) fil2(1+k:length(fil2))];
    
end
signfiltf = fil'.*fourier;

% complete the congugate
for ij=1:NFFT/2-1
   signfiltf(NFFT-ij+1)=conj(signfiltf(ij+1));     
end
    
sft = ifft(signfiltf);

stfinal=sft(1:length(signal));%/2);

return;
