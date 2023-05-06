% This function simulates the phase shifter (delay line) of each MZI.
%
% - Lmzi = length of each MZI arm
% - fis = additional phase shift in the sup MZI arm
% - fii = additional phase shift in the inf MZI arm
% - ag = arm dB losses



function [DL]=Delay_Line_MZI(delta_L,Lmzi,freq,nr,ag,fis,conta)

c0=299792458;


DL=[ exp((ag -1i*2*pi*freq*nr/c0)*Lmzi -1i*fis)               0 ;
    
    0                  exp((ag -1i*2*pi*freq*nr/c0)*Lmzi )  ];

end


