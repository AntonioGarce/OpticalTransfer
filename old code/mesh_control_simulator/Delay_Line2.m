% This function simulates the heater in the input lower arm of each MZI of
% the mesh.
%
% - phasei = additional phase shift of the heater



function [DL]=Delay_Line2(L,freq,nr,ag,phasei)

c0=299792458;


DL=[ 1                      0 ;
    
     0             exp(-1i*phasei) ];
 
 
%  
%  DL=[ exp((ag -1i*2*pi*freq*nr/c0)*L)                      0 ;
%     
%      0             exp((ag -1i*2*pi*freq*nr/c0)*L -1i*phasei) ];

end