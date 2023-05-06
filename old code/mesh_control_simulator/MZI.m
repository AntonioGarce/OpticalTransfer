% This function produces the transfer matrix of each step.
% One step is made by a phase shifter (lower heater) as a Delay Line + MZI (coupler+phase
% shifter+coupler).
% Hence, the whole step matrix is the product of each single component.



function [Matrix] = MZI(k,ag,f,nr,Lmzi,delta_L,fis,phasei,L,conta)

 % coupler transfer matrix
 
 Tc = [     sqrt((1-k))              -1i*(sqrt(k))  ;
        -1i*(sqrt(k))             sqrt(1-k)         ];
 
    
 % delay line transfer matrix

 DL = Delay_Line_MZI(delta_L,Lmzi,f,nr,ag,fis,conta);
 
 
 prod1 = Tc*DL;
 
 prod2 = prod1*Tc;
 
 [Matrix] = prod2*Delay_Line2(L,f,nr,ag,phasei); % consider also the heater at the MZI input
 
end