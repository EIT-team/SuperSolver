function [Zc] = contact_impedance(Data_comp_fn, freqs,Fem)
% Usage: [Zc] = contact_impedance(Data_comp_fn, freqs)
%
% General: estimate contact impedance out of the 2 terminal measurements data file. this function is temporarly, until contact impedance recovary would be embedded in the reconstruction procedure
%
% Input:
% Data_comp_fn (global) - compaliance data file name {string}
% Output:
% Zc (global) - contact impedance estimate {no of electrodes};
%
% LH 15/12/05
%==============================================================================%

% in case there is no compliance data file and the contact impedance is given as a scalar
if isempty(Data_comp_fn)
    if (length(Fem.zc)<2)
        Zc = Fem.zc * ones(size(Fem.pos,1),1);
    end
    return
end