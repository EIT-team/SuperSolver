function [Zc] = contact_impedance(Data_comp_fn, freqs);
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
global Fem Data

% in case there is no compliance data file and the contact impedance is given as a scalar
if isempty(Data_comp_fn)
    if (length(Fem.zc)<2)
        Zc = Fem.zc * ones(size(Fem.pos,1),1);
    end
    return
end


% fig_title='cable movement long leads 2 terminals gain 1 pure resistor 100 ohm';
% f1=figure('Name',fig_title);
% mesh(abs(squeeze(double(DATA))))
% valid_samples=[1:175];
% %semilogx(freq, abs(mean(squeeze(double(DATA(:,:,:,valid_samples))),2)),'o-')
% semilogx(freq, abs(mean(squeeze(double(DATA(:,:,:,valid_samples))),2)),freq, abs(squeeze(mean(double(t2_comp(:,:,5)),2))),'o-')
% axis tight
% ax=axis;
% axis([ax(1:2) 350 750])
% title(fig_title, 'FontName', 'Comic Sans MS','FontWeight','bold','FontSize',10)
% xlabel('frequency {Hz}', 'FontName', 'Comic Sans MS')
% ylabel('voltage {DSP units}', 'FontName', 'Comic Sans MS')