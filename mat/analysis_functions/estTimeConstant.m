function [tau,rb,resnorm] = estTimeConstant(pws,mt,t,wvfrms,varargin)
%ESTTIMECONSTANT Estimate time constant using method of Peterchev et al.
%2013 (lsqcurvefit)
%  
%   Inputs 
%   ------ 
%   pws : vector
%         pulse-widths in ms, same order as waveform recordings in wvfrms
%   mt : vector
%        threshold values, same size as pws
%   t : vector
%       time vector in seconds
%   wvfrms : array
%            num_timepoints x num_pws array of waveform recordings 
%   Optional Inputs 
%   --------------- 
%   tau_m0 : double
%            starting point for time constant estimation in sec
%   tau_m_lb : double
%              lower bound for time constant estimation in sec
%   tau_m_ub : double
%              upper bound for time constant estimation in sec
%   mt_b0 : double
%           starting point for estimation of rheobase (same units as mt)
%   mt_b_lb : double
%             lower bound for rheobase estimation
%   mt_b_ub : double
%             upper bound for rheobase estimation
%   Outputs 
%   ------- 
%   tau : double
%         time constant estimation in µs
%   rb : double
%        rheobase amplitude estimation in same units as mt
%   resnorm : double
%             squared norm of residual, output of lsqcurvefit
%   Examples 
%   --------------- 

% AUTHOR    : Aman Aberra 

% For lsqcurve fitting - input parameters
in.tau_m0 = 200e-6; % sec
in.tau_m_lb = 10e-6; %sec
in.tau_m_ub = 5e-3;
in.mt_b0 = 30;
in.mt_b_lb = 1; 
in.mt_b_ub = 5000; 
in = sl.in.processVarargin(in,varargin); 
options = optimoptions(@lsqcurvefit,'MaxFunEvals', 1000,'Algorithm','trust-region-reflective','Display','off');
[param,resnorm,~,~,~,~,~] = lsqcurvefit(@mt_mtcalc,...
                [in.tau_m0,in.mt_b0],pws,mt./mt,[in.tau_m_lb,in.mt_b_lb],...
                [in.tau_m_ub,in.mt_b_ub],options);
% [param,resnorm,~,~,~,~,~] = lsqcurvefit(@mt_mtcalc,...
%                 [in.tau_m0,in.mt_b0],pws,mt,[in.tau_m_lb,in.mt_b_lb],...
%                 [in.tau_m_ub,in.mt_b_ub],options);
tau = param(1)*1e6; % time constant in us
rb = param(2); % rheobase
       
    function mt_mt = mt_mtcalc(param,pws)
        % param(1) = membrane time const. (tau_m)
        % param(2) = depolarization threshold volt. relative to pulse amplitude (%)
        %           (rheobase)
        tau_m = param(1);   % membrane time const.
        mt_b = param(2);    % depolarization threshold volt. relative to pulse amplitude (%)
        n_pws = length(pws);
        mt_mt = zeros(n_pws,1);
        fs = 1e3/(t(2)-t(1));  % assumes waveform sampling frequency same for all waveforms
        b = 1/(1+2*tau_m*fs).*[1 1];
        a = [1 (1-2*tau_m*fs)/(1+2*tau_m*fs)];        
        for i = 1:n_pws            
            mt_mt(i) = mt_b/mt(i)/max(filter(b,a,wvfrms(:,i)));
%             mt_mt(i) = mt_b/max(filter(b,a,wvfrms(:,i)));
        end
    end

end