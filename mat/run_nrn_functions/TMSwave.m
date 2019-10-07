% Outputs normalized TMS waveform E and time vector t
% Input waveform data structure containing:
% frequency f (kHz), mode (1 or 2 for monophasic or biphasic, respectively, 3 for recorded biphasic),
% tau, for damping constant (ms), dt for time step (ms), and del for delay (ms)
% Ex:
% wave.mode = 1; 
% wave.dur = 30; 
% wave.del = 0;
% wave.dt = 0.001;
% [t,E] = TMSwave(wave);
function [tvec,Evec] = TMSwave(dt,tstop,waveform,plot_wave)
    if nargin == 0                
        % modes
        % 1-3 = Koponen recorded MagProX100 waveforms (1. mono,2. bi,3. half)        
        % 4 - cTMS1 waveforms (30 µs - 160 µs)        
        mode = 1; 
        dt = 0.005;
        tstop = 1; % ms
        del = 0.005; % ms
        dur = 0.03; % ms - only affects cTMS1 pulses (mode = 4)
        plot_wave = 1; % 1 to plot waveform
    else
        dur = waveform.dur;
        mode = waveform.mode;        
        del = waveform.del;
        %amp = waveform.amp;
    end          
    switch mode        
        case 1 % MagproX100 Monophasic waveform - Recorded by Lari Koponen (6/12/18)
            % Get data
            data = load('TMSwaves.mat');
            trec = data.tm; Erec = data.Erec_m;                        
            title_string = sprintf('Normalized E - MagproX100 Monophasic, del = %.1f ms', del);            
        case 2 % MagproX100 Biphasic waveform - Recorded by Lari Koponen (6/12/18)
            % Get data
            data = load('TMSwaves.mat');
            trec = data.tb; Erec = data.Erec_b;            
            title_string = sprintf('Normalized E - MagproX100 Biphasic, del = %.1f ms', del);            
        case 3 % MagproX100 Half-sine waveform - Recorded by Lari Koponen (6/12/18)
            % Get data
            data = load('TMSwaves.mat');
            trec = data.th; Erec = data.Erec_h;            
            title_string = sprintf('Normalized E - MagproX100 Half-sine, del = %.1f ms', del);                    
        case 4 % cTMS1 waveforms - from ctms1_wvfrm_11_26_2008.mat used in Peterchev 2013
             % Get data
            data = load('cTMSwaves.mat');
            trec = data.t; 
            Erec_all = data.W;
            pw = data.pw;             
            Eind = value2ind(dur,pw); % index of waveform to use
            Erec = Erec_all(:,Eind);
            if dur ~= pw
               fprintf('Using closest recorded cTMS waveform: %.1f us instead of input: %.1f us\n',pw(Eind)*1e3,dur*1e3)
            end                       
            title_string = sprintf('Normalized E - Recorded cTMS wave, PW = %g us, del = %.1f ms', pw(Eind)*1000, del);                    
    end
    % Downsample for simulation specific time step
    DT = mean(diff(trec));
    sample_factor = int8(dt/DT);
    if sample_factor < 1 % dt chosen is smaller than DT of recording
       sample_factor = 1; % don't downsample 
       fprintf('dt chosen smaller than recording, set to %f ms\n',DT); 
    end
    % construct final vectors            
    if del >= dt
        tvec = [0;del-dt;downsample(trec,sample_factor)+del]; % include delay           
        Evec = [0;0;downsample(Erec,sample_factor);0];
    else
        tvec = downsample(trec,sample_factor); % pulse starts at t = 0
        Evec = [downsample(Erec,sample_factor);0];                
    end                      
    tvec = [tvec;(tvec(end)+dt:dt:tstop)']; % at end point to end at 0            
%             Evec = Evec/max(Evec); % renormalize to max after downsampling            
    if length(tvec) > length(Evec)
        Evec = [Evec;zeros(length(tvec)-length(Evec),1)]; % add trailing zeros
    else
        Evec = Evec(1:length(tvec));
    end
    if plot_wave
        figure('Color','w');
        hold on;
        plot(tvec,Evec,'Color','k','LineWidth',2)
        title(title_string);
        xlabel('time (ms)'); xlim([0, tstop]);
        box off;
        set(gca,'FontSize',20)
    end
end