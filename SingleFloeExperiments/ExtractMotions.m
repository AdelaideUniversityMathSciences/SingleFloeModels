% function ExtractMotions(ice,motion,runnos,do_plot)
%
% LJ YIEW
% Created on  Nov 2015
% Last edited Oct 2016
%
% Extracts the 6 degree-of-freedom motions for
%  Floe NB (floe without a barrier)
%  Floe B  (floe with a barrier)
%
% INPUTS:
%  ice:
%   1: Floe NB
%   2: Floe B
%  motion  : 'Surge' or 'Heave' or 'Pitch'
%  runnos  : run numbers
%  do_plot : flag for time series plots
%   1: plot
%   0: don't plot 
%
% FILES NEEDED:
%  AMC_InputData.xlsx
%  AMC_DataFilter.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExtractMotions(ice,motion,runnos,do_plot)

% clear all
close all
% clc

%% INPUTS

% SELECT WHICH FLOE (1:NO BARRIER, 2:BARRIER)
if ~exist('ice','var')
 ice = 2;
end

% SELECT TYPE OF MOTION TO EXTRACT
if ~exist('motion','var')
 motion = 'Surge'; % e.g. 'Surge'
end

% SELECT CONDITION NUMBER
% condition = 1;
% (OR)
% SELECT RUN NUMBERS
if ~exist('runnos','var')
 runnos = 62;%[6:35,37:39,41:62,67:71,74,77:83];
end


% SELECT PLOT/NOPLOT
if ~exist('do_plot','var')
 do_plot = 1; % 1:plot 0:don't plot
end

% USABLE RUNS:
% condition 1: waveheight = 0.02 m, freq = [0.5,2] Hz
% [6;7;8;9;10;11;12;13;14;33;34;35;37;38;67;68;69;70;71]'
% condition 2: waveheight = 0.04 m, freq = [0.5,2] Hz
% [15;16;17;18;19;20;21;22;23;24;39;41;42;43;44;74;77;78;79;80;81;82]'
% condition 3: waveheight = 0.08 m, freq = [0.5,1.5] Hz
% [25;26;27;28;29;30;31;32;83]'
% condition 4: waveheight = [0.01,0.08] m, freq = 1.5 Hz
% [45;46;47;48;49;50]'
% condition 5: waveheight = [0.01,0.1] m, freq = 1.25 Hz
% [51;52;53;54;55;56;57]'
% condition 6: waveheight = [0.005,0.04] m, freq = 1.8 Hz
% [58;59;60;61;62]'

%%

% READ MEASURED DATA FOR:
%  FREQUENCY, WAVE HEIGHT, WAVELENGTH, WAVE CELERITY (f*lambda)
run_dat = xlsread('AMC_InputData','A3:F85');

% READ STEADY STATE DATA RANGE (START/END TIME):
ss_dat = xlsread('AMC_DataFilter','A3:G85');

%%

if exist('condition')
 switch condition
  case 1
   runnos = [6;7;8;9;10;11;12;13;14;33;34;35;37;38;67;68;69;70;71]';
  case 2
   runnos = [15;16;17;18;19;20;21;22;23;24;39;41;42;43;44;74;77;78;79;80;81;82]';
  case 3
   runnos = [25;26;27;28;29;30;31;32;83]';
  case 4
   runnos = [45;46;47;48;49;50]';
  case 5
   runnos = [51;52;53;54;55;56;57]';
  case 6
   runnos = [58;59;60;61;62]';
 end
end

if ice == 1
 switch motion
  case 'Surge'
   M_no = 3;
   units = 'mm';
  case 'Heave'
   M_no = 5;
   units = 'mm';
  case 'Sway'
   M_no = 4;
   units = 'mm';
  case 'Roll'
   M_no = 6;
   units = 'deg';
  case 'Pitch'
   M_no = 7;
   units = 'deg';
  case 'Yaw'
   M_no = 8;
   units = 'deg';
 end
 barr = 'NoBarrier';
 
elseif ice == 2
 switch motion
  case 'Surge'
   M_no = 19;
   units = 'mm';
  case 'Heave'
   M_no = 21;
   units = 'mm';
  case 'Sway'
   M_no = 20;
   units = 'mm';
  case 'Roll'
   M_no = 22;
   units = 'deg';
  case 'Pitch'
   M_no = 23;
   units = 'deg';
  case 'Yaw'
   M_no = 24;
   units = 'deg';
 end
 barr = 'Barrier';
end

%%

c1 = 1; % counter for run number

for run = runnos
 
 run

 clear M_max M_min time_M_max time_M_min

 % ADJUST FORMAT FOR FILE NAME
 if run < 10
  runno = ['00',num2str(run)];
 elseif run > 99
  runno = num2str(run);
 else
  runno = ['0',num2str(run)];
 end

 % READ QUALISYS DATA
 filename = ['RawData/Qualisys/Run0',runno,'_6D.tsv'];
 qualysis = fopen(filename,'rt');
 data     = textscan(qualysis,'%f','HeaderLines',11);
 data     = data{1};
 fclose(qualysis);

 % RUN PARAMETERS
 lambda   = run_dat(str2num(runno),5);     % wavelength
 c        = run_dat(str2num(runno),6);     % wave celerity
 freq     = run_dat(str2num(runno),2);     % frequency
 w_height = run_dat(str2num(runno),3)*1e3; % wave height

 % ARRANGE AND EXTRACT DATA
 n = 34;               % corresponds to no. of columns in .tsv file
 time = data(2:n:end); % extract time
 % extract M position, normalise wrt initial position
 M = data(M_no:n:end)-data(M_no); % M_no 22-roll, 20-sway, 24-yaw

 % CALCULATE STEADY STATE PERIOD
 % approx dist between wave maker and floe (MTB length = 35m)
 dist = 20;
 ss_start = dist/c;          % steady state start time
 ss_end   = 35+(35-dist)/c;  % steady state end time

 %%
 if do_plot == 1
  % PLOT RAW DATA VS TIME
  figure
  set(gcf,'position',[100 100 1000 600]);
  set(gca,'FontSize',14)
  hold on
  plot(time,M)
  ylabel([motion,' [',units,']']),xlabel('Time [s]')
  title(['Run ',num2str(run),', Frequency = ',num2str(freq),...
         ' Hz, Wave Height = ',num2str(w_height),' ',units])     
 end

 %%
 % SELECT START/END TIME FOR DATA PROCESSING
 starttime = ss_dat(run,2)
 endtime   = ss_dat(run,3)
%  starttime = input('start time: '); % start time in seconds
%  endtime   = input('end time:   '); % end time
 if endtime > 60
  endtime = 60;
 end
 % counter for start and end time i.e. corresponding to frame no.
 t_start = round(length(time)/time(end)*starttime);
 t_end   = round(length(time)/time(end)*endtime);

 %% 
 % SMOOTHEN DATA (FROM START TO END TIME)
 s_factor = 0.01; % default
 M_smooth = smooth(time(t_start:t_end),M(t_start:t_end),s_factor,'lowess');
 
 if do_plot == 1
  % PLOT SMOOTH DATA
  plot(time(t_start:t_end),M_smooth,'g')
 end
 
 
 %%
 if strcmp(motion,'Surge')
  
  % FIND GRADIENT (DRIFT)
  c_factor = ss_dat(run,6)
  M_grad   = smooth(linspace(time(t_start),time(t_end),...
              length(M_smooth)),M_smooth,...
              c_factor,'lowess'); % <== adjust smoothing factor (if necessary)

  % REMOVE DRIFT (SURGE ONLY)
  M_norm = M_smooth - M_grad;   

  % CALCULATE GRADIENT (DRIFT VELOCITY)
  grad_M(c1) = (M_grad(1) - M_grad(end))/(endtime - starttime)
  
  if do_plot == 1
   % PLOT DRIFT, SURGE
   plot(linspace(time(t_start),time(t_end),length(M_grad)),M_grad,'r')
   plot(linspace(time(t_start),time(t_end),length(M_norm)),M_norm,'c')
  end
  
 else
  M_norm = M_smooth;
  
 end

 %%
 % FIND LOCAL MAX AND MIN FOR HEAVE
 c2 = 2; % counter for ss range
 d1 = 1; % counter for max surge
 d2 = 1; % counter for min surge
 for time_c = [t_start:t_end-2];
  % find local max
  if M_norm(c2) > M_norm(c2-1) && M_norm(c2) > M_norm(c2+1)
   M_max(d1) = M_norm(c2);
   time_M_max(d1) = time(c2-1+t_start);
   d1 = d1+1;
  % find local min
  elseif M_norm(c2) < M_norm(c2-1) && M_norm(c2) < M_norm(c2+1)
   M_min(d2) = M_norm(c2);
   time_M_min(d2) = time(c2-1+t_start);
   d2 = d2+1;
  end
   c2 = c2+1;
 end
 
 %%
 if strcmp(motion,'Yaw')
  M_max_all = max(M_norm);
  M_min_all = min(M_norm);
  M_g_max(c1,1) = M_max_all - M_min_all;
 end
 
 %%
 if exist('M_max') && exist('M_min')
 
 % ARRANGE MAX/MIN PEAK DATA AND CORRESPONDING TIME
 M_max = [M_max', time_M_max'];
 M_min = [M_min', time_M_min']; 
 
 M_max_size(c1,1) = length(M_max);
 
 %%
 % CALCULATE AVERAGE AMPLITUDES
 M_avg_max   = mean(M_max(:,1)); % average max
 M_avg_min   = mean(M_min(:,1)); % average min
 M_avg(c1,1) = mean(M_max(:,1)) - mean(M_min(:,1));
 
 %%
 % ERROR BARS (MAX & MIN)
 M_err_max = abs(max(M_max(:,1)) - min(M_min(:,1)));
 M_err_min = abs(min(M_max(:,1)) - max(M_min(:,1)));
 M_err(c1,:) = [M_err_max M_err_min];

 %%
 
 if do_plot == 1
  % PLOT LOCAL MAX/MIN AND AVG MAX/MIN
  plot(M_max(:,2),M_max(:,1),'rv','MarkerSize',6)
  plot(M_min(:,2),M_min(:,1),'r^','MarkerSize',6)
  plot([time(t_start) time(t_end)],[M_avg_max M_avg_max],'r--')
  plot([time(t_start) time(t_end)],[M_avg_min M_avg_min],'r--') 
  xlim([0 60])
  center = (max([max(M_max(:,1)) max(M)])+min([min(M_min(:,1)) min(M)]))/2;
  half   = (max([max(M_max(:,1)) max(M)])-min([min(M_min(:,1)) min(M)]))/2;
  ylim([center-half*1.2 center+half*1.2])

  % PLOT STEADY STATE WINDOW
  plot([ss_start ss_start],[center-half*1.2 center+half*1.2],'r')
  plot([ss_end ss_end],    [center-half*1.2 center+half*1.2],'r')

  % PRINT START/END TIMES AND HEAVE
  text(1,0,{...
            ['Start time: ',num2str(starttime)],...
            ['End time: ',num2str(endtime)],...
            [motion,' amplitude (x2): ',num2str(M_avg(c1)),units],...
            [motion,' min/max: ',num2str(M_err_min),'/',num2str(M_err_max)]})
 end
          
 end

 box on
 hold off
 
%  % SAVE FIGURE
%  saveas(gca,['../../Figures/Experiments/SingleFloe/',motion,'_',barr,'/Run',runno,'.eps'],'epsc')
%  saveas(gca,['../../Figures/Experiments/SingleFloe/',motion,'_',barr,'/Run',runno,'.fig'])
%  close

 c1 = c1+1;
end

% M_avg  = M_avg'   % display average heave
% n.b. data is then copied to 'Excel/ExpData_Ice1Ice2.xlsx'

clearvars -except M_max M_min

end





