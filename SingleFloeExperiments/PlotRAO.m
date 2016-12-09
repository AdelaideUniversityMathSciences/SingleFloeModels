% function PlotRAO
%
% LJ YIEW
% Created on  Nov 2015
% Last edited Oct 2016
%
% Plots RAOs for surge, heave and pitch from AMC single-floe experiments
%
% FILES NEEDED:
%  AMC_DataRAO.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PlotRAO

clear all
clc
% close all

L = 0.2;   % floe half-length [m]
h = 0.831; % water depth [m]

% READ EXPERIMENTAL DATA
dat_ice = xlsread('AMC_DataRAO','B5:AB54'); % conditons 1-3 only
% % conditions 4,5,6 only
% dat_floeNB_c4 = xlsread('Excel/ExperimentalData_Ice1','C66:AA71');
% dat_floeNB_c5 = xlsread('Excel/ExperimentalData_Ice1','C72:AA78');
% dat_floeNB_c6 = xlsread('Excel/ExperimentalData_Ice1','C79:AA83');
% dat_floeB_c4 = xlsread('Excel/ExperimentalData_Ice1','C66:AA71');
% dat_floeB_c5 = xlsread('Excel/ExperimentalData_Ice1','C72:AA78');
% dat_floeB_c6 = xlsread('Excel/ExperimentalData_Ice1','C79:AA83');

% ARRANGE DATA
wavelength  = dat_ice(:,6);       
wavenum     = 2*pi./wavelength;
waveheight  = dat_ice(:,5)*1e-3;          % meters
waveheightR = round(waveheight.*1e2)./1e2; % round waveheight to nearest 10mm

for j = 1:length(waveheightR)
 if waveheightR(j) <= 0.025
  whR(j) = 1;
 elseif waveheightR(j) > 0.025 && waveheightR(j) <= 0.04
  whR(j) = 2;
 elseif waveheightR(j) > 0.04
  whR(j) = 3;
 end
end
overwash    = dat_ice(:,18);

surge_floeNB  = dat_ice(:,9)*1e-3;    % meters
heave_floeNB  = dat_ice(:,10)*1e-3;   % meters
pitch_floeNB  = dat_ice(:,11)*pi/180; % radians
surge_err_max_floeNB = abs(dat_ice(:,12)*1e-3 - surge_floeNB);
surge_err_min_floeNB = abs(dat_ice(:,13)*1e-3 - surge_floeNB);
heave_err_max_floeNB = abs(dat_ice(:,14)*1e-3 - heave_floeNB);
heave_err_min_floeNB = abs(dat_ice(:,15)*1e-3 - heave_floeNB);
pitch_err_max_floeNB = abs(dat_ice(:,16)*pi/180 - pitch_floeNB);
pitch_err_min_floeNB = abs(dat_ice(:,17)*pi/180 - pitch_floeNB);

surge_floeB  = dat_ice(:,19)*1e-3;   % meters
heave_floeB  = dat_ice(:,20)*1e-3;   % meters
pitch_floeB  = dat_ice(:,21)*pi/180; % radians
surge_err_max_floeB = abs(dat_ice(:,22)*1e-3 - surge_floeB);
surge_err_min_floeB = abs(dat_ice(:,23)*1e-3 - surge_floeB);
heave_err_max_floeB = abs(dat_ice(:,24)*1e-3 - heave_floeB);
heave_err_min_floeB = abs(dat_ice(:,25)*1e-3 - heave_floeB);
pitch_err_max_floeB = abs(dat_ice(:,26)*pi/180 - pitch_floeB);
pitch_err_min_floeB = abs(dat_ice(:,27)*pi/180 - pitch_floeB);

% no overwash data for floe NB
c = 1;
for j = 1:length(surge_floeNB)
 if overwash(j) == 1
 else
 surge_floeNB_no(c) = surge_floeNB(j);
 heave_floeNB_no(c) = heave_floeNB(j);
 pitch_floeNB_no(c) = pitch_floeNB(j);
 surge_err_max_floeNB_no(c) = surge_err_max_floeNB(j);
 surge_err_min_floeNB_no(c) = surge_err_min_floeNB(j);
 heave_err_max_floeNB_no(c) = heave_err_max_floeNB(j);
 heave_err_min_floeNB_no(c) = heave_err_min_floeNB(j);
 pitch_err_max_floeNB_no(c) = pitch_err_max_floeNB(j);
 pitch_err_min_floeNB_no(c) = pitch_err_min_floeNB(j);
 waveheight_no(c) = waveheight(j);
 wavenum_no(c)    = wavenum(j);
 wavelength_no(c) = wavelength(j);
 c = c+1;
 end
end

% all data for floe NB
c = 1;
for j = 1:length(surge_floeNB)
 surge_floeNB_a(c) = surge_floeNB(j);
 heave_floeNB_a(c) = heave_floeNB(j);
 pitch_floeNB_a(c) = pitch_floeNB(j);
 surge_err_max_floeNB_a(c) = surge_err_max_floeNB(j);
 surge_err_min_floeNB_a(c) = surge_err_min_floeNB(j);
 heave_err_max_floeNB_a(c) = heave_err_max_floeNB(j);
 heave_err_min_floeNB_a(c) = heave_err_min_floeNB(j);
 pitch_err_max_floeNB_a(c) = pitch_err_max_floeNB(j);
 pitch_err_min_floeNB_a(c) = pitch_err_min_floeNB(j);
 waveheight_a(c) = waveheight(j);
 wavenum_a(c)    = wavenum(j);
 wavelength_a(c) = wavelength(j);
 c = c+1;
end
 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PLOT EXPERIMENTAL RESULTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SURGE

% hold on
figure(1)
hold on
set(gcf,'position',[100 400 500 400]);
% set(gca,'FontSize',20)
set(gca,'FontSize',28)

f.sg.floeB.d = plot(wavelength/2/L,...
                   surge_floeB./waveheight.*tanh(wavenum*h),...
                   'kx','MarkerSize',10);
f.sg.floeNB.d = plot(wavelength/2/L,surge_floeNB./waveheight.*tanh(wavenum*h),...
                   'ro','MarkerSize',8);

% % PLOT FLOE NB, INDICATE OVERWASH
% for j = 1:length(surge_floeNB)
%  if overwash(j) == 1
%   f.sg.floeNB.d(j) = plot(wavelength(j)/2/L,...
%                         surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),...
%                         'ro','MarkerSize',8);
%  else
% %   f.sg.floeNB.d(j) = plot(wavelength(j)/2/L,...
% %                         surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),...
% %                         'rx','MarkerSize',10);
%  end
% end

% % PLOT ACCORDING TO WAVE HEIGHTS
% % for j=1:length(waveheightR)
% %  if whR(j) == 1
% %   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'ks','MarkerSize',10)
% %  elseif whR(j) == 2
% %   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'ro','MarkerSize',10)
% %  elseif whR(j) == 3
% %   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'b^','MarkerSize',10)
% %  elseif whR(j) == 4
% %   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'+','MarkerSize',10,'Color',[0,0.5,0])
% %  end
% % end
% for j=1:length(waveheightR)
%  if whR(j) == 1
%   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'bs','MarkerSize',10)
%  elseif whR(j) == 2
%   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'m+','MarkerSize',10)
%  elseif whR(j) == 3
%   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'g^','MarkerSize',10)
% %  elseif whR(j) == 4
% %   plot(wavelength(j)/2/L,surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'+','MarkerSize',10,'Color',[0,0.5,0])
%  end
% end

% PLOT ERROR BARS
% FLOE NB
% % no overwash
% f.sg.floeNB.err = errorbar(wavelength_no/2/L,...
%                 surge_floeNB_no./waveheight_no.*tanh(wavenum_no*h),...
%                 surge_err_min_floeNB_no./waveheight_no.*tanh(wavenum_no*h),...
%                 surge_err_max_floeNB_no./waveheight_no.*tanh(wavenum_no*h),'r.');
% all
f.sg.floeNB.err = errorbar(wavelength_a/2/L,...
                surge_floeNB_a./waveheight_a.*tanh(wavenum_a*h),...
                surge_err_min_floeNB_a./waveheight_a.*tanh(wavenum_a*h),...
                surge_err_max_floeNB_a./waveheight_a.*tanh(wavenum_a*h),'ro','MarkerSize',8);
% for j = 1:length(surge_floeNB)
%  if overwash(j) == 1
% %  else
%  f.sg.floeNB.err(j) = errorbar(wavelength(j)/2/L,...
%           surge_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),...
%           surge_err_min_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),...
%           surge_err_max_floeNB(j)./waveheight(j).*tanh(wavenum(j)*h),'r.');
%  end
% end
% FLOE B
f.sg.floeB.err = errorbar(wavelength/2/L,...
                surge_floeB./waveheight.*tanh(wavenum*h),...
                surge_err_min_floeB./waveheight.*tanh(wavenum*h),...
                surge_err_max_floeB./waveheight.*tanh(wavenum*h),'kx','MarkerSize',10);


% set(f.sg.floeNB.err,'LineWidth',1)
% set(f.sg.floeB.err,'LineWidth',1)

ylabel('Surge RAO')
xlabel('Wavelength / Floe length [m/m]')
ylim([0 1.2])
set(gca,'YTick',0:0.2:1.2,'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1' '1.2'})
legend('Barrier','No Barrier','Location','SouthEast')
box on

%% HEAVE
figure(2)
hold on
set(gcf,'position',[600 400 500 400]);
% set(gca,'FontSize',20)
set(gca,'FontSize',28)

f.hv.floeB.d = plot(wavelength/2/L,...
                   heave_floeB./waveheight,'kx','MarkerSize',10);
f.hv.floeNB.d = plot(wavelength/2/L,heave_floeNB./waveheight,'ro','MarkerSize',8);

% % PLOT FLOE NB, INDICATE OVERWASH
% for j = 1:length(surge_floeNB)
%  if overwash(j) == 1
%   f.hv.floeNB.d(j) = plot(wavelength(j)/2/L,heave_floeNB(j)./waveheight(j),'ro','MarkerSize',8);
%  else
%   f.hv.floeNB.d(j) = plot(wavelength(j)/2/L,heave_floeNB(j)./waveheight(j),'rx','MarkerSize',10);
%  end
% end

% % PLOT ACCORDING TO WAVE HEIGHTS
% for j=1:length(waveheightR)
%  if whR(j) == 1
%   plot(wavelength(j)/2/L,heave_floeNB(j)./waveheight(j),'bs','MarkerSize',10)
%  elseif whR(j) == 2
%   plot(wavelength(j)/2/L,heave_floeNB(j)./waveheight(j),'m+','MarkerSize',10)
%  elseif whR(j) == 3
%   plot(wavelength(j)/2/L,heave_floeNB(j)./waveheight(j),'g^','MarkerSize',10)
% %  elseif whR(j) == 4
% %   plot(wavelength(j)/2/L,heave_floeNB(j)./waveheight(j),'+','MarkerSize',10,'Color',[0,0.5,0])
%  end
% end

% PLOT ERROR BARS
% FLOE NB
% % no overwash
% errorbar(wavelength_no/2/L,...
%           heave_floeNB_no./waveheight_no,...
%           heave_err_min_floeNB_no./waveheight_no,...
%           heave_err_max_floeNB_no./waveheight_no,'r.');
% all
errorbar(wavelength_a/2/L,...
          heave_floeNB_a./waveheight_a,...
          heave_err_min_floeNB_a./waveheight_a,...
          heave_err_max_floeNB_a./waveheight_a,'ro','MarkerSize',8);
% for j = 1:length(heave_floeNB)
%  if overwash(j) == 1
%  else
%   f.hv.floeNB.err(j) = errorbar(wavelength(j)/2/L,...
%           heave_floeNB(j)./waveheight(j),...
%           heave_err_min_floeNB(j)./waveheight(j),...
%           heave_err_max_floeNB(j)./waveheight(j),'r.');
%  end
% end
% FLOE B
errorbar(wavelength/2/L,...
          heave_floeB./waveheight,...
          heave_err_min_floeB./waveheight,...
          heave_err_max_floeB./waveheight,'kx','MarkerSize',10);

% set(f.hv.floeNB.err,'LineWidth',1)
% set(f.hv.floeB.err,'LineWidth',1)

ylabel('Heave RAO')
xlabel('Wavelength / Floe length [m/m]')
ylim([0 1.2])
set(gca,'YTick',0:0.2:1.2,'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1' '1.2'})
% set(gca,'YTick',0:0.2:1.2,'YTickLabel',{''})
legend('Barrier','No Barrier','Location','SouthEast')
box on

%% PITCH
figure(3)
hold on
set(gcf,'position',[350 0 500 400]);
% set(gca,'FontSize',20)
set(gca,'FontSize',28)

f.pt.floeB.d = plot(wavelength/2/L,...
                   pitch_floeB./wavenum./waveheight,'kx','MarkerSize',10);
f.pt.floeNB.d = plot(wavelength/2/L,pitch_floeNB./wavenum./waveheight,'ro','MarkerSize',8);

% % PLOT FLOE NB, INDICATE OVERWASH
% for j = 1:length(surge_floeNB)
%  if overwash(j) == 1
%   f.pt.floeNB.d(j) = plot(wavelength(j)/2/L,...
%                         pitch_floeNB(j)./wavenum(j)./waveheight(j),...
%                         'ro','MarkerSize',8);
%  else
%   f.pt.floeNB.d(j) = plot(wavelength(j)/2/L,...
%                         pitch_floeNB(j)./wavenum(j)./waveheight(j),...
%                         'rx','MarkerSize',10);
%  end
% end

% % PLOT ACCORDING TO WAVE HEIGHTS
% for j=1:length(waveheightR)
%  if whR(j) == 1
%   plot(wavelength(j)/2/L,pitch_floeNB(j)./wavenum(j)./waveheight(j),'bs','MarkerSize',10)
%  elseif whR(j) == 2
%   plot(wavelength(j)/2/L,pitch_floeNB(j)./wavenum(j)./waveheight(j),'m+','MarkerSize',10)
%  elseif whR(j) == 3
%   plot(wavelength(j)/2/L,pitch_floeNB(j)./wavenum(j)./waveheight(j),'g^','MarkerSize',10)
% %  elseif whR(j) == 4
% %   plot(wavelength(j)/2/L,pitch_floeNB(j)./wavenum(j)./waveheight(j),'+','MarkerSize',10,'Color',[0,0.5,0])
%  end
% end

% PLOT ERROR BARS
% FLOE NB
% % no overwash
% errorbar(wavelength_no/2/L,...
%           pitch_floeNB_no./wavenum_no./waveheight_no,...
%           pitch_err_min_floeNB_no./wavenum_no./waveheight_no,...
%           pitch_err_max_floeNB_no./wavenum_no./waveheight_no,'r.');
% all
errorbar(wavelength_a/2/L,...
          pitch_floeNB_a./wavenum_a./waveheight_a,...
          pitch_err_min_floeNB_a./wavenum_a./waveheight_a,...
          pitch_err_max_floeNB_a./wavenum_a./waveheight_a,'ro','MarkerSize',8);
% for j = 1:length(pitch_floeNB)
%  if overwash(j) == 1
%  else
%   f.pt.floeNB.err(j) = errorbar(wavelength(j)/2/L,...
%           pitch_floeNB(j)./wavenum(j)./waveheight(j),...
%           pitch_err_min_floeNB(j)./wavenum(j)./waveheight(j),...
%           pitch_err_max_floeNB(j)./wavenum(j)./waveheight(j),'r');
%  end
% end
% FLOE B
errorbar(wavelength/2/L,...
          pitch_floeB./wavenum./waveheight,...
          pitch_err_min_floeB./wavenum./waveheight,...
          pitch_err_max_floeB./wavenum./waveheight,'kx','MarkerSize',10);


% set(f.pt.floeNB.err,'LineWidth',1)
% set(f.pt.floeB.err,'LineWidth',1)
         
ylabel('Pitch RAO')
xlabel('Wavelength / Floe length [m/m]')
ylim([0 1.2])
set(gca,'YTick',0:0.2:1.2,'YTickLabel',{'0' '0.2' '0.4' '0.6' '0.8' '1' '1.2'})
% set(gca,'YTick',0:0.2:1.2,'YTickLabel',{''})
legend('Barrier','No Barrier','Location','SouthEast')
box on




