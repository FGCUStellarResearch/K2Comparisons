function [ls_period,wave_period,acf_period,hht_period,ls_snr,wave_snr,acf_snr,hht_snr] = calc_periods(ttime,fflux,star_no,pipeline)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure(199)
set(figure(199),'visible','off');
plot(ttime,fflux/median(fflux(~isnan(fflux))),'.k','markersize',2)
xlabel('Time (d)','Interpreter','latex')
ylabel('Relative Flux','Interpreter','latex')
title(strcat(num2str(star_no),{' '},pipeline),'Interpreter','latex')
saveas(gcf,strcat('Output/LC/',num2str(star_no),'_',pipeline,'.png'))

%label and save plots with star No & type of pipeline (need to add that to
%pass)
%label and save plot outputs from each fit
% return values for amplitude, noise + P for LS

%1. Lomb-Scargle calculation
[P,f,alpha] = lomb(detrend(fflux(:)),ttime);
min_freq = 2./(ttime(end)-ttime(1));
%[aa, ls_loc] = max(p)
[aa, ls_loc] = findpeaks(P); %find all the peaks
%find all the peaks above min_freq
%find max peak above min_freq
[x1, x2] = find(f(ls_loc)>min_freq); %x1 are locations within b that correspond to peaks
[j1, j2] = max(aa(x1)); %j2 is location in aa/bb array
ls_loc = (ls_loc(x1(j2))); %not working properly :(
fmax = f(ls_loc) %Using fmax later for HHT

%alpha(bb)
if (alpha(ls_loc) < 1e-4)
    ls_period = 1./f(ls_loc);
else
    ls_period = -99;
end
%X0 = dft(ttime,detrend(fflux(:))/mean(fflux),f(bb));
%ls_amp = 2*abs(X0)/numel(ttime);


figure(99)
set(figure(99),'visible','off');
plot(f,P)
text(f(ls_loc),0.95*max(P),strcat('\leftarrow',num2str(ls_period)))
if f(ls_loc) < 5
    xlim([0 max(2,5*f(ls_loc))])
end
title(strcat(num2str(star_no),{' '},pipeline),'Interpreter','latex')
xlabel('Frequency ($\rm c~d^{-1}$)','Interpreter','latex')
ylabel('P','Interpreter','latex')
saveas(gcf,strcat('Output/LS/',num2str(star_no),'_',pipeline,'_','ls','.png'))

%calculate LS "noise"
ls_noise = median(P);
ls_snr = max(P)/ls_noise;

%first, only those with small enough alpha (less than 1e-4?) are chosen
%loop through until only the largest one of those is around? or more simply
%just ensure that it has alpha < 1e-4 or report no period

%then use dft to get peak amplitude! No....just fit a sin curve with known
%frequency...


%%can use peaks and width to find FWHM error...
%%need to fix this so it only searches space of less than time length of
%%segment...and/or checks vs alpha
%%also fix so that NaNs aren't used! go back and remove from the input  TS

%2 Autocorrelation calculation
ta = acf(detrend(fflux(:))/mean(fflux),min(2500,numel(fflux)-1));

[pks,locs] = findpeaks(smoothdata(ta,31),'MinPeakDistance',100);
%maybe limit to smoothing at 31 (or gaussian smooth?) and 
%stop the 100 issue because it loses the long-period peaks that are real!
%pks(locs<100) = [];
%w(locs<100) = [];
%p(locs<100) = [];
%locs(locs<100) = [];
[aa bb] = max(pks); %locs(bb) will be the point of max, where bb is in units of lags
%go to lag length equal to length of ttime

acf_period = locs(bb)*mode(diff(ttime));
if(isempty(bb))
    acf_period = max(ttime)-min(ttime);
end

if (ta(bb)<0)
    acf_period = -99;
end

figure(101)
set(figure(101),'visible','off');
plot(ta)
title(strcat(num2str(star_no),{' '},pipeline),'Interpreter','latex')
xlabel('Delay','Interpreter','latex')
ylabel('Amplitude','Interpreter','latex')

if ta(bb)>0
    text(f(bb),0.95*max(ta),strcat('\leftarrow',num2str(acf_period)))
end
saveas(gcf,strcat('Output/ACF/',num2str(star_no),'_',pipeline,'_','acf','.png'))
%calculate ACF "noise"
acf_peak = abs(ta(locs(bb)));
acf_noise = median(abs(ta));
acf_snr = acf_peak/acf_noise;
if (isempty(acf_snr))
    acf_snr = 0;
end

%aa
%acf_period
%ta2 = acf(diff(detrend(fflux(:)))/mean(fflux),2500);
%figure(200)
%plot(ta2)
%mean(abs(ta2))
%std(ta2)


%need to work on getting amplitude
%max must be positive or there is no period! in that case maybe -99 or NaN?
%also need to handle the case where vector is empty....maybe use full
%length of TS as period?
%w gives uncertainty

%3. Wavelet calculation
[wave,period,scale,coi,dj,paramout,k] = contwt(detrend(fflux)/mean(fflux),mode(diff(ttime)));
awave = abs(wave);
for ii=1:1001
swave(ii) = sum(awave(ii,:));
end

[pks,locs] = findpeaks(swave);
if isempty(pks)
    [pks,locs] = max(swave);
end
%pks(locs<100) = [];
%w(locs<100) = [];
%p(locs<100) = [];
%locs(locs<100) = [];
[aa bb] = max(pks);%locs(bb) will be the point of max, and period(locs(bb)) will be period there
wave_period = period(locs(bb));
%wave_amp = abs(swave(locs(bb)))/(numel(swave)*numel(ttime));
wave_amp = abs(swave(locs(bb)));

%calculate background level
wave_noise = median(swave);
wave_snr = wave_amp/wave_noise;

%%wave amplitude not right yet... FIX

figure(1001)
set(figure(1001),'visible','off');
imagesc(awave/max(awave(:)))
%tempt = [ttime(1):ceil((max(ttime)-min(ttime))/6):max(ttime)];
%pert = [period(1):ceil((max(period)-min(period))/6):max(period)];
colorbar
%set(gca, 'XTick',tempt,'XTickLabel',ttime(tempt))
%set(gca, 'XTick',tempt,'XTickLabel',ttime(tempt)-min(ttime))
%set(gca, 'YTick',[0:ceil(numel(period)/10):numel(period)])
%set(gca, 'XTick',[min(ttime):ttime(end)-min(ttime)],'XTickLabel',[min(ttime):max(ttime)-min(ttime)])
%set(gca, 'YTick',[1:ceil(numel(period)/6):numel(period)],'YTickLabel',round(pert))
set(gca, 'YTick',(1:100:length(period)),'YTickLabel',round(period(1:100:length(period)),2))
%set(gca, 'XTick',[1:numel(tempt)],'XTickLabel',tempt)
%set(gca, 'XTick',[1:ceil(numel(ttime)/6):numel(ttime)], 'XTickLabel',round(tempt))
set(gca, 'XTick',(1:500:length(ttime)), 'XTickLabel',round(ttime(1:500:length(ttime))))
%xtickformat('%.0f')
xlabel('Time (d)')
ylabel('Period (d)')
title(strcat(num2str(star_no),{' '},pipeline),'Interpreter','latex')
colormap jet

saveas(gcf,strcat('Output/wave/',num2str(star_no),'_',pipeline,'_','wav','.png'))

%%maybe try a linear transformation and then go back to the nearest
%%integers? IDK how to do this so it looks ok...

%end
%dlmwrite('ktwo201087784-c102_llc.txt',[time mdata],'precision','%.6f','delimiter','\t');

%4. HHT calculation
%allmode2 = eemd(fflux,0.1,100);
allmode2 = eemd(fflux,0.1,20);
[fm,am] = fa(allmode2,mode(diff(ttime)));

%use results from LS for HHT fm
fm_mean = mean(fm,1);
[afv, fmode] = min(abs(fm_mean - fmax)); %fmode is the mean fm that is closest to fmax
imf_mean = mean(am,1);

%first trial to make hisrogram %%=(lines that were commented out before)
%%find IMF to use by means(am(:,n))
%imf_mean  = mean(am,1);
%[temp, max_imf] = max(imf_mean(5:end-1));
%%[temp max_imf] = findpeaks(max(imf_mean(1:end-1));
%max_imf = max_imf+5;

figure(11)
%set(figure(11),'visible','off');
h = histogram(fm(:,fmode),100,'visible','on','normalization','probability','DisplayStyle','stairs');
%h = histogram(fm(:,max_imf),100,'visible','on','normalization','probability','DisplayStyle','stairs');
%h = histcounts(fm(:,max_imf),100);
[pks, locs] = findpeaks(h.Values);

[temp mpk] = max(pks);
peak_pos = mean(h.BinEdges(locs(mpk):locs(mpk)+1));
hht_period = 1./peak_pos;
hht_amp = (imf_mean(fmode))/mean(fflux);
%hht_amp = imf_mean(max_imf)/mean(fflux);

xlabel('Frequency ($\rm c~d^{-1}$)','Interpreter','latex')
ylabel('Relative Occurrence','Interpreter','latex')
title(strcat(num2str(star_no),{' '},pipeline),'Interpreter','latex')

%text(fm(locs(mpk),max_imf),0.95*peak_pos,strcat('\leftarrow',num2str(hht_period)))
text(peak_pos,0.95*max(pks),strcat('\leftarrow',num2str(hht_period)))
saveas(gcf,strcat('Output/HHT/',num2str(star_no),'_',pipeline,'_','hht','.png'))

%calculate HHT "noise"
hht_noise = std(diff(allmode2(:,fmode)))/mean(fflux);
hht_noise = median(abs(am(:,fmode)))/mean(fflux);
hht_snr = hht_amp/hht_noise;


end

