clear all
import matlab.io.*
%enames= element names, elements meaning fits files
enames = dir('/home/zechariah/Documents/Solar Analogs/Temp_Ev/*.fits');
%numel= number of elements
for ii=1:numel(enames)
    temp = enames(ii).name(21:29);
    fname{ii} = temp;
end
%fname= file name
%unique function makes sure to not return duplicates
fname = unique(fname);

%ones creates an array of ones    
%numel is the number of rows and 1 column 
testc = ones(numel(fname),1);
for ii=1:numel(fname)
    elist(ii) = str2num(cell2mat(fname(ii)));
end
elist = unique(elist);

scnames = dir('/home/zechariah/Documents/Solar Analogs/Temp_Sc/*.fits');
for ii=1:numel(scnames)
    temp = scnames(ii).name(18:26);
    scname{ii} = temp;
end
scname = unique(scname);

testc = ones(numel(scname),1);
for ii=1:numel(scname)
    sclist(ii) = str2num(cell2mat(scname(ii)));
end
sclist = unique(sclist);

ffnames = dir('/home/zechariah/Documents/Solar Analogs/Temp_Sff/*.fits');
for ii=1:numel(ffnames)
    temp = ffnames(ii).name(26:34);
    ffname{ii} = temp;
end
ffname = unique(ffname);    

testc = ones(numel(ffname),1);
for ii=1:numel(ffname)
    fflist(ii) = str2num(cell2mat(ffname(ii)));
end
fflist = unique(fflist);

prnames = dir('/home/zechariah/Documents/Solar Analogs/Temp_Pr/LC/*.fits');
for ii=1:numel(prnames)
    temp = prnames(ii).name(5:13);
    prname{ii} = temp;
end
prname = unique(prname);

testc = ones(numel(prname),1);
for ii=1:numel(prname)
    prlist(ii) = str2num(cell2mat(prname(ii)));
end
prlist = unique(prlist);

fgnames = dir('/home/zechariah/Documents/Solar Analogs/K2Project/outputs/ts/*target1_ts.out');
for ii=1:numel(fgnames)
    temp = fgnames(ii).name(13:21);
    fgname{ii} = temp;
end
fgname = unique(fgname);

testc = ones(numel(fgname),1);
for ii=1:numel(fgname)
    fglist(ii) = str2num(cell2mat(fgname(ii)));
end
fglist = unique(fglist);

z = mintersect(elist,sclist,fflist,prlist,fglist);


%z is now the EPIC list of all overlaps, but now need to cross that list
%with each fname list to figure out which files to bring! [c ia ib] =
%intersect(elist,z); is the approach to use to get the numbers for each.
%enames(ia).name is the list in enames

%ok, start with making file lists for each
[c ia1 ib] = intersect(elist, z); %enames(ia).name are names in the list
for ii = 1:numel(ia1)
    exfiles{ii} = enames(ia1(ii)).name;
end
%exfiles = sort(exfiles);

[c ia2 ib] = intersect(sclist, z); %enames(ia).name are names in the list
for ii = 1:numel(ia2)
    scxfiles{ii} = scnames(ia2(ii)).name;
end
%scxfiles = sort(scxfiles);

[c ia3 ib] = intersect(fflist, z); %enames(ia).name are names in the list
for ii = 1:numel(ia3)
    ffxfiles{ii} = ffnames(ia3(ii)).name;
end


[c ia4 ib] = intersect(prlist, z); %enames(ia).name are names in the list
for ii = 1:numel(ia4)
    prxfiles{ii} = prnames(ia4(ii)).name;
end

[c ia5 ib] = intersect(fglist, z'); %enames(ia).name are names in the list
for ii = 1:numel(ia5)
    fgxfiles{ii} = fgnames(ia5(ii)).name;
end

%OK, success at last! But this only works through C08, because after that
%time stars get observed more than once, so we need to find a way to deal
%with that...

%Next step is to read them all in? or do one at a time? probably the
%latter, and run each step through all the methods

%Step through all 1174 targets...read each in, apply each method, dump out
%results

%Also will need to produce a list of the stars so that we can collect all
%the EPIC catalog information on each one...

%Project Light Curves
%prxfiles has names of files
%path is dir('/home/derek/Documents/Solar Analogs/K2Project/LC/*.fits');

fclose('all');

if z~=0
    for ii=1:numel(prxfiles)
    %for ii=611:900
    %ff_ls_period = ones(1,100);
    %ff_wave_period = ones(1,100);
    %ff_acf_period = ones(1,100);
    %ff_hht_period = ones(1,100);

    %for ii=1:100
    

    data = fitsread(strcat('/home/zechariah/Documents/Solar Analogs/Temp_Pr/LC/',prxfiles{ii}),'binarytable');
    %remove flagged data points
    mdata = data{8};
    flag = data{10};
    time = data{1};
    mdata(flag>1) = [];
    time(flag>1) = [];
    %remove Nan data points
    time(isnan(mdata)) = [];
    mdata(isnan(mdata)) = [];

    %plot(time,mdata,'.')
    %plot(time,mdata/mean(mdata),'.')


    ttime = [min(time):mode(diff(time)):max(time)];
    fflux = interp1(time,mdata,ttime,'pchip');

    k2project_time{ii} = ttime;
    k2project_flux{ii} = fflux;


%%%%%%%%%%%%%%%%%%%%%%%%

    [pr_ls_period(ii),pr_wave_period(ii),pr_acf_period(ii),pr_hht_period(ii),pr_ls_snr(ii),pr_wave_snr(ii),pr_acf_snr(ii),pr_hht_snr(ii)] = calc_periods(ttime,fflux,z(ii),'TEMP_PR');
    pr_mean(ii) = mean(fflux);
    pr_range(ii) = (prctile(detrend(fflux),95)-prctile(detrend(fflux),5))/pr_mean(ii);
    disp(pr_range(ii));
    
%%%%%%%%%%%%%%%%%%%%%%%%

    %OK, next need to set this up so that we loop through for each star
    %also need to return some plot information: wavelet transform, Lomb
    %periodogram, etc.

    %could also include uncertainties on the next go-round
    %need to save 20 data points for each target...4 periods for each pipeline
    %also need to include names of targets and other(?) input information
    %get similar search results from Huber et al. catalog...

    %for ii=1:numel(exfiles)

    %Everest light curves
    data = fitsread(strcat('/home/zechariah/Documents/Solar Analogs/Temp_Ev/',exfiles{ii}),'binarytable');
    %remove flagged data points
    time = data{6};
    mdata = data{2};
    flag = data{5};
    mdata(flag>1) = [];
    time(flag>1) = [];
    %remove Nan data points
    time(isnan(mdata)) = [];
    mdata(isnan(mdata)) = [];

    ttime = [min(time):mode(diff(time)):max(time)];
    fflux = interp1(time,mdata,ttime,'pchip');

    k2everest_time{ii} = ttime;
    k2everest_flux{ii} = fflux;


%%%%%%%%%%%%%%%%%%%%%%%%

    [e_ls_period(ii),e_wave_period(ii),e_acf_period(ii),e_hht_period(ii),e_ls_snr(ii),e_wave_snr(ii),e_acf_snr(ii),e_hht_snr(ii)] = calc_periods(ttime,fflux,z(ii),'TEMP_EV');
    e_mean(ii) = mean(fflux);
    e_range(ii) = (prctile(detrend(fflux),95)-prctile(detrend(fflux),5))/e_mean(ii);
    disp(e_range(ii));
%%%%%%%%%%%%%%%%%%%%%%%%



    %plot(time,mdata,'.')
    %plot(time,mdata/mean(mdata),'.')

    %end
    %dlmwrite('hlsp_everest_k2_llc_201162999-c01_kepler_v2.0_lc.txt',[time mdata],'precision','%.6f','delimiter','\t');

    %for ii=1:numel(scxfiles)
    %K2SC light curves
    %data = fitsread(strcat('/home/derek/Documents/Solar Analogs/K2SC/',scxfiles{ii}),'binarytable');
    %remove flagged data points

    fptr = fits.openFile(strcat('/home/zechariah/Documents/Solar Analogs/Temp_Sc/',scxfiles{ii}),'readonly');
    %fits.getNumHDUs(fptr)
    fits.movAbsHDU(fptr,3); %HDU 3 is the SAP output; HDU 2 is the PDCMAP output
    %imgdata = readImg(fptr)
    %imgdata = fits.readImg(fptr)
    %fits.getColName(fptr,'time')
    %fits.getColName(fptr,'flux')
    time = fits.readCol(fptr,1);
    mdata = fits.readCol(fptr,9);
    flag = fits.readCol(fptr,3);

    fits.closeFile(fptr);

    %time = data{1};
    %mdata = data{6};
    %flag = data{3};

    mdata(flag>1) = [];
    time(flag>1) = [];
    %remove Nan data points
    time(isnan(mdata)) = [];
    mdata(isnan(mdata)) = [];

    %figure(1998)
    %plot(time,mdata)

    ttime = [min(time):mode(diff(time)):max(time)];
    fflux = interp1(time,mdata,ttime,'pchip');

    k2sc_time{ii} = ttime;
    k2sc_flux{ii} = fflux;


%%%%%%%%%%%%%%%%%%%%%%%%

    [sc_ls_period(ii),sc_wave_period(ii),sc_acf_period(ii),sc_hht_period(ii),sc_ls_snr(ii),sc_wave_snr(ii),sc_acf_snr(ii),sc_hht_snr(ii)] = calc_periods(ttime,fflux,z(ii),'TEMP_SC');
    sc_mean(ii) = mean(fflux);
    sc_range(ii) = (prctile(detrend(fflux),95)-prctile(detrend(fflux),5))/sc_mean(ii);
    disp(sc_range());
%%%%%%%%%%%%%%%%%%%%%%%%

    %K2SC light curves
    %data = fitsread(strcat('/home/derek/Documents/Solar Analogs/K2SC/',scxfiles{ii}),'binarytable');
    %remove flagged data points

    fptr = fits.openFile(strcat('/home/zechariah/Documents/Solar Analogs/Temp_Sc/',scxfiles{ii}),'readonly');
    %fits.getNumHDUs(fptr)
    fits.movAbsHDU(fptr,2); %HDU 3 is the SAP output; HDU 2 is the PDCMAP output
    %imgdata = readImg(fptr)
    %imgdata = fits.readImg(fptr)
    %fits.getColName(fptr,'time')
    %fits.getColName(fptr,'flux')
    time = fits.readCol(fptr,1);
    mdata = fits.readCol(fptr,9);
    flag = fits.readCol(fptr,3);

    fits.closeFile(fptr);

    %time = data{1};
    %mdata = data{6};
    %flag = data{3};

    mdata(flag>1) = [];
    time(flag>1) = [];
    %remove Nan data points
    time(isnan(mdata)) = [];
    mdata(isnan(mdata)) = [];

    %figure(1998)
    %plot(time,mdata)

    ttime = [min(time):mode(diff(time)):max(time)];
    fflux = interp1(time,mdata,ttime,'pchip');

    k2scp_time{ii} = ttime;
    k2scp_flux{ii} = fflux;


%%%%%%%%%%%%%%%%%%%%%%%%

    [scp_ls_period(ii),scp_wave_period(ii),scp_acf_period(ii),scp_hht_period(ii),scp_ls_snr(ii),scp_wave_snr(ii),scp_acf_snr(ii),scp_hht_snr(ii)] = calc_periods(ttime,fflux,z(ii),'TEMP_SC-P');
    scp_mean(ii) = mean(fflux);
    scp_range(ii) = (prctile(detrend(fflux),95)-prctile(detrend(fflux),5))/scp_mean(ii);
    disp(scp_range(ii));
%%%%%%%%%%%%%%%%%%%%%%%%

    %plot(time,mdata,'.')
    %plot(time,mdata/mean(mdata),'.')

    %end
    %dlmwrite('hlsp_k2sc_k2_llc_229174815-c102_kepler_v2_lc.txt',[time mdata],'precision','%.6f','delimiter','\t');

    fclose('all');

    %K2SFF light curves
    %for ii=1:numel(ffxfiles)
    data = fitsread(strcat('/home/zechariah/Documents/Solar Analogs/Temp_Sff/',ffxfiles{ii}),'binarytable');
    %remove flagged data points
    time = data{1};
    mdata = data{3};
    flag = data{5};
    mdata(flag>1) = [];
    time(flag>1) = [];
    %remove Nan data points
    time(isnan(mdata)) = [];
    mdata(isnan(mdata)) = [];

    ttime = [min(time):mode(diff(time)):max(time)];
    fflux = interp1(time,mdata,ttime,'pchip');

    k2sff_time{ii} = ttime;
    k2sff_flux{ii} = fflux;



%%%%%%%%%%%%%%%%%%%%%%%%

    [ff_ls_period(ii),ff_wave_period(ii),ff_acf_period(ii),ff_hht_period(ii),ff_ls_snr(ii),ff_wave_snr(ii),ff_acf_snr(ii),ff_hht_snr(ii)] = calc_periods(ttime,fflux,z(ii),'TEMP_SFF');
    ff_mean(ii) = mean(fflux);
    ff_range(ii) = (prctile(detrend(fflux),95)-prctile(detrend(fflux),5))/ff_mean(ii);
    disp(ff_range(ii));
    %[a1,a2,a3,a4] = calc_periods(ttime,fflux);
    %ff_ls_period(ii) = a1;
    %ff_wave_period(ii) = a2;
    %ff_acf_period(ii) = a3;
    %ff_hht_period(ii) = a4;

%%%%%%%%%%%%%%%%%%%%%%%%

    %plot(time,mdata,'.')
    %plot(time,mdata/mean(mdata),'.')
    %end
    %dlmwrite('hlsp_k2sff_k2_lightcurve_201087784-c102_kepler_v1_llc.txt',[time mdata],'precision','%.6f','delimiter','\t');

    %K2SP light curves
    %for ii=1:numel(fgxfiles)
    % strcat('/home/derek/Documents/Solar Analogs/K2Project/outputs/ts/',fgxfiles{ii})
    data=dlmread(strcat('/home/zechariah/Documents/Solar Analogs/Temp_Pr/outputs/ts/',fgxfiles{ii}));
    %remove flagged data points
    time = data(:,1);
    mdata = data(:,3);
    %remove Nan data points
    time(isnan(mdata)) = [];
    mdata(isnan(mdata)) = [];

    ttime = [min(time):mode(diff(time)):max(time)];
    fflux = interp1(time,mdata,ttime,'pchip');

    k2sp_time{ii} = ttime;
    k2sp_flux{ii} = fflux;


%%%%%%%%%%%%%%%%%%%%%%%%

    [fg_ls_period(ii),fg_wave_period(ii),fg_acf_period(ii),fg_hht_period(ii),fg_ls_snr(ii),fg_wave_snr(ii),fg_acf_snr(ii),fg_hht_snr(ii)] = calc_periods(ttime,fflux,z(ii),'TEMP_SP');
    fg_mean(ii) = mean(fflux);
    fg_range(ii) = (prctile(detrend(fflux),95)-prctile(detrend(fflux),5))/fg_mean(ii);

%%%%%%%%%%%%%%%%%%%%%%%%

    %append data to file...

    close all
    save '3Jan_outputB.mat'

    %plot(time,mdata,'.')
    %plot(time,mdata/mean(mdata),'.')
    end
    %dlmwrite('pipeout_ktwo201087784-c101_target1_ts.txt',[time mdata],'precision','%.6f','delimiter','\t');
end
fclose('all');

%also to do: plot results for each method for each light curve (lots of
%plots!! do as both fig and prf and make sure you run with graphics emulation), estimate errors (maybe not for poster...?)

%need to calculate amplitudes and decide if they are "good enough" for each
%fit...some kind of SNR...easy for L-S, ACF?, maybe try whitening by
%differencing, running again, and taking that as WN level?
%or maybe just put out a noise level for each? or get SNR and then
%empirically determine what "good" SNR is for each method based on what we
%see for those measurements that ARE common across all techniques...
