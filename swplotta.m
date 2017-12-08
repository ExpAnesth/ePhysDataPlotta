function swplotta(ds_in,gr_in)
% ** function swplotta(ds_in,gr_in) 
% a high-level plotting function which creates plots of 'sweeps' of raw
% time series plus optionally derived data with the following properties
% - one or more files
% - ONE channel (except for job 'synConn')
% - derived parameters:
% For information on the input arguments ds_in and gr_in please see
% Part 0, 'define defaults', in this m-file

% *************************************************************************
% This function is a member of a family of functions which perform limited,
% specialized tasks on whole data sets (experiments, each of them composed
% of several files). In contrast to e.g. tslbatch and perievdeal, which are
% called once to work on a whole list of experiments, the functions are
% called for each experiment and place the results (if any) in global
% variables for further use, with varying degrees of safeguards against
% potential problems (duplicate calls on same data, etc.). See the specific
% functions. For orientation, the rough order of evolution is
% [gfplotta,swplotta] -> periburc -> peribuexcit -> pscdeal -> stimdeal
% *************************************************************************

% TO DO
% - maybe compute slope of I-F curve, peak firing rates, etc.
% - there is an inconsistency: ds.fList{:,3} is used for specification of
% the analysis interval in job synconn, whereas ds.chList{:,4} is used for
% the same purpose in job depCurInj

persistent filtC delay fitModel

% ------------ PART 0: define defaults
% --- line plot appearance 
% width of line plot representing data trace
gr.lineW=.5;
% font size (scale bar, title)
gr.fontSz=8;
% scaling factor for post processing, see help of labelscale.m
gr.scaleFac=1;
% marker size
gr.markSz=3;

% --- layout of figure and plot
% layout of figure (usual matlab choices 'tall', 'landscape', etc.)
gr.ornt='tall';

% --- post-plot action settings
% graphics format in which to save plot (the usual input arguments into 
% matlab print function; set to [] for no plot)
% gr.printas='-dpsc2';
gr.printas=[];
% directory to save figures in
gr.fDir=[''];

% --- file specifics
% data directory
ds.dDir=[''];
% addition to results file name
ds.fnPart='';
% unique identifier of cell (last four letters: s=slice, c=cell)
ds.cellID='yyyy_mm_dd_sxcx';
% path & name of file into which computed parameters will be written via
% cellMaster
ds.cellTableFn='';
% file list; column order:
% - file names WITHOUT EXTENSION 
% - sweep selection (list of integers or 'a' for all in file)
% - the interval in POINTS on which analysis (e.g. fit) shall work
% - start and step values of current injection (job 'depCurInj')
ds.fList={...
  'file name', 'a', [nan nan], [nan nan];...
  };
% channel; column order:
% - name
% - type of analysis to perform
% - threshold for spike detection based on slope, mV/ms
% - spx analysis interval in ms 
% - y limits in plot
ds.chList={'channel name','analysis type', nan, nan, [nan nan]};
% name of image file(s) (up to two) depicting neuron
ds.imFile={''};
% downsampling factor
ds.sampFac=1;
% offset that should be added to data
ds.offset=0;
% factor by which data should be multiplied
ds.mFac=1;
% the minimal time (ms) the signal SLOPE has to be above threshold to be
% considered a spike
ds.minSpxSlopeWid=nan;

% ------------ PART 1: check input, adjust parameters, preparations
ds=checkFields(ds,ds_in);
gr=checkFields(gr,gr_in);  
% **
job=ds.chList{1,2};

% ** define a few fixed parameters:
% order of differentiator filter
diffFiltOrd=20;
% - cutout interval for spx (ms) relative to ts 
intvCutout=[-1.5 5];
% - interval for spx slope cutouts (ms) to be used for validation (must
% include complete rising phase of spx)
shortIntvCutout=[0 1.0];
% 'pre' cutout interval for spx in pts (used for determining exact
% beginning of spike and threshold Em)
preIntvCutout=[-.5 0];
% slope threshold value defining beginning of spike (mV/ms), regardless of detection
% threshold
spxTZeroSlopeThresh=20;
% - interval for calculation of base line (important for half width)
intvBaseLine=[-.2 0.0];
% - minimal number of spx in a train for an exponential fit of 1/ISI decay
minNSpxISIFit=8;
% - minimal fit quality of fit explained above expressed in adjusted R square
minIfrFitQuality=.6;
% - fraction of spx at end of train from which adapted/basal firing rate is
% to be computed
baseFrFract=.3;
% - limits of spx magnitude (mV), slope (mV/ms) and half width (ms) which
% shall be represented in plots
limitSpxAmp=[-75 35];
limitSpxSlope=[-200 600];
limitSpxHalfWid=[.2 3];



% ** parameters to be extracted and written to cellTable
Em=[];
spikeThresh=[];
spikeHalfWid=[];
spikeAmp=[];
rheobase=[];
IRData=[];

% --- graphics
labelscale('fontSz',gr.fontSz,'scaleFac',gr.scaleFac,'lineW',gr.lineW,'markSz',gr.markSz); 
nFile=size(ds.fList,1);
cm=coma('blueblackorange','n',16);
cm=cat(3,cm(1:2:end,:),flipud(cm(2:2:end,:)));
cm=permute(cm,[2 3 1]);
cm=reshape(cm,[3 16]);
cm=cm';

% ------------ PART 2: load, process and plot data
for g=1:nFile
  % determine number of sweeps and load as many as possible
  sweepList=ds.fList{g,2};
  [d,si,fi]=abfload([ds.dDir ds.fList{g,1} '.abf'],'info');
  nAllSweep=numel(fi.sweepStartInPts);
  disp(['**** number of available sweeps: ' int2str(nAllSweep)]);
  if isnumeric(sweepList)
    % last N sweeps requested
    if any(sweepList<0)
      if numel(sweepList)>1
        error('sweep list: only single negative values allowed ');
      else
        sweepList=nAllSweep+sweepList+1:nAllSweep;
      end
    end          
    if ~isempty(setdiff(sweepList,1:nAllSweep))
      warning('at least some of the requested sweeps do not exist - adjusting list')
      sweepList=intersect(sweepList,1:nAllSweep);
      % if no overlap between requested and existent sweeps sth is really
      % foul, so give up
      if isempty(sweepList)
        error('no overlap between requested and existing sweeps');
      end
    end
  else
    sweepList=1:nAllSweep;
  end
  [d,si]=abfload([ds.dDir ds.fList{g,1} '.abf'],'channels',ds.chList(:,1),'sweeps',sweepList);
  d=d(1:ds.sampFac:end,:,:)*ds.mFac+ds.offset;
  si=si*ds.sampFac;
  [n1 nChan nSweep]=size(d);
  % don't forget to adjust analysis interval to downsampling factor
  ds.fList{g,3}=floor(ds.fList{g,3}/ds.sampFac);
  % next thing to do: sort channels according to order in ds.chList (we
  % have to do it in a loop because intersect sorts the results) 
  tmpIx=nan(nChan,1);
  requestChNames=fi.recChNames(ismember(fi.recChNames,ds.chList(:,1)));
  for chIx=1:nChan
    tmpIx(chIx)=find(ismember(requestChNames,ds.chList(chIx,1)));
  end
  d=d(:,tmpIx,:);
  
  % initialize figure(s)
  fh=mkfig(g,'v');
  clf, orient(gr.ornt);
  % figure name: file name & first (two) channel name(s)
  tmp=[ds.chList{1:min(2,size(ds.chList,1)),1}];
  figName=[ds.fList{g,1} '_' tmp(~isspace(tmp))];

  switch job
    case 'synConn'
      pcol={'b','r'};
      % the factor by which presynaptic spx will be multiplied so that they
      % are approximately to scale with PSPs
      sFac=.02;
      % y offset of PSP traces (mV; this is also the limit of the PSP
      % traces)
      dy=-4;

%       sFac=.06;
%       dy=-12;
     
      intv=ds.fList{g,3};
      % exchange columns of 1st half of sweep so that presynaptic traces
      % are in one column
      d(1:intv(1,2)+1,:,:)=flipdim(d(1:intv(1,2)+1,:,:),2);
      % also, as pllplot plots traces downwards, flip along 3rd dim
      d=flipdim(d,3);
      for chIx=1:2
        subplot(nFile,2,(g-1)*2+chIx), hold on
        % baseline-subtracted postsynaptic trace
        tmpD=permute(d(intv(chIx,1):intv(chIx,2),2,:),[1 3 2]);
        baseL=prctile(tmpD,10);
        tmpD=tmpD-repmat(baseL,size(tmpD,1),1);
        % set anything exceeding y offset (spx) to nan
        tmpD(tmpD>abs(dy(1)))=nan;
        % plot (dy fixed)
        [ylim,dy,yscaleFac,ph]=pllplot(tmpD,...
          'si',si,'spacing','fixed','dy',dy,'noscb',1);
        set(ph,'color',pcol{mod(chIx,2)+1});
        th=text(0*(1:numel(baseL)),cumsum([0 dy])+dy(1)/8,strvcat(int2str(round(baseL'))));
        % scaled, baseline-subtracted presynaptic trace
        tmpD=permute(d(intv(chIx,1):intv(chIx,2),1,:),[1 3 2])*sFac;
        baseL2=prctile(tmpD,10);
        tmpD=tmpD-repmat(baseL2,size(tmpD,1),1);
        % overplot 
        [ylim,dy,yscaleFac,ph]=pllplot(tmpD,...
          'si',si,'spacing','fixed','dy',dy,'ylim',ylim);
        set(ph,'color',pcol{chIx});
        title([ds.fList{g,1} ', ' ds.chList{mod(chIx,2)+1,1} '->' ds.chList{chIx,1}],'interpreter','none');
      end
      if g==nFile && ~isempty(gr.printas),
        print(gr.printas,'-r300','-f1',[gr.fDir figName '_' job '_' ds.fnPart]);
      end
      
    case 'depCurInj'
      if nChan>1
        error('job ''depCurInj'' works with only one channel')
      end
      % spike threshold
      spxThresh=ds.chList{:,3};
      if spxThresh<spxTZeroSlopeThresh
        warning('spike detection threshold is lower than spike tzero threshold - adjusting to spike tzero threshold')
        spxThresh=spxTZeroSlopeThresh;
      end
      % spike analysis interval
      anIntv=ds.chList{:,4};
      anIntvPts=cont2discrete(anIntv,si/1000,'intv',1);
      % interval from which to compute Em: an interval of 50 ms 10 ms
      % before stim, if possible
      EmIntv=max([0 -60+anIntv(1)]);
      EmIntv(2)=min([EmIntv+50 anIntv(1)]);
      EmIntvPts=cont2discrete(EmIntv,si/1000,'intv',1);
      % the minimal time the signal SLOPE has to be above threshold to be
      % considered a spike in pts
      minSpxSlopeWidPts=cont2discrete(ds.minSpxSlopeWid,si/1000,'intv',1);
      % list of current step amplitudes
      tmpC=ds.fList{g,4};  
      curStepAmp=tmpC(1):tmpC(2):tmpC(1)+(nAllSweep-1)*tmpC(2);
      curStepAmp=curStepAmp(sweepList);
      % generate x and y coordinates to be used for lines depicting current
      % injection steps
      curStepCoXPts=[1 anIntvPts([1 1 2 2]) n1]';
      curStepCoYPts=[zeros(2,nSweep); curStepAmp; curStepAmp; zeros(2,nSweep)];
      % model for fitting inverse of ISI
      if isempty(fitModel)
        fitModel=fittype('a1*exp(-t/tau1)+o',...
          'independent',{'t'},...
          'coefficients',{'a1','tau1','o'});
      end
      
      % --- plot sweeps, one stacked on top of the other
      d=permute(d,[1 3 2]);
      sph_wf=subplot('position',[0.10 0.25 0.24 0.7]);
      % plot hyperpolarizing responses at bottom
      if mean(diff(curStepAmp))>0
        [ylim,dy,~,~,xUnitFac]=pllplot(fliplr(d),'si',si);
        % y offset for text (firing rates) to be placed above waveforms
        dy=fliplr(cumsum([0 dy]))+prctile(d,95);
      else
        [ylim,dy,~,~,xUnitFac]=pllplot(d,'si',si);
        % y offset for text (firing rates) to be placed above waveforms
        dy=cumsum([0 dy]+prctile(d,95));
      end
      title([ds.fList{g,1} ', ' ds.chList{1,1}],'interpreter','none');
      % plot current steps
      sph=subplot('position',[0.10 0.05 0.24 0.15]);
      ph=plot(curStepCoXPts,curStepCoYPts,'k');
      set(ph,'linewidth',.8);
      niceyax;
      axis off
      
      % --- compute Em:
      Em=cat(1,Em,median(d(EmIntvPts,:))');
      
      % --- extract spx waveforms:
      % cutout interval for spx in pts
      intvCutoutPts=cont2discrete(intvCutout,si/1000,'intv',1);
      % 'short' cutout interval for spx in pts (used for differentiating
      % between artifacts and spikes)
      shortIntvCutoutPts=cont2discrete(shortIntvCutout,si/1000,'intv',1);
      % 'pre' cutout interval for spx in pts (used for determining exact
      % beginning of spike and threshold Em)
      preIntvCutoutPts=cont2discrete(preIntvCutout,si/1000,'intv',1);
      helperArr=(preIntvCutoutPts(2):-1:preIntvCutoutPts(1))';
      % cutout interval for base line
      intvBaseLinePts=cont2discrete(intvBaseLine,si/1000,'intv',1);
      % detect spx based on slope 
      if isempty(filtC)
        % implement differentiator filter for less noise compared to diff
        filtC=designfilt('differentiatorfir','FilterOrder',diffFiltOrd, ...
          'PassbandFrequency',6000,'StopbandFrequency',8000, ...
          'SampleRate',1e6/si);
        delay=mean(grpdelay(filtC));
      end
      diffD=filter(filtC,d)/(si/1000);
      % delete samples to get rid of delay
      diffD(1:delay,:)=[];
      % substitute samples to get rid of transient
      diffD(1:diffFiltOrd,:)=0;
      
      % --- detect spx, refine time stamp
      tsl=tcd(diffD,'idx',repmat(spxThresh,1,nSweep));
      % first thing to do: define start of spike as the first point above
      % a fixed slope threshold, looking backwards from current point of
      % detection
      for sweepI=1:nSweep
        if ~isempty(tsl{sweepI})
          % cut out waveform around point of detection
          [tmpCut,isCutout,tsl(sweepI)]=tsl2exc(diffD(:,sweepI),'idx',tsl(sweepI),'win',preIntvCutoutPts);
          % going back from point of detection, identify last point below
          % fixed slope thresh and kick any spike not fulfilling this
          % criterion (note the unique(c))
          [r,c]=find(flipud(tmpCut)<=spxTZeroSlopeThresh);
          % adjust tsl
          tsl{sweepI}=tsl{sweepI}(unique(c))+helperArr(r(logical([1; diff(c)])));
        end
      end
      % remove spx outside analysis interval
      tsl=cellfun(@(x) x(x>=anIntvPts(1) & x<=anIntvPts(2)),tsl,'UniformOutput',false);

      % --- compute firing properties:
      %       % parameters listed in columns of IRData
      %       IRDataPar={...
      %         'curStepAmp',
      %         'avFr',
      %         'peakFr',
      %         'baseFr',
      %         'tauFiringFreq'
      %         };
      % array for all spx
      allSpx=[];
      allSpxSlope=[];
      allBase=[];
      % ...and corresponding sweep indexes
      allSweepIx=[];
      curIRData=[curStepAmp' nan(nSweep,4)];
      % container for fits
      collFitD=cell(nSweep,1);
      % loop over sweeps
      for sweepI=1:nSweep
        % first, remove spurious spx: 
        if ~isempty(tsl{sweepI})
          % any spike slope which remains above threshold for less than
          % ds.minSpxSlopeWid ms is likely an artifact or a low-amplitude
          % spike
          [tmpSpxSlope,~,tsl(sweepI)]=tsl2exc(diffD(:,sweepI),'idx',tsl(sweepI),'win',shortIntvCutoutPts);
          % § note that this is the total time, not the width
          badIx=sum(tmpSpxSlope>spxThresh)<minSpxSlopeWidPts;
          if any(badIx)
            tsl{sweepI}=tsl{sweepI}(~badIx);
          end
        end
        if isempty(tsl{sweepI})
          fr=0;
          curIRData(sweepI,2)=fr;
        else
          tsl{sweepI}=tsldeadt(tsl{sweepI},cont2discrete(1,si/1000),'include_1st',1);
          nSpx=numel(tsl{sweepI});
          % average firing rate in Hz 
          fr=numel(tsl{sweepI})*1000/diff(anIntv);
          curIRData(sweepI,2)=fr;
          if nSpx>=2
            % tsl in ms
            tmpTsl=tsl{sweepI}*(si/1e3);
            % starting at time zero (may be needed for fitting below)
            tsOffs=tmpTsl(1);
            tmpTsl=tmpTsl-tsOffs;
            % 'instantaneous firing rate' (1/ISI) in Hz (!)
            ifr=1000./diff(tmpTsl);
            % peak firing rate in Hz:
            curIRData(sweepI,3)=max(ifr);
            % adapted/basal firing rate: compute as average of last baseFrFract
            % of ISI
            tmpN=max(1,round((nSpx-1)*baseFrFract));
            curIRData(sweepI,4)=mean(ifr(end-tmpN+1:end));
          end
          % fit ifr (inverse of ISI) to exponential + offset
          if nSpx>=minNSpxISIFit
            fitX=tmpTsl(1:end-1);
            fitY=ifr;
            fitOpt=fitoptions('method','NonLinearLeastSquares',...
              'Startpoint',[max(ifr)-min(ifr) 50 min(ifr)],...
              'lower',[0 0 0],...
              'upper',[1000 10000 1000]);
            [respFit,gof]=fit(fitX,fitY,fitModel,fitOpt);
            % if fit is sufficiently good and tau not longer than the
            % duration of the curren step, extract tau
            tmp=coeffvalues(respFit);
            if gof.adjrsquare>=minIfrFitQuality && tmp(2)<= diff(anIntv)
              curIRData(sweepI,5)=tmp(2);
              % save eval'd fit for plot - don't forget to add time offset
              fitEvalX=linspace(fitX(1),fitX(end),100);
              collFitD{sweepI}=[fitEvalX'+tsOffs feval(respFit,fitEvalX)];
            end
          end

          % --- cutouts of spx, spx slope  and base line (making sure that
          % all match), the latter already averaged across time
          [spx,isCutout,tsl(sweepI)]=tsl2exc(d(:,sweepI),'idx',tsl(sweepI),'win',intvCutoutPts);
          [spxSlope]=tsl2exc(diffD(:,sweepI),'idx',tsl(sweepI),'win',intvCutoutPts);
          base=mean(tsl2exc(d(:,sweepI),'idx',tsl(sweepI),'win',intvBaseLinePts));
          % append
          allSpx=cat(2,allSpx,spx);
          allSpxSlope=cat(2,allSpxSlope,spxSlope);
          allBase=cat(2,allBase,base);
          % sweep indexes
          allSweepIx=cat(2,allSweepIx,repmat(sweepI,1,nSpx));
          % positive firing rates as text
        end
        if fr>0
          subplot(sph_wf)
          text(0,dy(sweepI),[int2str(round(fr)) ' Hz']);
        end
      end
      % concatenate IRData
      IRData=cat(1,IRData,curIRData);
      
      % --- compute properties of spx: 
      % allSpxSlope=diff(allSpx)/(si/1000);
      % subtract base line: 
      tmpAllSpx=allSpx-repmat(allBase,diff(intvCutoutPts)+1,1);
      % - min, max, amplitude, half-width
      r=evdeal(tmpAllSpx,'idx',{'minmax','fwidth'},'pTy','maxPeak','rH',.5);
      % - in each file, pick spx with the shortest half width (below the 5th
      % percentile) among those with an amplitude above the median
      threshWidth=prctile(r.fwidth,5);
      threshAmp=prctile(r.totampl,50);
      goodSpxIx=find(r.fwidth<=threshWidth & r.totampl>=threshAmp);
      if isempty(goodSpxIx)
        % in case there is no overlap of sets of spikes as defined above
        % simply pick the thinnest spike
        [~,goodSpxIx]=min(r.fwidth);
      end
      % spike thresh is identical to base line of 'good' spx
      spikeThresh=cat(1,spikeThresh,allBase(goodSpxIx)');
      % spike width: scale to ms and append to spikeHalfWid
      spikeHalfWid=cat(1,spikeHalfWid,r.fwidth(goodSpxIx)'*(si/1000));
      % spike amplitude
      spikeAmp=cat(1,spikeAmp,r.maxPeak(goodSpxIx)');
      
      % --- plots
      % time axis in ms
      t_ms=discrete2cont(intvCutoutPts(1):intvCutoutPts(end),si/1000);
      % spx waveform
      subplot(4,3,2), hold on
      plot(t_ms,allSpx,'k');
      % overplot selected ones
      plot(t_ms,allSpx(:,goodSpxIx),'r');
      axis tight
      set(gca,'xgrid','on','xtick',-5:10,'ylim',limitSpxAmp);
      xlabel('time (ms)');
      ylabel('E_m (mV)');
      
      % spx amplitude vs slope 
      subplot(4,3,3), hold on
      plot(allSpxSlope,allSpx,'k');
      % overplot selected ones
      plot(allSpxSlope(:,goodSpxIx),allSpx(:,goodSpxIx),'r');
      axis tight
      set(gca,'xgrid','on','ylim',limitSpxAmp,'xlim',limitSpxSlope);
      xlabel('dE_m/dt (mV/ms)');      
      ylabel('E_m (mV)');

      % spx slope
      subplot(4,3,5), hold on
      plot(t_ms,allSpxSlope,'k');
      % overplot selected ones
      plot(t_ms,allSpxSlope(:,goodSpxIx),'r');
      axis tight
      set(gca,'xgrid','on','xtick',-5:10,'ylim',limitSpxSlope);
      xlabel('time (ms)');
      ylabel('dE_m/dt (mV/ms)');
      
      % spx amplitude vs half-width
      subplot(4,3,6), hold on
      plot(r.fwidth*(si/1000),r.totampl,'ko')
      % overplot selected ones
      plot(r.fwidth(goodSpxIx)*(si/1000),r.totampl(goodSpxIx),'ro')
      axis tight
      set(gca,'xgrid','on','xtick',0:.5:10,'xlim',limitSpxHalfWid,'ylim',[30 diff(limitSpxAmp)]);
      xlabel('width @half height (ms)');
      ylabel('magnitude (mV)');

      % 1/isi vs ts
      subplot(4,3,8), hold on
      for sweepI=1:nSweep
        if numel(tsl{sweepI})>1
          tslMs=tsl{sweepI}*(si/1e3);
          ph=plot(tslMs(1:end-1),1000./diff(tslMs),'o');
          tmpD=collFitD{sweepI};
          if ~isempty(tmpD)
            col=get(ph,'color');
            ph=plot(tmpD(:,1),tmpD(:,2),'-');
            set(ph,'color',col,'linewidth',0.5);
          end
        end
      end
      nicexy0ax;
      grid on
      xlabel('spike time (ms)');
      ylabel('1/ISI (s^{-1})');

      % I-F curve
      subplot(4,3,9), hold on
      % lines & symbols for positive fr values
      ix=curIRData(:,2)>0;
      plot(curIRData(ix,1),curIRData(ix,2),'ko-');      
      % only symbols for zero fr values
      plot(curIRData(~ix,1),curIRData(~ix,2),'ko');
      % ignore negative current injections for this plot
      axis([0 max(curIRData(:,1))*1.05 0 max(curIRData(:,2))*1.05])
      grid on
      xlabel('current (pA)');
      ylabel('avg fr (Hz)');

      % up to two image files
      for iIx=1:numel(ds.imFile)
        subplot(4,3,10+iIx);
        if isempty(ds.imFile{iIx})
          smarttext('no image file',.08,.5);
          axis off
        else
          try
            img=imread(ds.imFile{iIx});
            image(img);
            colormap(gray(256))
            axis image off
          catch
            smarttext('image file not found',.08,.5);
            axis off
          end
        end
      end
      if ~isempty(gr.printas),
        % print(gr.printas,'-r300',[gr.fDir figName '_' job '_' ds.fnPart]);
        print(gr.printas,'-r300',[gr.fDir figName '_' job]);
      end
  end
  drawnow
end

% after-work jobs:
if strcmp(job,'depCurInj')
  
  % 1. consolidate IRData by combining data sets abd averaging across
  % repetitions of specific current values
  nCol=size(IRData,2);
  newIRData=unique(IRData(:,1));
  newIRData(:,2:nCol)=nan;
  for g=1:size(newIRData,1)
    ix=IRData(:,1)==newIRData(g,1);
    newIRData(g,2:nCol)=nanmean(IRData(ix,2:nCol),1);
  end
  % to make 100% sure, sort explicitly
  IRData=sortrows(newIRData,1);
  clear newIRData;
  
  % 2. compute maximal firing rate, rheobase
  maxFiringRate=max(IRData(:,2));
  % determine minimal stimulation strength leading to spiking (upper bound
  % of rheobase) 
  ix=find(IRData(:,2)>0,1);
  % only if firing rate at that stim strength is below 20% of the maximal
  % firing rate accept the value
  if IRData(ix,2)/maxFiringRate<=.2
    rheobase=IRData(ix,1);
  else
    rheobase=nan;
  end
  
  % ** write parameters
  cellmaster_writepar(ds.cellTableFn,ds.cellID,...
    {'restEm','spikeThresh','spikeAmp','spikeHalfWid','rheobase','maxFiringRate','IRData'},...
    {mean(Em),mean(spikeThresh),mean(spikeAmp),mean(spikeHalfWid),rheobase,maxFiringRate,IRData});
end


% ------------ LOCAL FUNCTIONS ------------------------------------------
% ------------ LOCAL FUNCTIONS ------------------------------------------
function f=checkFields(f,f_in)
s=fieldnames(f);
s_in=fieldnames(f_in);
sect=setdiff(s_in,s);
if ~isempty(sect)
  errordlg({'nonmatching field names(s): '; char(sect)});
end
% assign input values of input struct to current one
for g=1:numel(s_in)
  f.(s_in{g})=f_in.(s_in{g});
end

