function timeCoursePlotta(ds,gr,doOverwrite,tDrug)
% function timeCoursePlotta(ds,gr,doOverwrite,tDrug)
% is a q&d overview-of-tonic-currents-in-voltage-clamp-data plotter which
% is louzily programmed and needs a thorough upgrade in case it shall be
% used further

doOverlayPlot=false;

mergeFn=abfmerge2mat(ds.fList(:,1),'dDir',ds.dDir,'noMerge',1);
if ~exist([ds.dDir mergeFn],'file') || doOverwrite
  [mergeFn,fi,IN1]=abfmerge2mat(ds.fList(:,1),'dDir',ds.dDir,'cFreq',[nan 10],...
    'sampFreq',500,'gapVal',nan);
else
  load([ds.dDir mergeFn])
end
if 1
  % curb outliers (artifacts)
  outl=prctile(IN1,[.1 99.9]);
  IN1(IN1<outl(1) | IN1>outl(2))=nan;
else
  % better
  etslexcsubst
end

figure(2)
if doOverlayPlot
  hold on
  % subtract median of minutes 8-9
  ix=cont2discrete([8 9]*60*1000,fi.si/1000);
  IN1=IN1-nanmedian(IN1(ix(1):ix(2)));
  % downsample & filter again 
  IN1=IN1(1:10:end);
  fi.si=fi.si*10;
  IN1=medfilt1(double(IN1),50);
  pCol=['kbgrcym'];
  pCol=pCol(randi(7));
else
  clf
  pCol='k';
end


% minutes on x axis
plot((1:numel(IN1))*fi.si/60e6,IN1,pCol);
niceyax
% lines indicating drug applic times
yl=get(gca,'ylim');
lh=line([tDrug;tDrug],[yl;yl]','color','b','linewidth',2);
grid on
xlabel('time (min)')
ylabel('current (pA)');
if ~doOverlayPlot
  title(mergeFn,'interpreter','none')
  if ~isempty(gr.printas),
    print([gr.fDir strrep(mergeFn,'.mat','_tc')],gr.printas);
  end
end

