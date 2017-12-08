function rawplot_PSC_mSeries(dataPath,sumResultDataPath,ap)
% plots of raw data traces and PSC cutouts 


[nDsRow,nDsCol]=size(ds);

nCol=nDsCol+1;


for h=1:nDsCol
  for g=1:nDsRow
    % ---------------------------------------------
    [d,si]=abfload([ds(g,h).dDir '\' ds(g,h).fn],'channels',ds(g,h).ch,'start',ds(g,h).intv(1),'stop',ds(g,h).intv(2));
    % filter
    [d,si]=lofi(d,si,3000,'pickf',sampFac);
    %     d=hifi(d,si,200);
    % determine size
    [n1,n2]=size(d);

    
    % notch filter
    d=elim_hum(d,si,[50 100]);
    
    % get rid of offset
    d=d-median(d);
    
    % reshape & drop points at end
    tmpN1=nSweepPp*floor(n1/nSweepPp);
    if n1-tmpN1>0
      disp([int2str(n1-tmpN1) ' points have to be dropped at end']);
      d=d(1:tmpN1,:);
    end
    n1=tmpN1;
    d=reshape(d(1:n1,:),floor(n1/nSweepPp),nSweepPp);
    % subplot number
    spn=(g-1)*nCol+h;
    subplot(nDsRow,nCol,spn);
    rexy('ax',gca,'xfac',1.1,'yfac',1.15);
    if g==nDsRow && h==nDsCol
      pllplot(d,'si',si,'spacing','fixed','dy',dy,'ylim',yl);
    else
      pllplot(d,'si',si,'spacing','fixed','dy',dy,'ylim',yl,'noscb',1);
    end
  end
end


% print?
if ~isempty(printas)
  print(printas,'-r300',[plotPath mfilename],'-painters');
end
