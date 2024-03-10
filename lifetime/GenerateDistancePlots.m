function GenerateDistancePlots(edgeDist, dist, centroid)


%%How to calculate edgeDist

micronsppixel=1; % need to figure out what size each pixel is; is in xml 
% files but couldnt connect to box

% assumes your mask is 1 and NaN=background if mask is 1/0 then bwdist((~mask))

% figure;imagesc(edgeDist,'AlphaData',mask);

%now you have edgeDist Map and corresponding dist or alpha1/2


OPVals=[];
catRange=[];

totRange=[0	15;15 30;30	45;45 60;60	75; 75	200]; % tune to exact ranges and size that you want
figure

for(dd=1:6)
    rangeIdx=[edgeDist(:)>totRange(dd,1)] & [edgeDist(:)<=(totRange(dd,2))];
    OPVals=[OPVals;dist(rangeIdx)];
    catRange=[catRange;repmat(dd,[length(dist(rangeIdx)),1])];
end
labels={'0-15'	'15-30'	'30-45'	'45-60'	'60-75'	'Core'};

xrange=2*[1:dd]-.7;
b=boxchart(2*catRange-1,OPVals,'MarkerStyle','none');
hold on
% xlim([0 13])
xlabel('Distance from Edge (${\mu}$m)','Interpreter','latex')
set(gca,'FontSize',30,'TickLength',[0 0],'Color','none')
xticks(xrange);
xticklabels(labels);
ylim([0 .5])



x = [];
y = [];

for i = 1:6

    x(i) = 2*i-.7;
    y(i) = mean(OPVals(catRange == i), 'omitnan');


end
plot(x, y)












%% how to calculate distance from the center of image
% %need to first decide how you will determine center of image 
% % Option 1: Centroid of largest continous block, prob useful when
% migrating cells could distort center

% stats=regionprops(~isnan(mask),'all');
% [val,idx]=max([stats.Area]);
% refPt=round(stats(idx).Centroid);

% % Option 2: Center of mass of entire mask
% [r, c] = find(mask == 1);
% refPt = [mean(r), mean(c)];

% % Option 3 user defined 


% 
% figure
% refPt=round(centroid);
% 
% canvas=zeros(512);
% canvas(refPt(1),refPt(2))=1;
% centroidDist=bwdist(canvas);
% 
% OPVals=[];
% catRange=[];
% 
% totRange=[0	15;15 30;30	45;45 60;60	75; 75	200]; % tune to exact ranges and size that you want
% for(dd=1:6)
%     rangeIdx=[centroidDist(:)>totRange(dd,1)] & [centroidDist(:)<=(totRange(dd,2))];
%     OPVals=[OPVals;dist(rangeIdx)];
%     catRange=[catRange;repmat(dd,[length(dist(rangeIdx)),1])];
% end
% labels={'0-15'	'15-30'	'30-45'	'45-60'	'60-75'	'75-200'};
% 
% xrange=2*[1:dd]-.7;
% b=boxchart(2*catRange-1,OPVals,'MarkerStyle','none');
% hold on
% %        xlim([0 13])
% xlabel('Distance from Center (${\mu}$m)','Interpreter','latex')
% set(gca,'FontSize',30,'TickLength',[0 0],'Color','none')
% xticks(xrange);
% xticklabels(labels);


