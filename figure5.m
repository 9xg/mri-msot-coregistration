tumorFiles = {'K81_1L1R_bigtum_OE_DCE.mat','K81_1L1R_smalltum_OE_DCE.mat','K81_1R_OE_DCE.mat'}

for k=1:length(tumorFiles)
   load(['dceOe-comparison/Pixel-distributions/' tumorFiles{k}]) 
   
   x = [AuC_inNonResp;AuC_inResp]
   g1 = repmat({'Non-responding'},length(AuC_inNonResp),1);
   g2 = repmat({'Responding'},length(AuC_inResp),1);
   g = [g1; g2];
   
   %figure;boxplot(x,g,'Colors',[0 0 0;0 0 0])
   figure;distributionPlot({AuC_inNonResp;AuC_inResp},'xNames',{'Non-responding','Responding'},'xValues',[-1,1],'distWidth',0.8)
    set(findobj(gca,'type','line'),'linew',3)
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    set(lines, 'Color', 'r');
    plt = Plot();
    plt.ShowBox = false;
    plt.BoxDim = [4 4];
    plt.YLabel = 'DCE value';
    plt.XLim = [-2 2];
    plt.YLim = [-5000 100000];
    plt.XMinorTick = false;
    plt.TickLength = [0.01 0.01]
    plt.export([tumorFiles{k} '_dceInOeMaps.pdf']);
 
end