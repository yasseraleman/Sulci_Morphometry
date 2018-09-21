load('/media/COSAS/scripts/Sulcal_Processing/Final_Sulc_Processing/matlab_toplay.mat');


correctNodes = ['F.I.P._left'];
inds = find(ismember(cellstr(SulcLabels),cellstr(correctNodes)));


% h = Plot_Surf(PialSurfMat,'FigID',gcf)
% h = Plot_Surf(HullSurfMat,'FigID',gcf);h(1).FaceAlpha = 1;h(1).FaceColor = [1 1 1];



% Surfj = Compound_Surf(Surfsulci(TmtktriIds(inds)));



SulcSurfMat = Surfsulci(TmtktriIds(inds));

[Sulcmet, Surfo, SurfL] = Compute_Node_Metrics_toplay(SulcSurfMat, PialSurfMat, HullSurfMat);



Surfj = Compound_Surf(Surfsulci(TmtktriIds(inds)));

Graph = surface2graph(Surfj,1);
LabNet = Label_Graph_Components(Graph);





