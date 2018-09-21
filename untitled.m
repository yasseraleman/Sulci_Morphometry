



Surf = Read_Surface('/media/COSAS/Test/freesurfer/fsaverage/surf/lh.pial','curv');

S = ['/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/subjects/1Test_HCP_899885-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/segmentation/mesh/1Test_HCP_899885-20140807-T1wMPR1_LFS_hemi.mesh'];
Surf = load_mesh(S);
S = ['/media/MyDisk/PROCESSING_RESULTS/8-BrainVisaDataBase/subjects/1Test_HCP_899885-20140807-T1wMPR1/t1mri/default_acquisition/default_analysis/segmentation/mesh/1Test_HCP_899885-20140807-T1wMPR1_LFS_hemi_sulci.mesh'];
Surfsulci = load_mesh(S);
[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(Surf.SurfData.vertices,Surf.SurfData.faces);
Surf.Is = Cmean;

% ind  = find(Surf.Is > 0);
% temp = sum(ismember(Surf.SurfData.faces,ind),2);
% ind2rem = find(sum(ismember(Surf.SurfData.faces,ind),2) >0);
% 
% 
% 
% indexbound = Surf.SurfData.faces(ind2rem,:);
% indexbound = unique(indexbound(:));
% indexbound(find(ismember(indexbound,ind))) = [];
% 
% 
% Surft = Surf
% Surft.SurfData.faces(ind2rem,:) = [];
% Surft = Reorg_Surf(Surft);

opts.curvth = 0;

indp = find(Surf.Is > opts.curvth); % Positive Curvature: Sulci
indn = find(Surf.Is <= opts.curvth); % Negative Curvature: Gyri

% Curvature Map
curvmap = Surf.Is;

Surf.Is = zeros(length(curvmap),1); % Gyral regions are labeled with 0
Surf.Is(indn,1) = 1; % Sulci regions are labeled with 1

%% =================== Computing Neighoborhoods ======================== %%
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;


%% =================== Finding Local Curvature Minima ================= %%
% distmap = -1*curvmap;
% temp1 = distmap(temp);
% temp1(indz) =0;
% t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
% indminima = find(sum(ismember(t,-1),2) - sum(logical(t),2) == 0);
% hold on;
% plot3(Surf.SurfData.vertices(indminima,1),Surf.SurfData.vertices(indminima,2),Surf.SurfData.vertices(indminima,3),'.k','MarkerSize',20);


%% =================== Readjusting Curvature threshold ================= %%
temp1 = Surf.Is(temp);
temp1(indz) =  max(temp1(:))+1;
NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
NewMat(indz) = 0;
boundPoints = logical(sum(logical(NewMat)')');
indtemp = find((Surf.Is == 0)&(boundPoints == 1));
readCurvth = mean(curvmap(indtemp));

indp = find(curvmap > readCurvth); % Positive Curvature: Sulci
indn = find(curvmap <= readCurvth); % Negative Curvature: Gyri

Surf.Is = zeros(length(curvmap),1); % Gyral regions are labeled with 0
Surf.Is(indn,1) = 1; % Gyral regions are labeled with 1


%% ================= Labelling Sulcal Clusters ======================== %%
% % % % % % % % set(0,'RecursionLimit',1000);
% % % % % % % % [labid] = Recur_Corr(Surf,0,zeros(size(Surf.Is)),1);
% % % % % % % % Surf.Is = labid;
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % sts = unique(nonzeros(Surf.Is));
% % % % % % % % 
% % % % % % % % for i  = 1:length(sts);
% % % % % % % %     
% % % % % % % %     ind = find(Surf.Is == sts(i));
% % % % % % % %     %     temp = sum(ismember(Surf.SurfData.faces,ind),2);
% % % % % % % %     %     indexbound = find(temp <3&temp>0);
% % % % % % % %     %     indexbound = Surf.SurfData.faces(indexbound,:);
% % % % % % % %     %     indexbound = unique(indexbound(:));
% % % % % % % %     %     indexbound(find(ismember(indexbound,ind))) = [];
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     Surft = Surf;
% % % % % % % %     Surft.Is = curvmap;
% % % % % % % %     
% % % % % % % %     % Selecting only the faces in the sulci
% % % % % % % %     %indneg = find(Surft.Is <0);
% % % % % % % %     a = ismember(Surft.SurfData.faces,ind);
% % % % % % % %     faces2rem = find(sum(a,2) ~= 3);
% % % % % % % %     Surft.SurfData.faces(faces2rem,:) = [];
% % % % % % % %     a = ismember(ind,Surft.SurfData.faces(:));
% % % % % % % %     ind(find(a == 0)) = [];
% % % % % % % %     
% % % % % % % %     % Reorganicing surface
% % % % % % % %     [Surft] = Reorg_Surf(Surft);
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     distmap = -1*Surft.Is;
% % % % % % % %     [Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
% % % % % % % %     Temp = sum(Trip);
% % % % % % % %     Trip(:,Temp==0) = [];
% % % % % % % %     temp = Trip(:,3:end);
% % % % % % % %     indz = find(temp == 0);
% % % % % % % %     temp(indz) = 1;
% % % % % % % %     
% % % % % % % %     temp1 = distmap(temp);
% % % % % % % %     temp1(indz) =0;
% % % % % % % %     t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
% % % % % % % %     indminima = find(sum(ismember(t,-1),2) - sum(logical(t),2) == 0);
% % % % % % % % %     hold on;
% % % % % % % % %     plot3(Surft.SurfData.vertices(indminima,1),Surft.SurfData.vertices(indminima,2),Surft.SurfData.vertices(indminima,3),'.y','MarkerSize',20);
% % % % % % % % Surfo = Surft;
% % % % % % % % cont = 0;
% % % % % % % %     while (cont < 300)% ~isempty(indminima)||
% % % % % % % %         cont = cont + 1
% % % % % % % %         temp(indz) = 1;
% % % % % % % %         X = Surft.SurfData.vertices(:,1);
% % % % % % % %         X = X(temp);
% % % % % % % %         X(indz) =0;
% % % % % % % %         X = X(indminima,:);
% % % % % % % %         tempX = sum(X,2);
% % % % % % % %         Xo = tempX./sum(logical(X),2);
% % % % % % % %         
% % % % % % % %         
% % % % % % % %         X = Surft.SurfData.vertices(:,2);
% % % % % % % %         X = X(temp);
% % % % % % % %         X(indz) =0;
% % % % % % % %         X = X(indminima,:);
% % % % % % % %         tempX = sum(X,2);
% % % % % % % %         Yo = tempX./sum(logical(X),2);
% % % % % % % %         
% % % % % % % %         X = Surft.SurfData.vertices(:,3);
% % % % % % % %         X = X(temp);
% % % % % % % %         X(indz) =0;
% % % % % % % %         X = X(indminima,:);
% % % % % % % %         tempX = sum(X,2);
% % % % % % % %         Zo = tempX./sum(logical(X),2);
% % % % % % % %         
% % % % % % % %         Surft.SurfData.vertices(indminima,:) = [Xo Yo Zo];
% % % % % % % %         
% % % % % % % %         X = Surft.Is;
% % % % % % % %         X = X(temp);
% % % % % % % %         X(indz) =0;
% % % % % % % %         X = X(indminima,:);
% % % % % % % %         tempX = sum(X,2);
% % % % % % % %         Cmean = tempX./sum(logical(X),2);
% % % % % % % %         
% % % % % % % %         
% % % % % % % %         %Surft.SurfData=smoothpatch(Surft.SurfData,1,8);
% % % % % % % %         %[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(Surft.SurfData.vertices,Surft.SurfData.faces);
% % % % % % % %         Surft.Is(indminima) = Cmean;
% % % % % % % %         
% % % % % % % %         distmap = -1*Surft.Is;
% % % % % % % %         temp1 = distmap(temp);
% % % % % % % %         temp1(indz) =0;
% % % % % % % %         t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
% % % % % % % %         indminima = find(sum(ismember(t,-1),2) - sum(logical(t),2) == 0);
% % % % % % % %     end
% % % % % % % %     if i == 1
% % % % % % % %         Surftemp.SurfData.vertices = Surft.SurfData.vertices;
% % % % % % % %         Surftemp.SurfData.faces = Surft.SurfData.faces;
% % % % % % % %     else
% % % % % % % %         Surftemp.SurfData.faces = [Surftemp.SurfData.faces; Surft.SurfData.faces + size(Surftemp.SurfData.vertices,1)];
% % % % % % % %         Surftemp.SurfData.vertices = [Surftemp.SurfData.vertices; Surft.SurfData.vertices];
% % % % % % % %         
% % % % % % % %     end
% % % % % % % %     
% % % % % % % %     Surf.SurfData.vertices(ind,:) = Surft.SurfData.vertices;
% % % % % % % %     
% % % % % % % %     
% % % % % % % %     
% % % % % % % % end




Surf.Is(indp,1) = 2; % Gyral regions are labeled with 2
sts = unique(nonzeros(Surf.Is));


ind = find(Surf.Is == 1);
%     temp = sum(ismember(Surf.SurfData.faces,ind),2);
%     indexbound = find(temp <3&temp>0);
%     indexbound = Surf.SurfData.faces(indexbound,:);
%     indexbound = unique(indexbound(:));
%     indexbound(find(ismember(indexbound,ind))) = [];


Surft = Surf;
Surft.Is = curvmap;

% Selecting only the faces in the sulci
%indneg = find(Surft.Is <0);
a = ismember(Surft.SurfData.faces,ind);
faces2rem = find(sum(a,2) ~= 3);
Surft.SurfData.faces(faces2rem,:) = [];
a = ismember(ind,Surft.SurfData.faces(:));
ind(find(a == 0)) = [];

% Reorganicing surface
[Surft] = Reorg_Surf(Surft);


distmap = -1*Surft.Is;
[Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

temp1 = distmap(temp);
temp1(indz) =0;
t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
indminima = find(sum(ismember(t,-1),2) - sum(logical(t),2) == 0);
%     hold on;
%     plot3(Surft.SurfData.vertices(indminima,1),Surft.SurfData.vertices(indminima,2),Surft.SurfData.vertices(indminima,3),'.y','MarkerSize',20);
Surfo = Surft;
cont = 0;
while (cont < 500)% ~isempty(indminima)||
    cont = cont + 1;
    temp(indz) = 1;
    X = Surft.SurfData.vertices(:,1);
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Xo = tempX./sum(logical(X),2);
    
    
    X = Surft.SurfData.vertices(:,2);
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Yo = tempX./sum(logical(X),2);
    
    X = Surft.SurfData.vertices(:,3);
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Zo = tempX./sum(logical(X),2);
    
    Surft.SurfData.vertices(indminima,:) = [Xo Yo Zo];
    
    X = Surft.Is;
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Cmean = tempX./sum(logical(X),2);
    
    
    %Surft.SurfData=smoothpatch(Surft.SurfData,1,8);
    %[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(Surft.SurfData.vertices,Surft.SurfData.faces);
    Surft.Is(indminima) = Cmean;
    
    distmap = -1*Surft.Is;
    temp1 = distmap(temp);
    temp1(indz) =0;
    t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
    indminima = find(sum(ismember(t,-1),2) - sum(logical(t),2) == 0);
end


Surf.SurfData.vertices(ind,:) = Surft.SurfData.vertices;







ind = find(Surf.Is == 2);
%     temp = sum(ismember(Surf.SurfData.faces,ind),2);
%     indexbound = find(temp <3&temp>0);
%     indexbound = Surf.SurfData.faces(indexbound,:);
%     indexbound = unique(indexbound(:));
%     indexbound(find(ismember(indexbound,ind))) = [];


Surft = Surf;
Surft.Is = curvmap;

% Selecting only the faces in the sulci
%indneg = find(Surft.Is <0);
a = ismember(Surft.SurfData.faces,ind);
faces2rem = find(sum(a,2) ~= 3);
Surft.SurfData.faces(faces2rem,:) = [];
a = ismember(ind,Surft.SurfData.faces(:));
ind(find(a == 0)) = [];

% Reorganicing surface
[Surft] = Reorg_Surf(Surft);


distmap = Surft.Is;
[Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

temp1 = distmap(temp);
temp1(indz) =0;
t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
indminima = find(sum(ismember(t,1),2) - sum(logical(t),2) == 0);

Surfo = Surft;
cont = 0;
while (cont < 500)% ~isempty(indminima)||
    cont = cont + 1;
    temp(indz) = 1;
    X = Surft.SurfData.vertices(:,1);
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Xo = tempX./sum(logical(X),2);
    
    
    X = Surft.SurfData.vertices(:,2);
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Yo = tempX./sum(logical(X),2);
    
    X = Surft.SurfData.vertices(:,3);
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Zo = tempX./sum(logical(X),2);
    
    Surft.SurfData.vertices(indminima,:) = [Xo Yo Zo];
    
    X = Surft.Is;
    X = X(temp);
    X(indz) =0;
    X = X(indminima,:);
    tempX = sum(X,2);
    Cmean = tempX./sum(logical(X),2);
    
    
    %Surft.SurfData=smoothpatch(Surft.SurfData,1,8);
    %[Umin,Umax,Cmin,Cmax,Cmean,Cgauss,Normal] = compute_curvature(Surft.SurfData.vertices,Surft.SurfData.faces);
    Surft.Is(indminima) = Cmean;
    
    distmap = -1*Surft.Is;
    temp1 = distmap(temp);
    temp1(indz) =0;
    t = sign(repmat(distmap,[1 size(temp1,2)]) - temp1);
    indminima = find(sum(ismember(t,1),2) - sum(logical(t),2) == 0);
end


Surf.SurfData.vertices(ind,:) = Surft.SurfData.vertices;























