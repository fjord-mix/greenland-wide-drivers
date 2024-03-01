function check_glacier_clustering(datasets,fjords_compilation)
    
if isfile(datasets.glaciers_file)
    glaciers = load(datasets.glaciers_file).glaciers;
else
    disp('Processed glaciers file not found. Pre-processing glaciers structure...')
    glaciers = initialise_glaciers(datasets);
end


glaciers_matrix = NaN([length(glaciers),3]); % qsg, d, gldepth
for i=1:length(glaciers)
    glaciers_matrix(i,1)=mean(glaciers(i).runoff.q,'omitnan');
    glaciers_matrix(i,2)=mean(glaciers(i).iceberg.discharge,'omitnan');
    glaciers_matrix(i,3)=abs(glaciers(i).gldepth);
    glaciers_matrix(i,4)=glaciers(i).regionID;
end

fjords_matrix = NaN([length(fjords_compilation),4]); %maxdepth, silldepth, length, width
for i=1:length(fjords_compilation)
    fjords_matrix(i,1) = abs(fjords_compilation(i).maxdepth);
    fjords_matrix(i,2) = abs(fjords_compilation(i).silldepth);
    fjords_matrix(i,3) = abs(fjords_compilation(i).length);
    fjords_matrix(i,4) = abs(fjords_compilation(i).width);
end

%% k-means cluster analysis
% n_clusters=6;
% [cidx,cmeans] = kmeans(glaciers_matrix(:,1:3),n_clusters,'dist','sqeuclidean','replicates',10,'display','final');
% figure;
% [silh,h] = silhouette(glaciers_matrix(:,1:3),cidx,'sqeuclidean');
% 
% lcolor=lines(n_clusters);
% markers={'o','+','x','v','square','pentagram','_'};
% figure;
% for i = 1:n_clusters
%     clust = find(cidx==i);
%     scatter3(glaciers_matrix(clust,1),glaciers_matrix(clust,2),glaciers_matrix(clust,3),'marker',markers{i},'color',lcolor(glaciers_matrix(clust,4),:));
%     hold on
%     endering
% scatter3(cmeans(:,1),cmeans(:,2),cmeans(:,3),'ko');
% scatter3(cmeans(:,1),cmeans(:,2),cmeans(:,3),'kx');
% hold off
% xlabel('Qsg');
% ylabel('D');
% zlabel('Zg');
% view(-137,10);
% grid on

%% DBSCAN clustering analysis - did not work? only 1 cluster always...
% minpts = 5; % needs to be 3+1 or 4+1 for glaciers or fjords, respectively
% 
% % determine the value for epsilon
% kD = pdist2(glaciers_matrix(:,1:3),glaciers_matrix(:,1:3),'euc','Smallest',minpts);
% plot(sort(kD(end,:)));
% title('k-distance graph'); xlabel(sprintf('Points sorted with %dth nearest distances',minpts)); ylabel('50th nearest distances')
% grid
% 
% epsilon=0.5; % epsilon could be 2 or 3
% 
% labels = dbscan(glaciers_matrix(:,1:3),epsilon,minpts);
% numGroups = length(unique(labels));
% figure;
% scatter3(glaciers_matrix(:,1),glaciers_matrix(:,2),glaciers_matrix(:,3),labels,hsv(numGroups));
% title('epsilon = 2 and minpts = 50')
% xlabel('Qsg'); ylabel('D'); zlabel('Zg');
% view(-137,10);
% grid

%% Hierarchical clustering
eucD = pdist(glaciers_matrix,'euclidean');
clustTree = linkage(eucD,'average');
cophenet(clustTree,eucD) % Large values indicate that the tree fits the distances well
figure; [h,nodes] = dendrogram(clustTree,0);

end