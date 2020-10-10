addpath(genpath('\\your\directory'));


listing = dir('Folder containing the WT tubule DGE files\')
figure(1);
stageout = zeros(1,length(listing)-2);
for i =3:length(listing)
   
name = listing(i).name;
name = ['Folder containing the WT tubule DGE files\' name];

t = readtable(name);

x = t.xcoord;
y = t.ycoord;
celltype = t.max_cell_type;
knn_idx = knnsearch([x y],[x y],'K',11)

M_out = zeros(max(celltype),max(celltype));

for j = 1:max(celltype)
    idx = find(celltype == j)
    N_out = zeros(1,max(celltype));
    for k =1:length(idx)
        neighbors = celltype(knn_idx(idx(k), 2:end));
        [N,edges] = histcounts(neighbors,[0.5:1:max(celltype)+0.5]);
        N_out = N + N_out;
    end
    M_out(j,:) = N_out;
end

subplot(10,11,i-2)
imagesc(M_out);

M_cell{i-2} = M_out;
stageout(i-2) =t.stage(1);


end

%diabetes

listing = dir('Folder containing the diabetic tubule DGE files\')
figure(2);
for i =3:length(listing)
   
name = listing(i).name;
name = ['Folder containing the diabetic tubule DGE files\' name];

t = readtable(name);

x = t.xcoord;
y = t.ycoord;
celltype = t.max_cell_type;
knn_idx = knnsearch([x y],[x y],'K',11)

M_out = zeros(max(celltype),max(celltype));
for j = 1:max(celltype)
    idx = find(celltype == j)
    N_out = zeros(1,max(celltype));
    for k =1:length(idx)
        neighbors = celltype(knn_idx(idx(k), 2:end));
        [N,edges] = histcounts(neighbors,[0.5:1:max(celltype)+0.5]);
        N_out = N + N_out;
    end
    M_out(j,:) = N_out;
end

subplot(14,14,i-2)
imagesc(M_out);

diabetesM_cell{i-2} = M_out;


end

%here we are dropping the 9th cell tyope
wtarray = zeros(45, length(M_cell))
for i = 1:length(M_cell)
    M_out = M_cell{i};
   
    m  = triu(true(9,9));
    
    if length(M_out) <9
        wtarray(:,i) = NaN(1,45)
        
    else
    v  = M_out(m).';
    wtarray(:,i) = v;
    end
end


%here we are dropping the 9th cell tyope
obarray = zeros(45, length(diabetesM_cell))
for i = 1:length(diabetesM_cell)
    M_out = diabetesM_cell{i};
   
    m  = triu(true(9,9));
    
    if length(M_out) <9
        obarray(:,i) = NaN(1,45)
        
    else
    v  = M_out(m).';
    obarray(:,i) = v;
    end
end


obarray = obarray(:,all(~isnan(obarray)));   % for nan - columns

stageout = stageout(all(~isnan(wtarray)));

wtarray = wtarray(:,all(~isnan(wtarray)));   % for nan - columns

%checking out some clustering
concat = [obarray wtarray];
normconcat = normc(concat);
[a_pca,b_pca,c_pca,d_pca,e_pca] = pca(normconcat','NumComponents',10);

y1 = tsne(b_pca);

figure;gscatter(y1(:,1),y1(:,2),[ones(1,length(obarray)), 2*ones(1,length(wtarray))]);title('TSNE by condition')
figure; gscatter(b_pca(:,1),b_pca(:,2),[ones(1,length(obarray)), 2*ones(1,length(wtarray))]);title('TSNE by condition')
xlabel('PC1')
ylabel('PC2')
axis([-.5 1 -.5 1])
yticks([-.5 0 .5 1])
printFig('outName','PC WT OB')

figure; gscatter(b_pca(length(obarray)+1:end,1),b_pca(length(obarray)+1:end,2),stageout);title('PCA by stage')
xlabel('PC1')
ylabel('PC2')
axis([-.5 1 -.5 1])
yticks([-.5 0 .5 1])
printFig('outName','PC WT Stage')

figure; gscatter(y1(length(obarray)+1:end,1),y1(length(obarray)+1:end,2),stageout);title('TSNE by stage')


[k,idx] = kmeans(b_pca,2,'Distance','correlation');
cgo = clustergram(idx,'RowPDist','correlation','ColumnPDist','correlation')
k_new = my_changem(k,1:2,cellfun(@str2num,get(cgo,'RowLabels')));
figure;gscatter(y1(:,1),y1(:,2),k_new);title('Tsne colored by cluster')
figure; gscatter(b_pca(:,1),b_pca(:,2),k_new);title('TSNE by condition')

figure; imagesc([normconcat(:,find(k_new ==1)) normconcat(:,find(k_new ==2)) normconcat(:,find(k_new ==3)) normconcat(:,find(k_new ==4))])



concat = [wtarray];
normconcat = normc(concat);
[a_pca,b_pca,c_pca,d_pca,e_pca] = pca(normconcat','NumComponents',10);

y1 = tsne(b_pca);

figure
scatter(y1(:,1),y1(:,2))



concat = [obarray];
normconcat = normc(concat);
[a_pca,b_pca,c_pca,d_pca,e_pca] = pca(normconcat','NumComponents',10);

y1 = tsne(b_pca);

figure
scatter(y1(:,1),y1(:,2))




%looking at some means
obmean = nanmean(normc(obarray)');
wtmean = nanmean(normc(wtarray)');

figure
bar([obmean;wtmean]')

%looking at some means

group1 = find(k_new==1);

obmean = nanmean(normc(obarray(:, group1(find(group1<=length(obarray)))))');
wtmean = nanmean(normc(wtarray(:, group1(find(group1>length(obarray))) - length(obarray) ))');

 out_array = obarray(:, group1(find(group1<=length(obarray))));
   out_array = normc(out_array);
   for j = 1:length(out_array)
        
        outmat = triu(ones(9));

        outmat(outmat > 0) = out_array(:,j);
        csvwrite(['ob_cluster1_' num2str(j) '.csv'], outmat)
   end
     out_array = wtarray(:, group1(find(group1>length(obarray))) - length(obarray) )
   out_array = normc(out_array);
   for j = 1:size(out_array,2)
        
        outmat = triu(ones(9));

        outmat(outmat > 0) = out_array(:,j);
        csvwrite(['wt_cluster1_' num2str(j) '.csv'], outmat)
   end




figure
bar([obmean;wtmean]')



obmeanmat = triu(ones(9));

obmeanmat(obmeanmat > 0) = obmean;

wtmeanmat = triu(ones(9));

wtmeanmat(wtmeanmat > 0) = wtmean;
figure
subplot(2,1,1)
imagesc(wtmeanmat);
title('wt')
subplot(2,1,2)
imagesc(obmeanmat);
title('ob')
figure();
imagesc(log2(obmeanmat./wtmeanmat));

group2 = find(k_new==2);

obmean = nanmean(normc(obarray(:, group2(find(group2<=length(obarray)))))');
wtmean = nanmean(normc(wtarray(:, group2(find(group2>length(obarray))) - length(obarray) ))');

   out_array = obarray(:, group2(find(group2<=length(obarray))));
   out_array = normc(out_array);
   meanob2 = zeros(9,9);
   for j = 1:length(out_array)
        
        outmat = triu(ones(9));

        outmat(outmat > 0) = out_array(:,j);
        csvwrite(['ob_cluster2_' num2str(j) '.csv'], outmat)
        meanob2 = outmat+meanob2;
   end
     out_array = wtarray(:, group2(find(group2>length(obarray))) - length(obarray) )
   out_array = normc(out_array);
   for j = 1:size(out_array,2)
        
        outmat = triu(ones(9));

        outmat(outmat > 0) = out_array(:,j);
        csvwrite(['wt_cluster2_' num2str(j) '.csv'], outmat)
   end



figure
bar([obmean;wtmean]')


obmeanmat = triu(ones(9));

obmeanmat(obmeanmat > 0) = obmean;

wtmeanmat = triu(ones(9));

wtmeanmat(wtmeanmat > 0) = wtmean;
figure
subplot(2,1,1)
imagesc(wtmeanmat);
title('wt')
subplot(2,1,2)
imagesc(obmeanmat);
title('ob')

figure();
imagesc(log2(obmeanmat./wtmeanmat));


