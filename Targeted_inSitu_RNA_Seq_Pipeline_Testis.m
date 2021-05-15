% Based on code from Zachary Chiang, Buenrostro Lab, Harvard University

clc;    % Clear the command window.
workspace;  % Make sure the workspace panel is showing.


%% define parameters

FOVs = [1]; 
FOV = 12; % this is the index of the FOV 
%% set up environment

tic
home_dir = '\\home\directory\containing\all\the\code';
data_dir = '\\data\directory\containing\all\the\images';
cd(home_dir);
addpath(genpath('scripts/')) %Add subsidiary code to the home direcotry
%% Load hyb DAPI

% all files have 1 series, 1 timepoint for now

series = 1;
timepoint = 1;

num_FOVs = size(FOVs,2);

%Here we use Round1 DAPI image as anchor for registration.
hyb_reader = bfGetReader(sprintf('%s/FOV%d/R1/*.tif',data_dir,FOV));

num_hyb_channels = hyb_reader.getSizeC;
hyb_dapi_channel = 5;

xlen(1) = hyb_reader.getSizeX; ylen(1) = hyb_reader.getSizeY; zlen(1) = hyb_reader.getSizeZ;

dapi_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
    
for z=1:zlen(1)
    dapi_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_dapi_channel,timepoint);
end



%% set testis fov bounds, manual for now

xmin = 1;
xmax = 2047;
ymin = 1;
ymax = 2047;
zmin = 1;
zmax = 17;

width = xmax-xmin;
height = ymax-ymin;
depth = zmax;
FOV_bounds = [xmin ymin zmin width height depth]; 



%% make FOV directories

for i=1:num_FOVs

    FOV_dir = sprintf('%s/FOV%d/Analysis',data_dir,FOV);
    mkdir(FOV_dir)

    fig_dir = sprintf('%s/figure',FOV_dir);
    mkdir(fig_dir)
    
    processed_dir = sprintf('%s/processed',FOV_dir);
    mkdir(processed_dir)    
    
    offset_dir = sprintf('%s/processed/offset',FOV_dir);
    mkdir(offset_dir)

    reg_dir = sprintf('%s/reg',FOV_dir);
    if ~exist(reg_dir, 'dir') mkdir(reg_dir), end
    
end


for i=1:num_FOVs
    dlmwrite(sprintf('%s/processed/bounds.txt',FOV_dir),FOV_bounds(i,:));
end



%% save hyb images

for i=1:num_FOVs
    
    write_3d_tif(sprintf('%s/processed/offset/hyb_dapi.tif',FOV_dir),imcrop_xyz(dapi_stacks{1},FOV_bounds(i,:)));
    
end


disp(sprintf('%s: wrote hyb and assoc. images',sec2time(toc)))

%% load data

% set channel and cycle info

dapi_channel = 5;
num_cycles = 5; %number of sequencing rounds

dapi_offset_stacks{1} = dapi_stacks{1};

for cycle=1:num_cycles
    
    % Read in sequencing stack from each round.
    seq_reader = bfGetReader(sprintf('%s/FOV%d/R%d/*.tif',data_dir,FOV,cycle));
    
    xlen(cycle+1) = seq_reader.getSizeX; ylen(cycle+1) = seq_reader.getSizeY; zlen(cycle+1) = seq_reader.getSizeZ; 

    max_dims = [max(xlen) max(ylen) max(zlen)];
    num_channels = seq_reader.getSizeC;
    
    % load DAPI
    
    dapi_stacks{cycle+1} = zeros(xlen(cycle+1),ylen(cycle+1),zlen(cycle+1),'uint16');
    for z=1:zlen(cycle+1)
        dapi_stacks{cycle+1}(:,:,z) = readPlane(seq_reader,series,z,dapi_channel,timepoint);
    end
    
    
    disp(sprintf('%s: loaded DAPI for cycle %d',sec2time(toc),cycle))
    
    % offset DAPI

    [dapi_offset_stacks{cycle+1} offsets{cycle+1}] = get_offset_xyz(dapi_offset_stacks{1},dapi_stacks{cycle+1},max_dims);

    
    disp(sprintf('%s: offset DAPI for cycle %d',sec2time(toc),cycle))
    
    % write DAPI
    
    for i=1:num_FOVs
        write_3d_tif(sprintf('%s/R%02d_dapi.tif',offset_dir,cycle),imcrop_xyz(dapi_offset_stacks{cycle+1},FOV_bounds(i,:)));
    end
    
    disp(sprintf('%s: wrote DAPI for cycle %d',sec2time(toc),cycle))
    
    % loop through seq channels
    
    channels = 1:num_channels; channels(channels==dapi_channel) = [];
    for channel=[channels]
        
        % load seq
        
        seq_stack = zeros(xlen(cycle+1),ylen(cycle+1),zlen(cycle+1),'uint16');
        for z=1:zlen(cycle+1)
            seq_stack(:,:,z) = readPlane(seq_reader,series,z,channel,timepoint);
        end
        
        disp(sprintf('%s: loaded cycle %d, channel %d',sec2time(toc),cycle, channel))
        
        % offset seq
        
        offset_seq_stack = apply_offset_xyz(seq_stack,offsets{cycle+1},max_dims);
        figure; imshowpair(capImage(max(dapi_offset_stacks{2},[],3),99,'prc'),capImage(max(offset_seq_stack,[],3),99,'prc'))
        
       
        
        disp(sprintf('%s: offset cycle %d, channel %d',sec2time(toc),cycle, channel))
        
        % write seq
        
        for i=1:num_FOVs
            write_3d_tif(sprintf('%s/R%02d_Channel%02d.tif',offset_dir,cycle,channel),imcrop_xyz(offset_seq_stack,FOV_bounds(i,:)));
        end
        
        disp(sprintf('%s: wrote cycle %d, channel %d',sec2time(toc),cycle, channel))
        
    end
    
end


min_overlap = 0.5;
num_channels = 4;
num_cycles = 5;

bounds = dlmread(sprintf('%s/processed/bounds.txt',FOV_dir));

stack = zeros(bounds(5),bounds(4),bounds(6),num_channels,num_cycles,'uint16');

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/processed/offset/R%02d_Channel%02d.tif',FOV_dir,cycle,channel);
        stack(:,:,:,channel,cycle) = read_3d_tif(filename,bounds(5),bounds(4),bounds(6));
    end
end

reg_stack2 = stack;

%% Deconvolve with Gaussian filter

deconv_stack = zeros(size(reg_stack2));

for cycle=1:num_cycles
    for channel=1:num_channels
        for z=1:size(deconv_stack,3)
            high_pass_filter = imgaussfilt(reg_stack2(:,:,z,channel,cycle),2);
            high_pass_image = reg_stack2(:,:,z,channel,cycle) - high_pass_filter;
     
            deconv_stack(:,:,z,channel,cycle) = high_pass_image;
        end
    end
end

deconv_stack(deconv_stack<0) = 0;

disp(sprintf('%s: Deconvolved images with Gaussian filter',sec2time(toc)))

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/decon_R%02d_Channel%02d.tif',reg_dir,cycle,channel);
        write_3d_tif(filename, deconv_stack(:,:,:,channel,cycle));
    end
end


disp(sprintf('%s: Saved deconvolved stacks',sec2time(toc)))

%% Load DAPI stack from round 1. 
reg_dapi_stack = read_3d_tif(sprintf('%s/hyb_dapi.tif',offset_dir),bounds(5),bounds(4),bounds(6));

disp(sprintf('%s: Loaded normalized stacks',sec2time(toc)))
%% Peak calling

all_peaks = [];
channel_peaks = {};
cycle = 1;

for channel=1:num_channels
        
    [Maxima,MaxPos,Minima,MinPos]=MinimaMaxima3D(deconv_stack(:,:,:,channel,cycle),1,0);
    
    channel_peaks_pos{channel} = MaxPos;
    channel_peaks_max{channel} = Maxima;

    disp(sprintf('%s: Found %d 3D peaks in cycle %d, channel %d',sec2time(toc),length(channel_peaks_pos{channel}),cycle,channel));
          
end

peaks_table = table;
peaks_table.x = [channel_peaks_pos{1}(:,1); channel_peaks_pos{2}(:,1); channel_peaks_pos{3}(:,1); channel_peaks_pos{4}(:,1)];
peaks_table.y = [channel_peaks_pos{1}(:,2); channel_peaks_pos{2}(:,2); channel_peaks_pos{3}(:,2); channel_peaks_pos{4}(:,2)];
peaks_table.z = [channel_peaks_pos{1}(:,3); channel_peaks_pos{2}(:,3); channel_peaks_pos{3}(:,3); channel_peaks_pos{4}(:,3)];
peaks_table.val = [channel_peaks_max{1}; channel_peaks_max{2}; channel_peaks_max{3}; channel_peaks_max{4}];
peaks_table.channel = [repmat(1,length(channel_peaks_max{1}),1); repmat(2,length(channel_peaks_max{2}),1); repmat(3,length(channel_peaks_max{3}),1); repmat(4,length(channel_peaks_max{4}),1);];

writetable(peaks_table,sprintf('%s/peaks.txt',FOV_dir))

%% Read peaks

head(peaks_table)

ch1 = (peaks_table(peaks_table.val == '1',:));


%% Peak call video

out_dir = sprintf('%s/peak_calls',fig_dir);
if ~exist(out_dir, 'dir') mkdir(out_dir), end

video = VideoWriter(sprintf('%s/FOV%d',out_dir,FOV),'MPEG-4');
video.FrameRate = 1;
open(video);

z_buffer = 2;

%Need to plot the peak distribution to determine the thresh
thresh = [32 32 16 32];

for z=1:size(deconv_stack,3)
    
    f = figure('visible','off');
    p = tight_subplot(1,4,[0.001 0.001],[0.001 0.001],[0.001 0.001]);
    
    for channel=1:num_channels
    
        peaks = peaks_table{peaks_table.channel == channel & peaks_table.val > thresh(channel),1:3};
      
        z_channel_peaks = peaks(ismember(peaks(:,3),(z+[-z_buffer:+z_buffer])),:);

        axes(p(channel));
        imshowpair(capImage(deconv_stack(:,:,z,channel,1),99,'prc'),capImage(reg_dapi_stack(:,:,z),95,'prc')); hold on;
        plot(z_channel_peaks(:,2),z_channel_peaks(:,1),'Marker','.','MarkerEdgeColor','r','LineStyle','none');
        title(sprintf('Peak calls (n = %d)',size(peaks,1)))
        
    end
   
    f.Position = [0 0 size(deconv_stack,2).*0.5.*(4) size(deconv_stack,1).*0.75];    
    writeVideo(video,getframe(gcf))
    
end
close(video); %close the file
close all

disp(sprintf('%s: Saved peak call video',sec2time(toc)));

%this part was for trying to do fine registration using peaks with high intensities. 
thresh = [200 200 200 200];
all_peaks_reg = [];
for channel=1:num_channels
    all_peaks_reg = [all_peaks_reg; peaks_table{peaks_table.channel == channel & peaks_table.val > thresh(channel),1:3}];
end

norm_stack_decon = zeros(size(reg_stack2));

for cycle=1:num_cycles
  
    tmp_stack = reshape(deconv_stack(:,:,:,:,cycle),bounds(5)*bounds(4)*bounds(6),num_channels);
    quantile_norm = quantilenorm(double(tmp_stack));
    norm_stack_decon(:,:,:,:,cycle) = reshape(quantile_norm,bounds(5),bounds(4),bounds(6),num_channels);
    
end


max_stacks = squeeze(max(norm_stack_decon, [],3));
max_cyc = squeeze(max(max_stacks,[],3));
filename = sprintf('%s/before_fine_reg.tif',reg_dir);
write_3d_tif(filename, max_cyc);


%fine registration by imrecorr in 2d. 
refImg = im2double(mat2gray(capImage(max_cyc(:,:,1),99,'prc')));
postfine_deconv = zeros(size(deconv_stack));
postMax = zeros(size(max_cyc));
postMax(:,:,1) = refImg;

for cycle = 2:num_cycles
moving = im2double(mat2gray(capImage(max_cyc(:,:,cycle),99,'prc')));

tform = imregcorr(moving,refImg);

Rfixed = imref2d(size(refImg));
movingReg = imwarp(moving,tform,'OutputView',Rfixed);

postMax(:,:,cycle) = movingReg;
for channel = 1:4
    for z = 1:size(deconv_stack,3)
        postfine_deconv(:,:,z,channel,cycle) = imwarp(deconv_stack(:,:,z,channel,cycle),tform,'OutputView',Rfixed);
    end
end

end

filename = sprintf('%s/after_fine_reg.tif',reg_dir);
write_3d_tif(filename, uint16(postMax*65536));
postfine_deconv(:,:,:,:,1) = deconv_stack(:,:,:,:,1);


thresh = [32 32 16 32];
all_peaks = [];
for channel=1:num_channels
    all_peaks = [all_peaks; peaks_table{peaks_table.channel == channel & peaks_table.val > thresh(channel),1:3}];
end

pad = 2;
maxz =2;
mat = zeros(size(all_peaks,1),num_channels,num_cycles);
purity = zeros(size(all_peaks,1),num_cycles);
consensus = zeros(size(all_peaks,1),num_cycles);

xleft = max(all_peaks(:,1)-pad,1); xright = min(all_peaks(:,1)+pad,size(deconv_stack,1));
yleft = max(all_peaks(:,2)-pad,1); yright = min(all_peaks(:,2)+pad,size(deconv_stack,2));
zleft = max(all_peaks(:,3)-maxz,1); zright = min(all_peaks(:,3)+maxz,size(deconv_stack,3));

for i=1:size(all_peaks,1); disp(i)
    
    peak = all_peaks(i,:);
    peak_mat = squeeze(sum(sum(sum(postfine_deconv(xleft(i):xright(i),yleft(i):yright(i),zleft(i):zright(i),:,:),1),2),3));  
    mat(i,:,:) = peak_mat;
    [tmp consensus(i,:)] = max(peak_mat,[],1);
    
end


%quantilenorm of peak intensities 

figure
histogram(mat(:,1,1));
hold on
histogram(mat(:,2,1));
histogram(mat(:,3,1));
histogram(mat(:,4,1));

Nmat = zeros(size(mat));
for j = 1:5
Nmattem = quantilenorm(squeeze(mat(:,:,j)));
for i = 1:4
Nmat(:,i,j) = Nmattem(:,i);
end
end

figure
histogram(Nmat(:,1,1))
hold on
histogram(Nmat(:,2,1))
histogram(Nmat(:,3,1))
histogram(Nmat(:,4,1))

 for i=1:size(all_peaks,1); disp(i)
    
    peak_mat =  squeeze(Nmat(i,:,:));
    purity(i,:) = max(peak_mat.^2,[],1)./sum(peak_mat.^2,1);
    [tmp consensus(i,:)] = max(peak_mat,[],1);
    
 end   

figure
histogram(purity(:,1))
hold on
histogram(purity(:,2))
histogram(purity(:,3))
histogram(purity(:,4))
histogram(purity(:,5))
legend('purity r1', 'purity r2','purity r3','purity r4','purity r5')


thresh = 0.6;
purity_score = sum(purity>thresh,2);

figure; histogram(sum(purity>thresh,2))
xlabel('rounds w/ purity score > 0.6'); xticks([0:4])
ylabel('spots')

video = VideoWriter(sprintf('%s/Purity_video_FOV%d',out_dir, FOV),'MPEG-4');
video.FrameRate = 1;
open(video);

z_buffer = 2;

for z=1:size(postfine_deconv,3)
    
    f = figure;
    
    z_peaks = all_peaks(ismember(all_peaks(:,3),(z+[-z_buffer:+z_buffer])),:);
    z_purity = purity_score(ismember(all_peaks(:,3),(z+[-z_buffer:+z_buffer])));
        
    imshowpair(capImage(max(postfine_deconv(:,:,z,channel,1),[],4),99,'prc'),capImage(reg_dapi_stack(:,:,z),95,'prc')); hold on;
    scatter(z_peaks(:,2),z_peaks(:,1),15,z_purity,'filled');
    colormap(jet); caxis([0 4])
  
    writeVideo(video,getframe(gcf))
    
end
close(video); %close the file
close all


%% Add in barcode info

valid_barcodes = [...
    "3003131","0113101","3030232","1301121",...
    "3312133","2230321","0210010","2013313","1130013","1200211","0012200","2121320","3102233",...
    "1313021","1002012","0021100","2302322",...
    "0330103","3123032","1233120","1332212","3232033",...
    ]';

tmp = char(valid_barcodes);

%The barcodes are designed based on the order of SeqN-1_1st, SeqN-1_2nd, SeqN_1ST, SeqN_2nd, SeqN-2_2nd,
%but our imaging order was SeqN, SeqN-1, and SeqN-2, hence the reorder
%below.
valid_barcodes = string([tmp(:,3) tmp(:,4) tmp(:,1) tmp(:,2) tmp(:,5)]);

genes = [...
    "Cyp11a1","Adgre1","Vwf","Acta2",...
    "Id4","Gfra1","Sox9","Stra8","Kit","Piwil1","Acrv1","Prm1","X1110018N20Rik",...
    "Bsdc1","Ccdc87","Hyal4","Kcnip4",...
    "Noxred1","Slc12a5","Smim23","Tmem183a","Txndc17",...
    ];

[sort_valid_barcodes sort_order] = sort(valid_barcodes);
sort_genes = genes(sort_order);

sel_consensus = consensus(purity_score>=4,:);

%This order is flipped b/c the seq channels we used is at the order of FITC, Cy3, TxRed
%and Cy5
sel_consensus_flipped = abs(sel_consensus -4)+1;
barcodes = strrep(string(num2str(4-sel_consensus_flipped)),' ','');
[C, ia, ic] = unique(barcodes);
counts = accumarray(ic,1);

on_target = ismember(C,sort_valid_barcodes);
on_target_prc = sum(counts(on_target))./sum(counts);
on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

figure;
hold on
for i = 1:length(C)
    h=bar(i,counts(i));
    if on_target(i) == 1
        set(h,'FaceColor','r');
    else
        set(h,'FaceColor','b');
    end
end
hold off

xlim([1 size(C,1)])

xticks(find(on_target))
xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
set(gca,'fontsize',8)

xtickangle(90)
title(sprintf("stringent on target %%: %.02f",on_target_prc*100))

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 5];

saveas(fig,sprintf('%s/Stringent_on_target_FOV%d.png',fig_dir,FOV));

genes_stringent = on_target_genes;
counts_stringent = counts(on_target);


sel_consensus = consensus(purity_score>=3,:);
sel_consensus_flipped = abs(sel_consensus -4)+1;
barcodes = strrep(string(num2str(4-sel_consensus_flipped)),' ','');
%barcodes = strrep(string(num2str(sel_consensus-1)),' ','');
[C, ia, ic] = unique(barcodes);
counts = accumarray(ic,1);

on_target = ismember(C,sort_valid_barcodes);
on_target_prc = sum(counts(on_target))./sum(counts);
on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

figure;
hold on
for i = 1:length(C)
    h=bar(i,counts(i));
    if on_target(i) == 1
        set(h,'FaceColor','r');
    else
        set(h,'FaceColor','b');
    end
end
hold off

xlim([1 size(C,1)])

xticks(find(on_target))
xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
set(gca,'fontsize',8)

xtickangle(90)
title(sprintf("lenient on target %%: %.02f",on_target_prc*100))

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 5];

saveas(fig,sprintf('%s/Linient_on_target_FOV%d.png',fig_dir,FOV));

genes_lenient = on_target_genes;
counts_lenient = counts(on_target);

%% Save gene counts

gene_counts = table;
gene_counts.gene = sort_genes';

gene_counts.counts_stringent = zeros(size(gene_counts,1),1);
gene_counts.counts_stringent(ismember(sort_genes',genes_stringent')) = counts_stringent;

gene_counts.counts_lenient = zeros(size(gene_counts,1),1);
gene_counts.counts_lenient(ismember(sort_genes',genes_lenient')) = counts_lenient;

writetable(gene_counts,sprintf('%s/gene_counts_FOV%d.txt',fig_dir,FOV))


sel_purity = purity_score(purity_score>=3);
on_target_purity = sel_purity(ismember(barcodes,sort_valid_barcodes));

sel_peaks = all_peaks(purity_score>=3,:);
on_target_peaks = sel_peaks(ismember(barcodes,sort_valid_barcodes),:);

barcode_dict = strings(4000,1);
barcode_dict(str2double(sort_valid_barcodes)) = sort_genes';
on_target_genes = barcode_dict(str2double(barcodes(ismember(barcodes,sort_valid_barcodes))));

gene_xyz_table = table;
gene_xyz_table.gene = on_target_genes;
gene_xyz_table.x = on_target_peaks(:,1);
gene_xyz_table.y = on_target_peaks(:,2);
gene_xyz_table.z = on_target_peaks(:,3);
gene_xyz_table.purity = on_target_purity;

writetable(gene_xyz_table,sprintf('%s/gene_xyz_table_FOV%d.txt',fig_dir,FOV))

%% Load segmentation use DAPI

%This file is from the 3D segmentation pipeline in CellProfiler using
%Round1 DAPI stack (i.e., hyb_dapi.tiff). The output from the CellProfiler
%needs to be opened in ImageJ first and be re-saved as a tiff file.
stain_seg_stack = read_3d_tif(sprintf('%s/FOV%d/Segmentation/hyb_dapi_testis_FOV%d_segmentedNuclei_5_ws_fiji.tif',data_dir,FOV,FOV),FOV_bounds(5),FOV_bounds(4),FOV_bounds(6));

seg_dapi_stack = read_3d_tif(sprintf('%s/hyb_dapi.tif',offset_dir),FOV_bounds(5),FOV_bounds(4),FOV_bounds(6));

max_dims = [FOV_bounds(4), FOV_bounds(5), FOV_bounds(6)];

offset_seg_stack = stain_seg_stack;
offset_seg_dapi_stack = seg_dapi_stack;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_3d_tif(sprintf('%s/seg.tif',offset_dir),offset_seg_stack);
write_3d_tif(sprintf('%s/offset_seg_dapi.tif',offset_dir),offset_seg_dapi_stack);


disp(sprintf('%s: Loaded stack',sec2time(toc)))


seg_stack = read_3d_tif(sprintf('%s/seg.tif',offset_dir),bounds(5),bounds(4),bounds(6));

%% Get assignments

sel = gene_xyz_table.purity >= 4; % 3 pr 4

gene_idx = sub2ind(size(seg_stack),gene_xyz_table.x(sel),gene_xyz_table.y(sel),gene_xyz_table.z(sel));
cell_assignment = seg_stack(gene_idx);

%% transcripts per cell

assigned = cell_assignment(cell_assignment ~= 0);
[C ia ic] = unique(assigned);
counts = accumarray(ic,1);

num_cells = max(seg_stack(:));
genes_per_cell = zeros(num_cells,1);
genes_per_cell(C) = counts;

fig = figure; hold on;
for i = 1:length(genes_per_cell)
    h=bar(i,genes_per_cell(i));
end

xlim([0 num_cells])
xlabel('cell index'); ylabel('# transcripts')
title(sprintf('%.2f%% of on-target transcripts assigned',sum(cell_assignment>0)./size(cell_assignment,1)*100))

saveas(fig,sprintf('%s/transcripts_per_cell_FOV%d.png',fig_dir,FOV));

%% create matrix

sel = gene_xyz_table.purity >= 4; 
num_genes = size(valid_barcodes,1);
gene_xyz_table.cell_assignment = zeros(size(gene_xyz_table,1),1);
gene_xyz_table.cell_assignment(sel) = cell_assignment;

cell_counts_mat = zeros(num_cells,num_genes);

for i=1:num_genes
    
    sel_cells = gene_xyz_table.cell_assignment(gene_xyz_table.gene == sort_genes(:,i));
    [C ia ic] = unique(sel_cells(sel_cells>0));
    counts = accumarray(ic,1);
    cell_counts_mat(C,i) = counts;
end

cell_counts_table = array2table(cell_counts_mat);
cell_counts_table.Properties.VariableNames = cellstr(sort_genes');
cell_counts_table.cell_index = (1:size(cell_counts_table,1))';

cell_info = regionprops(seg_stack,'Centroid','Area');
centroids = reshape([cell_info.Centroid]',3,[])';
cell_counts_table.x = centroids(:,2)*.17;
cell_counts_table.y = centroids(:,1)*.17;
cell_counts_table.z = centroids(:,3)*.4;
cell_counts_table.area = [cell_info.Area]';

cell_counts_table = [cell_counts_table(:,end-4:end) cell_counts_table(:,1:end-5)]; 

writetable(cell_counts_table,sprintf('%s/cell_counts_table_FOV%d.txt',fig_dir,FOV));

%% visualiztion

video = VideoWriter(sprintf('%s/FOV%d_cell_assignments',fig_dir,FOV),'MPEG-4');
video.FrameRate = 1;
open(video);
colors = distinguishable_colors(700); colors(4,:) = [];
zbuf = 1;
B= {}
fig = figure; 
for z=1:size(reg_dapi_stack,3)
    
    imshow(capImage(reg_dapi_stack(:,:,z),99,'prc'),[]); hold on;
    
    z_peaks = gene_xyz_table(any(gene_xyz_table.z == z+[-zbuf:zbuf],2) & gene_xyz_table.cell_assignment > 0,:);
    scatter(z_peaks.y,z_peaks.x,15,colors(z_peaks.cell_assignment,:),'filled');

    for i=1:num_cells
        [tmp L] = bwboundaries(seg_stack(:,:,z)==i,8,'noholes');
        if size(tmp,1)>0
            B{i} = tmp{1};
        end
    end

    for k = 1:length(B)
        if size(B{k},1)>0
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1),'Color',colors(round(k/2),:), 'LineWidth', 1); hold on;
        end
    end
    
    fig.Units = 'pixels';
    fig.Position = [0 0 size(reg_dapi_stack,2) size(reg_dapi_stack,1)];    
    
    writeVideo(video,getframe(gcf))
end

close(video)

%% visualize gene positions
 
for list = 1:length(genes)
    sel_gene = genes{list};
    sel_exp = cell_counts_table{:,find(sort_genes==sel_gene)+5};

    figure; 
    scatter3(cell_counts_table.x,cell_counts_table.y,cell_counts_table.z,50,sel_exp,'filled'); hold on;
    colorbar; %caxis([prctile(sel_exp,10) prctile(sel_exp,90)])
    title(sel_gene)

    %scatter3(gene_xyz_table.x(sel)*.17,gene_xyz_table.y(sel)*.17,gene_xyz_table.z(sel)*.4,10,sel_exp(gene_xyz_table.cell_assignment(sel)),'.');

    axis equal;
    saveas(gcf,sprintf('%s/Fig_FOV%d_%s',fig_dir,FOV, sel_gene),'epsc');
end
