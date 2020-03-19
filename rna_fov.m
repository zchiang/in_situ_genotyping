function[] = rna_exp(exp_id,fov)
%% Set manual parameters, initialize environment, create directories

% ------ Set manual parameters (comment out if running as function) ------ %

exp_id = 'italy_solid';
fov = 10;

% ------------------------------------------------------------------------ %

tic
clearvars -except exp_id fov
addpath(genpath('X:\Confocal Data June2018-\Zack\in_situ\rna\matlab_functions\'));
load(sprintf('%s/exp_parameters.mat',exp_id));
colors = distinguishable_colors(250);
disp(sprintf('%s: Started analysis of experiment %s, FOV %d',sec2time(toc),exp_id,fov));

% Create directories

fov_dir = sprintf('%s\\fov%03d',exp_dir,fov);
if ~exist(fov_dir, 'dir') mkdir(fov_dir), end

reg_dir = sprintf('%s\\registration',fov_dir);
if ~exist(reg_dir, 'dir') mkdir(reg_dir), end

disp(sprintf('%s: Output will be saved to %s',sec2time(toc),fov_dir));

%% Load all cycles and draq5 stain

stack = zeros(fov_xlen,fov_ylen,fov_zlen,num_channels,num_cycles);

for cycle=1:(length(exp_design.cycle)-1)
    
    % Read cycle parameters
    
    nd2 = exp_design.nd2{cycle};
    t = exp_design.t(cycle);
    
    reader = bfGetReader(sprintf('%s/%s', nd2_dir, nd2));
    disp(sprintf('%s: Processing cycle %d from file %s, t=%d',sec2time(toc),cycle,nd2,t));
    
    % Load all images
    
    for channel=1:num_channels
        for z=1:fov_zlen
            stack(:,:,z,channel,cycle) = readPlane(reader,fov,z,channel,t);
        end
    end
end

disp(sprintf('%s: Loaded all cycles for FOV %d',sec2time(toc),fov));

% Load draq5_stain

stain_stack = zeros(fov_xlen,fov_ylen,fov_zlen,stain_num_channels);

nd2 = exp_design.nd2{cycle+1};
t = exp_design.t(cycle+1);

reader = bfGetReader(sprintf('%s/%s',nd2_dir,nd2));
disp(sprintf('%s: Processing draq5 stain from file %s, t=%d',sec2time(toc),nd2,t));
    
for channel=1:stain_num_channels
    for z=1:fov_zlen
         stain_stack(:,:,z,channel) = readPlane(reader,fov,z,channel,1);
    end
end

disp(sprintf('%s: Loaded draq5 stain for FOV %d',sec2time(toc),fov));

%% Flatten images and register

flat_stack = squeeze(max(stack,[],3));
flat_stain_stack = squeeze(max(stain_stack,[],3));

reg_stack = zeros(size(flat_stack));
%reg_stack(:,:,:,1) = flat_stack(:,:,:,1);
%reg_stack(:,:,1,:) = flat_stack(:,:,1,:); % color correction is not robust for WGA channel

cy_offsets = zeros(num_cycles,2);
cc_offsets = zeros(num_channels,num_cycles,2);

bead_thresh = 20000;
bead_size = [10 10];
registration_channel = 2;

% find registration beads

[beads bead_locs] = find_beads_xy(flat_stain_stack(:,:,2),bead_thresh,bead_size);
disp(sprintf('%s: Identified %d beads',sec2time(toc),length(beads)));

%figure; imshow(flat_stain_stack(:,:,2)>bead_thresh,[]); hold on;
%plot(beads(:,1),beads(:,2),'rx')

for cycle=1:num_cycles

    [reg_stack(:,:,:,cycle) cy_offsets(cycle,:)] = fov_offset_xy(flat_stain_stack,flat_stack(:,:,:,cycle),1);
    disp(sprintf('%s: Applied offset of %d (x) by %d (y) to all channels in cycle %d',sec2time(toc),cy_offsets(cycle,1),cy_offsets(cycle,2),cycle));

    % calculate color correction
    
    for channel=2:num_channels
        [reg_stack(:,:,channel,cycle) cc_offsets(channel,cycle,:)] = calc_color_correction_xy(flat_stain_stack(:,:,2),reg_stack(:,:,channel,cycle),beads,bead_locs,bead_size);
        disp(sprintf('%s: Applied offset of %.02f (x) by %.02f (y) to channels %d, cycle %d',sec2time(toc),cc_offsets(channel,cycle,1),cc_offsets(channel,cycle,2),channel,cycle));
    end
     
end

%%
%figure; imshowpair(cap_image(max(reg_stack(:,:,2:5,1),[],3),10000),cap_image(flat_stain_stack(:,:,2),10000))
figure; imshowpair(cap_image(max(norm_stack(:,:,2:5,1),[],3),5000),cap_image(max(norm_stack(:,:,2:5,2),[],3),5000))
figure; imshowpair(cap_image(max(flat_stack(:,:,2:5,1),[],3),5000),cap_image(max(flat_stack(:,:,2:5,2),[],3),5000))

%% Two-pass Gaussian filter

deconv_stack = zeros(size(reg_stack));
deconv_stack(:,:,1,:) = reg_stack(:,:,1,:);

for cycle=1:num_cycles
    for channel=2:num_channels
        high_pass_filter = imgaussfilt(reg_stack(:,:,channel,cycle),4);
        high_pass_image = reg_stack(:,:,channel,cycle).*2 - high_pass_filter;
        deconv_stack(:,:,channel,cycle) = imgaussfilt(high_pass_image,1);
    end
end

deconv_stack(deconv_stack<0) = 0;

disp(sprintf('%s: Processed images with two-pass Gaussian filter',sec2time(toc)))

%% Fine registration

reg2_stack = zeros(size(deconv_stack));
reg2_stack(:,:,1,:) = deconv_stack(:,:,1,:);
reg2_stack(:,:,:,1) = deconv_stack(:,:,:,1);

metric = registration.metric.MeanSquares;
optimizer = registration.optimizer.RegularStepGradientDescent;
optimizer.MaximumIterations = 50;

for cycle=2:num_cycles
  
    [imreg_new, ~, tform_new] = imregister2(max(deconv_stack(:,:,2:5,cycle),[],3), max(deconv_stack(:,:,2:5,cycle-1),[],3), 'affine', optimizer, metric, 'PyramidLevels', 1);    
    for channel=1:num_channels
        reg2_stack(:,:,channel,cycle)=imwarp(deconv_stack(:,:,channel,cycle), tform_new, 'outputView', imref2d([2048 2048]));
    end
    
    disp(sprintf('%s: Registered cycle %d',sec2time(toc),cycle));
    
end

%% Cap high pixel values

cap = {};
cap_prctile = 99.95;
cap_stack = zeros(size(reg2_stack));
cap_stack(:,:,1,:) = reg2_stack(:,:,1,:);

for cycle=1:num_cycles
    for channel=2:num_channels
        cap{channel,cycle} = prctile(reshape(reg2_stack(:,:,channel,cycle),[],1),cap_prctile);
        tmp_stack = reg2_stack(:,:,channel,cycle);
        tmp_stack(tmp_stack>cap{channel,cycle}) = cap{channel,cycle};
        cap_stack(:,:,channel,cycle) = tmp_stack;
    end
end

disp(sprintf('%s: Capped high pixel values',sec2time(toc)));

%% Crop and quantile normalization

crop_bounds = [max([cy_offsets(:,1)' 1]) max([cy_offsets(:,2)' 1]) min(fov_xlen-max(cy_offsets(:,1))+[0 min(cy_offsets(:,1))]) min(fov_ylen-max(cy_offsets(:,2))+[0 min(cy_offsets(:,2))])];
crop_stack = imcrop_xy(cap_stack,crop_bounds);
crop_stain_stack = imcrop_xy(flat_stain_stack,crop_bounds);
dlmwrite(sprintf('%s/crop_bounds.txt',fov_dir),crop_bounds);

disp(sprintf('%s: Cropped FOV to size of %d (x) by %d (y)',sec2time(toc),crop_bounds(1,3),crop_bounds(1,4)));

norm_stack = zeros(size(crop_stack));
norm_stack(:,:,1,:) = crop_stack(:,:,1,:);

for cycle=1:num_cycles 
    tmp_stack = reshape(crop_stack(:,:,2:5,cycle),crop_bounds(1,4)*crop_bounds(1,3),num_channels-1);
    quantile_norm = quantilenorm(tmp_stack);
    norm_stack(:,:,2:5,cycle) = reshape(quantile_norm,crop_bounds(1,4),crop_bounds(1,3),num_channels-1);  
end

disp(sprintf('%s: Performed quantile normaliation on channels 2-5',sec2time(toc)))

%% Save visualization

figure('visible','off');
p = tight_subplot(1,num_cycles-1,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles-1
    axes(p(cycle)); imshowpair(max(norm_stack(:,:,2:5,cycle),[],3),max(norm_stack(:,:,2:5,cycle+1),[],3));
    title(sprintf('Cycle %d-%d',cycle,cycle+1))
end

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles-1)*8*size(norm_stack,2)/max(size(norm_stack,1),size(norm_stack,2)) 9*size(norm_stack,1)/max(size(norm_stack,1),size(norm_stack,2))];

fig_dir = sprintf('%s/figures/cycle_overlap',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/cycle_overlap.png',fov_dir));
saveas(fig,sprintf('%s/fov%03d.png',fig_dir,fov));

disp(sprintf('%s: Saved cycle overlap visualization',sec2time(toc)))

%% Save images

for cycle=1:num_cycles
    for channel=1:num_channels
        write_3d_tif(sprintf('%s/cy%02d_ch%01d.tif',reg_dir,cycle,channel),norm_stack(:,:,channel,cycle))
    end
end

disp(sprintf('%s: Wrote registered flat images to file',sec2time(toc)));

for channel=1:stain_num_channels
    write_3d_tif(sprintf('%s/st%01d.tif',reg_dir,channel),crop_stain_stack(:,:,channel))
end

disp(sprintf('%s: Wrote registered flat stain images to file',sec2time(toc)));

%% Nuclei segmentation

% watershed segmentation

draq5 = crop_stain_stack(:,:,2);
bw = (draq5>2500);
D = -bwdist(~bw);
D(~bw) = Inf;

L = watershed(imhmin(D,3));
L(~bw) = 0;
labels = label2rgb(L,'prism');

% filter cells by size
L = bwlabeln(L); % gets rid of isolated watershed pixels
cell_labels = regionprops(L, 'Area', 'BoundingBox');
cell_labels_filt = zeros(size(L));
cell_bounds = [];

for i=1:size(cell_labels,1)
    box = uint16(cell_labels(i).BoundingBox);
    if cell_labels(i).Area > 5000 %& box(1) > 2 & box(2) > 2 & box(1)+box(3) < crop_bounds(:,3)-1 & box(2)+box(4) < crop_bounds(:,4)-1
        cell_labels_filt(L==i)=i;
        cell_bounds = cat(1, cell_bounds, box);
    end
end

cell_labels_filt = bwlabeln(cell_labels_filt);
cells = regionprops(cell_labels_filt, 'BoundingBox','Image','Centroid','Area','PixelList');

% write segmentation to file

imwrite(uint16(cell_labels_filt),sprintf('%s/segmentation.tif',fov_dir));
dlmwrite(sprintf('%s/cell_bounds.txt',fov_dir), cell_bounds);
disp(sprintf('%s: Wrote 2D bounds for %d cells',sec2time(toc),size(cells,1)));

% display segmentation

figure; imshow(crop_stain_stack(:,:,2),[0 10000]); hold on;

[B L] = bwboundaries(cell_labels_filt,8,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'Color',colors(k,:), 'LineWidth', 2); hold on;
    text(cells(k).Centroid(1)-25,cells(k).Centroid(2),sprintf('%03d',k)); hold on;
end

fig_dir = sprintf('%s/figures/2d_segmentation',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(gcf,sprintf('%s/2d_segmentation.png',fov_dir));
saveas(gcf,sprintf('%s/fov%03d.png',fig_dir,fov));

%% Load intermediate data
%{
crop_bounds = dlmread(sprintf('%s/fov%03d/crop_bounds.txt',exp_dir,fov));

for cycle=1:num_cycles
    for channel=1:num_channels
        norm_stack(:,:,channel,cycle) = read_3d_tif(sprintf('%s/cy%02d_ch%01d.tif',reg_dir,cycle,channel),crop_bounds(4),crop_bounds(3),1);
    end
end

disp(sprintf('%s: Loaded registered flat images',sec2time(toc)));

for channel=1:stain_num_channels
    crop_stain_stack(:,:,channel) = read_3d_tif(sprintf('%s/st%01d.tif',reg_dir,channel),crop_bounds(4),crop_bounds(3),1);
end

disp(sprintf('%s: Loaded registered flat stain images',sec2time(toc)));

cell_labels_filt = imread(sprintf('%s/segmentation.tif',fov_dir));
cells = regionprops(cell_labels_filt,'BoundingBox','Image','Centroid','PixelList');

disp(sprintf('%s: Loaded cell segmentation',sec2time(toc)));

%}
%% Peak calling

all_peaks = {};
prc_thresh = 99.5;

for channel=2:5

    image = norm_stack(:,:,channel,1);
    peaks = FastPeakFind(image,prctile(image(:),prc_thresh));
    peaks = reshape(peaks,[],length(peaks)/2)'; peaks = [peaks(:,2) peaks(:,1)];
    all_peaks{channel} = peaks;

    maxima = image(find(sub2ind(size(image),peaks(:,1),peaks(:,2))));
    %figure; histogram(image,1:1:10000); set(gca,'YScale','log') 
    
    %figure; imshow(norm_stack(:,:,channel,1),[]); hold on;
    %plot(peaks(:,2),peaks(:,1),'ro','MarkerSize',5)
    
    disp(sprintf('%s: Found %d 3D peaks in cycle %d, channel %d',sec2time(toc),length(peaks),cycle,channel));
    
end

%% Combine peaks

combined_peaks = [];

for channel=2:5
    combined_peaks = cat(1,combined_peaks,all_peaks{channel});
end

combined_peaks = unique(combined_peaks,'rows');

figure; imshowpair(cap_image(max(norm_stack(:,:,2:5,1),[],3),10000),cap_image(flat_stain_stack(:,:,2),10000)); hold on;
plot(combined_peaks(:,2),combined_peaks(:,1),'ro','MarkerSize',5)

disp(sprintf('%s: Found %d total peaks across all channels in cycle 1',sec2time(toc),length(combined_peaks)));

%% Peak filtering

% filter peaks called on bright junk

junk_thresh = 10000;
junk_min_size = 100;
junk_image = crop_stain_stack(:,:,draq5_channel) > junk_thresh;
junk_labels = bwlabeln(junk_image);

junk_regions = regionprops(junk_labels, 'Area', 'BoundingBox');
junk_labels_filt = zeros(size(junk_labels));

for i=1:size(junk_regions,1)
    if junk_regions(i).Area > junk_min_size;
        junk_labels_filt(junk_labels==i)=i;
    end
end

junk_labels_filt = bwlabeln(junk_labels_filt);
junk_regions_filt = regionprops(junk_labels_filt,'PixelList');
junk_pixels = cat(1,junk_regions_filt.PixelList,zeros(0,2));
junk_pixels = [junk_pixels(:,2) junk_pixels(:,1)];

filtered_peaks = combined_peaks(~ismember(combined_peaks,junk_pixels,'rows'),:);

%figure; imshow(max(crop_stain_stack(:,:,2)>10000,[],3),[]); hold on;
%plot(filtered_peaks(:,2),filtered_peaks(:,1),'ro','MarkerSize',5)

% filter peaks within distance thresh from cell

dist_thresh = 100;

cell_pixels = cat(1,cells.PixelList);
cell_pixels = [cell_pixels(:,2) cell_pixels(:,1)];
pd_mat = pdist2(filtered_peaks,cell_pixels);
min_dist = min(pd_mat,[],2);

filtered_peaks = filtered_peaks(min_dist<dist_thresh,:);
figure; imshow(max(crop_stain_stack(:,:,2),[],3),[]); hold on;
plot(filtered_peaks(:,2),filtered_peaks(:,1),'ro','MarkerSize',5)

disp(sprintf('%s: Kept %d total peaks after filtering',sec2time(toc),length(filtered_peaks)));

%% Quantify spots

corr_thresh = 0.90;
spot_dim = [5 5];

fast_stack = norm_stack(:,:,2:5,:);
intensity_mat = zeros(length(filtered_peaks),num_channels-1,num_cycles+1);
    
for peak=1:length(filtered_peaks)
    intensity_mat(peak,:,:) = spot_caller_xy(filtered_peaks(peak,:),fast_stack,spot_dim,corr_thresh,[0 0]);
end

disp(sprintf('%s: Quantified all spots',sec2time(toc)));

%% Simple spot quant (PERFORMS WORSE)
%{
fast_stack = norm_stack(:,:,2:5,:);
intensity_mat = zeros(length(filtered_peaks),num_channels-1,num_cycles+1);
p = 2;

for i=1:length(filtered_peaks)
    peak = filtered_peaks(i,:);
    intensity_mat(i,1:2,1) = peak;
    intensity_mat(i,:,2:num_cycles+1) = reshape(sum(sum(fast_stack(peak(1)+[-p:p],peak(2)+[-p:p],:,:),1),2),1,num_channels-1,num_cycles);
end
%}

%% Sum of squares normalization

norm_intensity_mat = [];
for i=1:size(intensity_mat,1)
    norm_intensity_mat(i,:,1) = intensity_mat(i,:,1);
    norm_intensity_mat(i,:,2:num_cycles+1) = normc(squeeze(intensity_mat(i,:,2:num_cycles+1)));
end

disp(sprintf('%s: Performed sum of squares normalization',sec2time(toc)));

[purity_sort sort_order] = sort(prod(max((norm_intensity_mat(:,:,2:num_cycles+1).^2),[],2),3),'descend');
sort_intensity_mat = norm_intensity_mat(sort_order,:,:);
sort_peaks = filtered_peaks(sort_order,:);

disp(sprintf('%s: Sorted intensity matrix by spot purity',sec2time(toc)));

% Duplicate removal

[max_val seqs] = max(sort_intensity_mat(:,:,2:num_cycles+1),[],2);
seqs = seqs(:,:,1)*1000+seqs(:,:,2)*100+seqs(:,:,3)*10+seqs(:,:,4)*1;

pd = pdist2(sort_peaks,sort_peaks);
pd(pd==0) = 1000;

dup_thresh = 3;
[rows columns] = ind2sub(size(pd),find(pd<dup_thresh));
potential_dups = max(rows,columns);
dup_index = potential_dups(seqs(rows) == seqs(columns));

sort_intensity_mat(dup_index,:,:) = [];
sort_peaks(dup_index,:) = [];
purity_sort(dup_index) = [];
seqs(dup_index) = [];
seqs = seqs-1111;
%figure; histogram(purity_sort);

disp(sprintf('%s: Removed %d duplicate spots',sec2time(toc),size(unique(dup_index),1)));

%% Set variables based on exp

if exp_id == "italy_solid"
    probeset_1 = [0231 3300 2013 1122 1203 0110 3212 0323];
    labels_1 = ["ACTBref-G","TUBA1Bref-G","HSP90B1ref+C","L1TD1ref+C","ACTG1ref-C","PARP1ref-T","HSPD1ref-A","GSTP1ref+T"];
    probeset_2 = [1030 0002 3021 2332 1311 2101 2220 3133];
    labels_2 = ["ACTBalt-A","TUBA1Balt-A","HSP90B1alt+T","L1TD1alt+T","ACTG1alt-T","PARP1alt-C","HSPD1alt-G","GSTP1alt+C"];

    if fov < 6%4
        probeset = probeset_1;
        labels = labels_1;
    elseif fov < 11%7
        probeset = probeset_2;
        labels = labels_2;
    else
        probeset = cat(2,probeset_1,probeset_2);
        labels = cat(2,labels_1,labels_2);
    end
elseif exp_id == "italy_sbs"
    
    probeset_1 = [2301 0022 3210 1133 1320 2112 0313 2030];
    labels_1 = ["ACTBref-G","TUBA1Bref-G","HSP90B1ref+C","L1TD1ref+C","ACTG1ref-C","PARP1ref-T","HSPD1ref-A","GSTP1ref+T"];
    probeset_2 = [1202 2223 0231 3003 1011 3121 3332 0100];
    labels_2 = ["ACTBalt-A","TUBA1Balt-A","HSP90B1alt+T","L1TD1alt+T","ACTG1alt-T","PARP1alt-C","HSPD1alt-G","GSTP1alt+C"];

    if fov < 4
        probeset = probeset_1;
        labels = labels_1;
    elseif fov < 7
        probeset = probeset_2;
        labels = labels_2;
    else
        probeset = cat(2,probeset_1,probeset_2);
        labels = cat(2,labels_1,labels_2);
    end
    
end

%%

init_thresh = 0;

filter = sum(max(sort_intensity_mat(:,:,2:4).^2,[],2)>init_thresh,3)==3;
pos = squeeze(sort_intensity_mat(filter,1:2,1));
seqs = seqs(filter);


[C,ia,ic] = unique(seqs);
a_counts = accumarray(ic,1);
value_counts = [C,a_counts];
on_target = ismember(C,probeset);
    
figure; 
bar([1:size(on_target,1)],a_counts.*on_target,0.25,'r'); hold on;
bar([1:size(on_target,1)],a_counts.*~on_target,0.25,'b'); hold on;
ax = gca; ax.XTick = [1:length(C)]; 
set(gca,'xticklabel',pad(string(C),4,'left','0'),'fontsize',6); xtickangle(90)
title(sprintf('FOV %d, %02d%% on target before filtering',fov,round(sum(on_target.*a_counts)/sum(a_counts)*100)));

fig_dir = sprintf('%s/figures/pre_filter',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
saveas(gcf,sprintf('%s/pre_filter.png',fov_dir));
saveas(gcf,sprintf('%s/fov%03d.png',fig_dir,fov));

%fov_purity = min(max(sort_intensity_mat(:,:,2:num_cycles+1).^2,[],2),[],3);
fov_purity = log10(prod(max(sort_intensity_mat(:,:,2:num_cycles+1).^2,[],2),3));
fov_plot = length(seqs);
fov_color = ismember(seqs,probeset);

for i=1:length(C)
    fov_purity(seqs==C(i));
    fov_plot(seqs==C(i)) = i;
end
figure; scatter(fov_plot+((rand(length(fov_plot),1)-0.5)/5)',fov_purity,15,fov_color); hold on;
ax = gca; ax.XTick = [1:length(C)]; 
set(gca,'xticklabel',pad(string(C),4,'left','0'),'fontsize',6); xtickangle(90)
colormap(jet)

fig_dir = sprintf('%s/figures/purity_bar',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
saveas(gcf,sprintf('%s/purity_bar.png',fov_dir));
saveas(gcf,sprintf('%s/fov%03d.png',fig_dir,fov));

% filter

init_thresh = -0.5;
pass = fov_purity>init_thresh;
        
[C,ia,ic] = unique(seqs(pass));
a_counts = accumarray(ic,1);
value_counts = [C,a_counts];
on_target = ismember(C,probeset);
        
figure; 
bar([1:size(on_target,1)],a_counts.*on_target,0.25,'r'); hold on;
bar([1:size(on_target,1)],a_counts.*~on_target,0.25,'b'); hold on;
ax = gca; ax.XTick = [1:length(C)]; 
set(gca,'xticklabel',pad(string(C),4,'left','0'),'fontsize',6); xtickangle(90)
title(sprintf('FOV %d, %02d%% on target after filtering',fov,round(sum(on_target.*a_counts)/sum(a_counts)*100)));

fig_dir = sprintf('%s/figures/post_filter',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
saveas(gcf,sprintf('%s/post_filter.png',fov_dir));
saveas(gcf,sprintf('%s/fov%03d.png',fig_dir,fov));


%% Show filtered peaks

init_thresh = -1;
filtered_pos = pos(fov_purity>init_thresh & ismember(seqs,probeset),:);
filtered_seq = seqs(fov_purity>init_thresh & ismember(seqs,probeset));
[C,ic,ia] = unique(filtered_seq);
h = [];

figure; imshow(crop_stain_stack(:,:,2),[0 10000]); hold on;

[B L] = bwboundaries(cell_labels_filt,8,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'Color',colors(k,:), 'LineWidth', 2); hold on;
    text(cells(k).Centroid(1)-25,cells(k).Centroid(2),sprintf('%03d',k)); hold on;
end

for i=1:length(C)
    h(i) = scatter(filtered_pos(ia==i,2),filtered_pos(ia==i,1),50,colors(i,:),'filled'); hold on;
end

legend(h,pad(string(C),4,'left','0'))

%% Assign peaks to cells

% Create spot-cell distance matrix

sc_dist = zeros(length(cells),length(filtered_pos));

for i=1:length(cells)
    
    cell_pixels = [cells(i).PixelList(:,2) cells(i).PixelList(:,1)];
    dist_mat = pdist2(cell_pixels,filtered_pos);
    sc_dist(i,:) = min(dist_mat,[],1);
    
end

%% Show spots in cells

in_cell = (min(sc_dist,[],1)==0)';
%{
figure; imshow(crop_stain_stack(:,:,2),[0 10000]); hold on;

[B L] = bwboundaries(cell_labels_filt,8,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'Color',colors(k,:), 'LineWidth', 2); hold on;
    text(cells(k).Centroid(1)-25,cells(k).Centroid(2),sprintf('%03d',k)); hold on;
end

for i=1:length(C)
    h(i) = scatter(filtered_pos(ia==i&in_cell,2),filtered_pos(ia==i&in_cell,1),50,colors(i,:),'filled'); hold on;
end
%}
%% Infer cell assignments

assignment_mat = zeros(length(filtered_seq),5);
dist_factor = 2;

for i=1:length(filtered_seq)
   
    cell = -1;
    dist = -1;
    
    if min(sc_dist(:,i)) == 0
        cell = find(sc_dist(:,i)==0);
        dist = 0;
    else
        sort_dist = sort(sc_dist(:,i));
        if sort_dist(1) < sort_dist(2)/dist_factor
            cell = find(sc_dist(:,i)==sort_dist(1));
            dist = sort_dist(1);
        end
    end
    
    assignment_mat(i,1) = i;
    assignment_mat(i,2:3) = filtered_pos(i,:);
    assignment_mat(i,4) = filtered_seq(i,:);
    assignment_mat(i,5) = cell;
    assignment_mat(i,6) = dist;
    
end

figure; imshow(crop_stain_stack(:,:,2),[0 10000]); hold on;

[B L] = bwboundaries(cell_labels_filt,8,'noholes');
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'Color',colors(k,:), 'LineWidth', 2); hold on;
    text(cells(k).Centroid(1)-25,cells(k).Centroid(2),sprintf('%03d',k)); hold on;
end

probeset_count = [];
for i=1:length(probeset)
    h(i) = scatter(assignment_mat(assignment_mat(:,5)>0&assignment_mat(:,4)==probeset(i),3),assignment_mat(assignment_mat(:,5)>0&assignment_mat(:,4)==probeset(i),2),50,colors(i,:),'filled'); hold on;
    probeset_count(i) = sum(assignment_mat(:,4)==probeset(i));
end

%legend(h,pad(string(C),4,'left','0'))
legend(h,labels+" ("+string(probeset_count+")"))


high_conf = sum(assignment_mat(:,6)==0);
inferred =  sum(assignment_mat(:,6)>0);
unassigned =  sum(assignment_mat(:,6)==-1);
total = length(assignment_mat);
title(sprintf('%d high confidence cell assignments (%d%%), %d inferred assignments (%d%%), %d unassigned (%d%%)',high_conf,round(high_conf/total*100),inferred,round(inferred/total*100),unassigned,round(unassigned/total*100)))

fig_dir = sprintf('%s/figures/assigned_spots',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(gcf,sprintf('%s/assigned_spots.png',fov_dir));
saveas(gcf,sprintf('%s/fov%03d.png',fig_dir,fov));

%% Genotypes

if length(labels) == 8
    return
end

colors_special = [colors(4,:);colors(1,:);colors(3,:);colors(2,:)];

counts_table = zeros(length(cells),length(labels));
genotypes = zeros(length(cells),length(labels)/2);

for i=1:length(labels)/2
    
    gene = char(labels(i));
    gene = gene(1:end-5);
    
    ref_seq = probeset_1(i);
    alt_seq = probeset_2(i);
    
    for j=1:length(cells)
        cell_mat = assignment_mat(assignment_mat(:,5)==j,:);
        
        num_ref = sum(cell_mat(:,4) == ref_seq);
        num_alt = sum(cell_mat(:,4) == alt_seq);
        
        counts_table(j,(i-1)*2+1) = num_ref;
        counts_table(j,(i-1)*2+2) = num_alt;
        
        if num_ref == 0 & num_alt == 0
            genotypes(j,i) = -1;
        elseif num_alt/(num_ref+num_alt) < 0.25
            genotypes(j,i) = 0;
        elseif num_alt/(num_ref+num_alt) >= 0.25 & num_alt/(num_ref+num_alt) <= 0.75
            genotypes(j,i) = 1;
        else
            genotypes(j,i) = 2;
        end
    end
    
    figure; imshow(crop_stain_stack(:,:,2),[0 10000]); hold on;
    title(gene)
    
    [B L] = bwboundaries(cell_labels_filt,8,'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1),'Color',colors_special(genotypes(k,i)+2,:), 'LineWidth', 2); hold on;
        text(cells(k).Centroid(1)-25,cells(k).Centroid(2),sprintf('%03d',k)); hold on;
    end
   
    h = [];
    h(1) = scatter(assignment_mat(assignment_mat(:,4)==ref_seq&assignment_mat(:,5)>0,3),assignment_mat(assignment_mat(:,4)==ref_seq&assignment_mat(:,5)>0,2),50,colors(1,:),'filled'); hold on;
    h(2) = scatter(assignment_mat(assignment_mat(:,4)==alt_seq&assignment_mat(:,5)>0,3),assignment_mat(assignment_mat(:,4)==alt_seq&assignment_mat(:,5)>0,2),50,colors(2,:),'filled'); hold on;
    legend(h,["REF" "ALT"]);
    
    fig_dir = sprintf('%s/figures/gt_%s',exp_dir,gene);
    if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

    saveas(gcf,sprintf('%s/gt_%s.png',fov_dir,gene));
    saveas(gcf,sprintf('%s/fov%03d.png',fig_dir,fov));
    
    
end

dlmwrite(sprintf('%s/counts_table.txt',fov_dir),counts_table)
dlmwrite(sprintf('%s/genotypes.txt',fov_dir),genotypes)


%%

solid_counts = [];
for fov=11:15
    tmp_counts = dlmread(sprintf('%s/fov%03d/genotypes.txt',exp_dir,fov));
    solid_counts = cat(1,solid_counts,tmp_counts);
end
    

%%

pool_genotypes_file = '\\flynn-cifs/seq_Buenrostro/Users/Zack/italy_analysis/italy_genotypes.txt';
pool_genotypes = dlmread(pool_genotypes_file);

num_donors = size(pool_genotypes,1);
identity_scores = zeros(length(cells),num_donors);

for i=1:length(solid_counts)
    
    %counts_table(i,:)
    solid_counts(i,:);
    
    for j=1:num_donors
        
        pool_genotypes(j,:);
        
        score = 0;
        for k=1:length(probeset)/2
            
            %counts_table(i,(k-1)*2+1:(k-1)*2+2)
            gt_img = solid_counts(i,k);
            gt_seq = pool_genotypes(j,k);
            
            if gt_img == 0
                
                if gt_seq == 0
                    score = score + 1;
                elseif gt_seq == 1
                    score = score + 0;%0.5;
                else
                    score = score + 0;
                end
                
            elseif gt_img == 1
                
                if gt_seq == 0
                    score = score + 0;
                elseif gt_seq == 1
                    score = score + 1;
                else
                    score = score + 0;
                end
                
            elseif gt_img == 2
                
                if gt_seq == 0
                    score = score + 0;
                elseif gt_seq == 1
                    score = score + 0;%0.5;
                else
                    score = score + 1;
                end
                
            else
                score = score + 0.5;
            end
            
            %score
            
        end
        
        identity_scores(i,j) = score;
        
        
    end
end


%%
figure; imagesc(identity_scores)
colormap(redblue)

%%

figure;
cgo_all = clustergram(identity_scores,'Colormap',redbluecmap,'Standardize','Row')

%%
Y = pdist(pool_genotypes);
Z = linkage(Y)
figure; [H,T,outperm] = dendrogram(Z);

%%

gt_dist_mat = zeros(length(pool_genotypes));
good_probes = 1:8%[3 4 5 6 8];
%clustering_order = [16 3 9 17 8 5 13 12 10 11 15 14 7 2 1 4 6 18];
%pool_genotypes2 = pool_genotypes(clustering_order,:);

for i=1:length(pool_genotypes)
    for j=1:length(pool_genotypes)
        
        pool_genotypes(i,:);
        pool_genotypes(j,:);
        gt_dist_mat(i,j) = sum(abs(pool_genotypes(i,good_probes)-pool_genotypes(j,good_probes)));
        
        
    end
end

figure; imagesc(gt_dist_mat)
colormap(redblue)
caxis([0 12])
xlabel('Donor')
ylabel('Donor')
%figure;
%cgo_all = clustergram(gt_dist_mat,'Colormap',redbluecmap,'Standardize','Row')

%% Compare solid vs. SBS

solid_counts = zeros(1,16);
solid_num_cells = 0;

solid_dir = 'italy_solid';
for fov=11:15
    tmp_counts = dlmread(sprintf('%s/fov%03d/counts_table.txt',solid_dir,fov));
    solid_counts = solid_counts+sum(tmp_counts,1);
    solid_num_cells = solid_num_cells + size(tmp_counts,1);
end

sbs_counts = zeros(1,16);
sbs_num_cells = 0;

sbs_dir = 'italy_sbs';
for fov=7:9
    tmp_counts = dlmread(sprintf('%s/fov%03d/counts_table.txt',sbs_dir,fov));
    sbs_counts = sbs_counts+sum(tmp_counts,1);
    sbs_num_cells = sbs_num_cells + size(tmp_counts,1);
end

figure; scatter(solid_counts./solid_num_cells,sbs_counts./sbs_num_cells)
 axis([0 5 0 5]); hold on;

%%


%figure; imshowpair(max(norm_stack(:,:,3,1),[],3),max(norm_stack(:,:,4,1),[],3)); hold on;

%plot(sort_peaks(seqs==1122,2),sort_peaks(seqs==1122,1),'ro','MarkerSize',5)
%plot(sort_peaks(seqs==1222,2),sort_peaks(seqs==1222,1),'bo','MarkerSize',5)
%plot(sort_peaks(seqs==2122,2),sort_peaks(seqs==2122,1),'go','MarkerSize',5)



















return

%%

%% Quantify each cell

for cell=1:size(cells,1)
%for cell=2:2
    mito_cell(exp_id,fov,cell,cells,nonneg_stack,crop_stain_stack,otsu_thresh)
end

return

%%























%%








%%
pd = pdist2(combined_peaks,combined_peaks);


pd(pd==0) = 1000;
figure; histogram(min(pd),0:0.5:10)
matches = min(pd)<1.5;
matched_peaks = combined_peaks(matches,:);

figure; imshow(max(norm_stack(:,:,2:5,1),[],3),[0 5000]); hold on;
plot(matched_peaks(:,2),matched_peaks(:,1),'rx');
%%








%%
err_pos = pos(seq==333,:);
figure; imshowpair(max(norm_stack(:,:,5,1),[],3)>10000,max(norm_stack(:,:,3,1),[],3)>10000); hold on;
plot(err_pos(:,2),err_pos(:,1),'rx')
%% Visualize high confidence calls


pos = squeeze(sort_intensity_mat((purity_sort > 0.5),1:2,1));
[max_val seq] = max(sort_intensity_mat((purity_sort > 0.5),:,2:4),[],2);
seq = seq(:,:,1)*100+seq(:,:,2)*10+seq(:,:,3);
uniq = unique(seq);
colors = distinguishable_colors(size(uniq,1));

%for i=1:size(seq,1)
%    markers(i) = find(uniq==seq(i));
%end


%figure; imshow(labels,[]); hold on;
%scatter(pos(:,2),pos(:,1),10,colors(markers',:),'filled')
%plot(pos(:,2),pos(:,1),'rx')


%%

cells = regionprops(cell_labels_filt,'Centroid','Area','BoundingBox','PixelList');
matrix = [];
count = 0;

for i=1:size(cells,1)
    
    pixels = [cells(i).PixelList(:,2) cells(i).PixelList(:,1)];
    [inside ia ib] = intersect(pixels,pos,'rows');
    seq(ib);
    for j=1:size(ib)
        count = count + 1;
        matrix(count,:) = [i seq(ib(j)) pos(ib(j),:)];
    end
    
end

%%

calls_matrix = zeros(size(cells,1),size(colorspace,1));
counts_matrix = zeros(size(cells,1),size(colorspace,1)*2);

for cell=1:size(cells,1)

%figure('visible','off');
%p = tight_subplot(1,size(colorspace,1),[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for gene=1:size(colorspace,1)
  
    cell_matrix = matrix(matrix(:,1) == cell,:);
    
    num_ref = sum(cell_matrix(:,2) == colorspace{gene,2});
    num_alt = sum(cell_matrix(:,2) == colorspace{gene,3});
    
    counts_matrix(cell,(gene-1)*2+[1:2]) = [num_ref num_alt];

    %pos_ref = cell_matrix(cell_matrix(:,2) == colorspace{gene,2},:);
    %pos_ref(:,3) = pos_ref(:,3)-(box(2)-pad);
    %pos_ref(:,4) = pos_ref(:,4)-(box(1)-pad);
    
    %pos_alt = cell_matrix(cell_matrix(:,2) == colorspace{gene,3},:);
    %pos_alt(:,3) = pos_alt(:,3)-(box(2)-pad);
    %pos_alt(:,4) = pos_alt(:,4)-(box(1)-pad);
    
    if num_ref == 0 & num_alt == 0
        gt = -1;
    elseif num_alt/(num_ref+num_alt) < 0.33
        gt = 0;
    elseif num_alt/(num_ref+num_alt) >= 0.33 & num_alt/(num_ref+num_alt) <= 0.66
        gt = 1;
    else
        gt = 2;
    end
        
    calls_matrix(cell,gene) = gt;
    
    %box = round(cells(cell).BoundingBox);
    %xy_dim = max(box(3),box(4));
    %axes(p(gene)); imshow(reg_stain_stack(box(2)-pad:box(2)+xy_dim+pad,box(1)-pad:box(1)+xy_dim+pad,3),[]); hold on;
    
    %scatter(pos_ref(:,4),pos_ref(:,3),5,'r'); hold on;
    %scatter(pos_alt(:,4),pos_alt(:,3),5,'b'); hold on;
    %title(sprintf('%s (%d)\nR=%d, A = %d ',colorspace{gene,1}{1},gt,num_ref,num_alt))
    
end

%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 36 2.5];
%saveas(fig,sprintf('cell%03d.png',cell));

end

dlmwrite('fov037_counts.txt',counts_matrix);

%%

test2 = [];
test = calls_matrix + 2;
for i=1:size(test,1)
    test2(i,:) = reshape(full(ind2vec(test(i,:),4)),[],1);
end

%%


k = 11;
colors = distinguishable_colors(k+1);
colors(4,:) = [];
Y = tsne(test2);
T = kmeans(test2,k);
%[centers,U] = fcm(test,k);
figure; scatter(Y(:,1),Y(:,2),10,colors(T,:))
colormap(prism)

%%
Y = pdist(test2);
squareform(Y);
Z = linkage(Y,'average');
figure; dendrogram(Z,0);
T = cluster(Z,'cutoff',0.9);
colors = distinguishable_colors(max(T));

X = tsne(test2);
figure; scatter(X(:,1),X(:,2),10,colors(T,:))
colormap(prism)

%%
%%

figure; imshow(reg_stain_stack(:,:,3),[0 10000]); hold on;

[B L] = bwboundaries(cell_labels_filt,'noholes');

for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1),'Color',colors(T(k),:), 'LineWidth', 2); hold on;
end


%%
colorspace = readtable('probe_colorspace.txt');
colorspace{:,2} = floor((colorspace{:,2} + 1111)/10);
colorspace{:,3} = floor((colorspace{:,3} + 1111)/10);

for gene=1:size(colorspace,1)
%for gene=1:1
    
    colorspace{gene,1};
    
    num_ref = sum(seq == colorspace{gene,2});
    num_alt = sum(seq == colorspace{gene,3});
    
    pos_ref = pos(seq == colorspace{gene,2},:);
    pos_alt = pos(seq == colorspace{gene,3},:);
    
    figure; imshow(reg_stain_stack(:,:,3),[0 10000]); hold on;
    title(sprintf('%s, REF=%d, ALT=%d',colorspace{gene,1}{1},num_ref,num_alt))
    scatter(pos_ref(:,2),pos_ref(:,1),10,'r'); hold on;
    scatter(pos_alt(:,2),pos_alt(:,1),10,'b'); hold on;
    
    ref_cells = unique(matrix((matrix(:,2) == colorspace{gene,2}),1));
    alt_cells = unique(matrix((matrix(:,2) == colorspace{gene,3}),1));
    
    
    for i=1:size(cells,1)
        
        if any(ref_cells==i) & any(alt_cells==i)
            [B L] = bwboundaries(cell_labels_filt==i);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
            end
        elseif any(ref_cells==i)
            [B L] = bwboundaries(cell_labels_filt==i);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
            end
        elseif any(alt_cells==i)
            [B L] = bwboundaries(cell_labels_filt==i);
            for k = 1:length(B)
               boundary = B{k};
               plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2)
            end
        end
        
    end
    
    saveas(gcf,sprintf('ref_vs_alt2/%s.png',colorspace{gene,1}{1}));
    
end

%%
thresh = 0.5;
ts = sum(max(sort_intensity_mat(:,:,2:4).^2,[],2)>thresh,3)==3;
pos = squeeze(sort_intensity_mat(ts,1:2,1));
[max_val seq] = max(sort_intensity_mat(ts,:,2:4),[],2);
seq = seq(:,:,1)*100+seq(:,:,2)*10+seq(:,:,3);

a = length(seq)
b = sum(ismember(seq,all_barcodes,'rows'))
b/a

bad_spots = ~ismember(seq,all_barcodes,'rows');
pos_bad_spots = sort_intensity_mat(bad_spots,1:2,1);
figure; imshow(max(reg_stack(:,:,2:5,1),[],3),[0 1000]); hold on;
plot(pos_bad_spots(:,2),pos_bad_spots(:,1),'rx')

%%
%{
colorspace = readtable('probe_colorspace.txt');
colorspace{:,2} = floor((colorspace{:,2} + 1111)/10);
colorspace{:,3} = floor((colorspace{:,3} + 1111)/10);

all_barcodes = cat(1,colorspace{:,2},colorspace{:,3});
figure;
tp = 27437;
fp = 35930 - 27437;


for thresh=0.01:0.01:1

    ts = sum(max(sort_intensity_mat(:,:,2:4).^2,[],2)>thresh,3)==3;
    pos = squeeze(sort_intensity_mat(ts,1:2,1));
    [max_val seq] = max(sort_intensity_mat(ts,:,2:4),[],2);
    %pos = squeeze(sort_intensity_mat((purity_sort > thresh),1:2,1));
    %[max_val seq] = max(sort_intensity_mat((purity_sort > thresh),:,2:4),[],2);
    seq = seq(:,:,1)*100+seq(:,:,2)*10+seq(:,:,3);

    a = length(seq);
    b = sum(ismember(seq,all_barcodes,'rows'));
    scatter((a-b)/fp,b/tp); hold on;
end
%%

%%

%%





%%

cells = regionprops(cell_labels_filt,'PixelList','Centroid','BoundingBox');
matrix = [];
count = 0;

maf_matrix = zeros(size(cells,1),size(colorspace,1));

for i=1:size(cells,1)
%for i=50:50
    
    i;
    pixels = [cells(i).PixelList(:,2) cells(i).PixelList(:,1)];
    [inside ia ib] = intersect(pixels,pos,'rows');
    seqs = seq(ib);
    
    
    for gene=1:size(colorspace,1)
        num_ref = sum(seqs == colorspace{gene,2}) + 1;
        num_alt = sum(seqs == colorspace{gene,3}) + 1;
        maf_matrix(i,gene) = num_ref/(num_ref+num_alt);
    end
    
end








%}
