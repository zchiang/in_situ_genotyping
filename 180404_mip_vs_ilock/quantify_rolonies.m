function[] = quantify_rolonies(nd2)
%% Start fresh

clearvars -except nd2
nd2 = char(nd2);
tic
addpath('D:\MicroscopyShare\zack\in_situ_toolkit\matlab_functions');
nd2_dir = 'D:\MicroscopyShare\TT\180404_mipvsiLock\iLock'; 
%nd2 = 'iLock_5a.nd2'; 

disp(sprintf('%s: Starting processing of %s',sec2time(toc),nd2));

%% Set parameters

exp_id = nd2(1:end-4);

%exp_dir = sprintf('%s',exp_id);
%if ~exist(exp_dir, 'dir') mkdir(exp_dir), end

%disp(sprintf('%s: Output files will be saved to %s',sec2time(toc),exp_dir));

%%

reader = bfGetReader(sprintf('%s/%s', nd2_dir, nd2));
fov_xlen = reader.getSizeX;
fov_ylen = reader.getSizeY;
fov_zlen = reader.getSizeZ;
num_channels = reader.getSizeC;
num_cycles = reader.getSizeT;
stack = zeros(fov_xlen,fov_ylen,fov_zlen,num_channels);

rolony_channel = 1;
wga_channel = 2;

disp(sprintf('%s: X = %d, Y = %d, Z = %d, C = %d, T = %d',sec2time(toc),fov_xlen,fov_ylen,fov_zlen,num_channels,num_cycles));

%% Load images

t = 1;
series = 1;
for channel=1:num_channels
    for z=1:fov_zlen
        stack(:,:,z,channel) = readPlane(reader,series,z,channel,t);
    end
end

flat_stack = squeeze(max(stack,[],3));

disp(sprintf('%s: Loaded all images',sec2time(toc)));

%% Peak calling

all_peaks = {};
prc_thresh = 95;
rolony_channel = 1;

image = flat_stack(:,:,rolony_channel);

[counts,x] = imhist(image,max(image(:)));
[T,EM] = otsuthresh(counts);
otsu_thresh = T*65536;

peaks = FastPeakFind(image,prctile(image(:),prc_thresh));
peaks = reshape(peaks,[],length(peaks)/2)'; peaks = [peaks(:,2) peaks(:,1)];
maxima = image(find(sub2ind(size(image),peaks(:,1),peaks(:,2))));
    
disp(sprintf('%s: Identified %d 2D peaks',sec2time(toc),size(maxima,1)));

%% Visualize peaks

figure; imshow(flat_stack(:,:,rolony_channel),[0 prctile(image(:),prc_thresh)]); hold on;
%figure; imshowpair(flat_stack(:,:,rolony_channel),flat_stack(:,:,wga_channel)); hold on;
plot(peaks(:,2),peaks(:,1),'rx')
title(sprintf('%d rolonies identified',size(maxima,1)))

saveas(gcf,sprintf('%s.png',exp_id));

%%

statistics = fopen('rolony_counts.txt','a');
fprintf(statistics,'%s\t%d\n',exp_id,size(maxima,1));
%fprintf(statistics,'%s\t%d\t%d\n',exp_id,size(nuclei_3d_bounds,1),size(filtered_peaks,1));

disp(sprintf('%s: Wrote statistics for %s to file',sec2time(toc),exp_id));

end