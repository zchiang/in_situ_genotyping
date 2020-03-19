function[] = rna_exp(exp_id,nd2_dir,exp_design_file)
% Author: Zack Chiang, Harvard University & Broad Institute
% Analysis of in situ ATAC-seq experiment

%% Set manual parameters, clear environment, add function paths, create directories

% ------ Set manual parameters (comment out if running as function) ------ %

exp_id = 'italy_sbs';
nd2_dir = 'X:\Confocal Data June2018-\Tongtong\180802_ITALYpool\SBS\';
exp_design_file = 'italy_sbs_exp_design.txt';

% ------------------------------------------------------------------------ %

tic
clearvars -except exp_id nd2_dir exp_design_file
addpath(genpath('X:\Confocal Data June2018-\Zack\in_situ\rna\matlab_functions\'));
disp(sprintf('%s: Started analysis of experiment %s',sec2time(toc),exp_id));

% Create directories

home_dir = 'X:\Confocal Data June2018-\Zack\in_situ\rna';
cd(home_dir)

exp_dir = sprintf('%s\\%s',home_dir,exp_id);
if ~exist(exp_dir, 'dir') mkdir(exp_dir), end

fig_dir = sprintf('%s\\figures',exp_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

disp(sprintf('%s: Output will be saved to %s',sec2time(toc),exp_dir));

%% Load .nd2 reader, set experimental parameters

exp_design = readtable(sprintf('%s/exp_design/%s',home_dir,exp_design_file));
nd2_reader = bfGetReader(sprintf('%s/%s',nd2_dir,exp_design.nd2{1}));

fov_xlen = nd2_reader.getSizeX;
fov_ylen = nd2_reader.getSizeY;
fov_zlen = nd2_reader.getSizeZ;

num_channels = nd2_reader.getSizeC;
num_cycles = length(exp_design.cycle)-1;
num_fov = nd2_reader.getSeriesCount;

stain_num_channels = 2;
wga_channel = 1;
draq5_channel = 2;

disp(sprintf('%s: X = %d, Y = %d, Z = %d, C = %d, T = %d, FOVs = %d',sec2time(toc),fov_xlen,fov_ylen,fov_zlen,num_channels,num_cycles,num_fov));

%% Save all parameters

clear nd2_reader
filename = sprintf('%s\\exp_parameters.mat',exp_dir); save(filename)
disp(sprintf('%s: Saved all parameters to %s',sec2time(toc),filename));

%% Per FOV

%for fov=1:num_fov
for fov=10:12%num_fov
    rna_fov(exp_id,fov);
end
