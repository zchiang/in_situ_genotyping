%% Load in data

seed = rng; % random seed

snp_by_donor_mat = dlmread('high_exp.012');
snp_by_donor_mat = snp_by_donor_mat(:,2:size(snp_by_donor_mat,2));
donors = readtable('high_exp.012.indv','FileType','text','delimiter',',','ReadVariableNames',false);
snps = readtable('high_exp.012.all_info','FileType','text','delimiter','\t');
snps.Properties.VariableNames = {'chr' 'pos' 'gene' 'exon' 'strand','d0_exp','d28_exp','ref','alt','ref_seq','alt_seq',...
    'probes5','probes3','sel','region5','sel_probe','region3','has_repeat','gc5','gc3'};

snp_by_donor_mat(:,string(snps.chr) == "chrX") = [];
snps(string(snps.chr) == "chrX",:) = [];

snp_by_donor_mat(:,string(snps.sel) == '0') = [];
snps(string(snps.sel) == '0',:) = [];

snp_by_donor_mat(:,string(snps.alt_seq) == "NA") = [];
snps(string(snps.alt_seq) == "NA",:) = [];

snp_by_donor_mat(:,snps.has_repeat == 1) = [];
snps(snps.has_repeat == 1,:) = [];

snp_by_donor_mat(:,snps.gc5 < 6 | snps.gc5 > 12 | snps.gc3 < 6 | snps.gc3 > 12) = [];
snps(snps.gc5 < 6 | snps.gc5 > 12 | snps.gc3 < 6 | snps.gc3 > 12,:) = [];

snps.sel =  str2num(char(string(snps.sel)));

num_probes = 16;
num_donors = size(donors,1);
num_snps = size(snps,1);

sel_snps = [];
sel_snp_by_donor_mat = [];
sel_dist_mat = diag(ones(num_donors,1)*10);

%%

for probe=1:num_probes
    
disp(sprintf("Selecting probe %d",probe));

% Create donor distance mat

abs_donor_mat = zeros(num_donors,num_donors);

for i=1:num_donors
    for j=1:num_donors
        missing_info = (snp_by_donor_mat(i,:) ~= -1 & snp_by_donor_mat(j,:) ~= -1);
        abs_donor_mat(i,j) = sum(abs(snp_by_donor_mat(i,~missing_info) - snp_by_donor_mat(j,~missing_info)));
    end
end

%abs_donor_mat = abs_donor_mat .* triu(ones(size(abs_donor_mat)));
abs_donor_mat(abs_donor_mat == 0) = 1000000;

% Create priority mat

min_dist = min(sel_dist_mat(:));
min_mat = sel_dist_mat == min_dist;

priority_mat = abs_donor_mat .* min_mat;

disp(sprintf("%d pairs of donors are %d alleles apart",sum(min_mat(:)),min_dist))

% Create distance mat for each SNP

snp_dist_mat = zeros(num_donors,num_donors,num_snps);

for i=1:num_snps  
    snp_dist_mat(:,:,i) = abs(repmat(snp_by_donor_mat(:,i),1,num_donors)-repmat(snp_by_donor_mat(:,i)',num_donors,1));
    snp_dist_mat(snp_by_donor_mat(:,i)==-1,:,i) = 0;
    snp_dist_mat(:,snp_by_donor_mat(:,i)==-1,i) = 0;
end

% Find best SNP

snp_var = nansum(nansum((snp_dist_mat.*min_mat)./priority_mat,1),2);
[max_var max_snp_pos] = max(snp_var);

% Add best SNP

max_gene = snps.gene{max_snp_pos};
sel_snps = [sel_snps; snps(max_snp_pos,:)];
sel_snp_by_donor_mat = [sel_snp_by_donor_mat; snp_by_donor_mat(:,max_snp_pos)'];
sel_dist_mat = sel_dist_mat + snp_dist_mat(:,:,max_snp_pos);

disp(sprintf("Selected SNP on gene %s at %s:%d",max_gene,snps.chr{max_snp_pos},snps.pos(max_snp_pos)))

% Remove all other SNPs in gene from matrix

snp_by_donor_mat(:,string(snps.gene) == max_gene) = [];
snps(string(snps.gene) == max_gene,:) = [];
num_snps = size(snps,1);


end
%% Choose donors based on selected distance mat

num_sel_donors = 10;

% choose random first donor

rng(1); % for pool_10
%rng(0); %for pool_20

sel_donors = [randi(num_donors)];
tmp_dist_mat = sel_dist_mat .* ~diag(ones(size(sel_dist_mat,1),1));

for i=2:num_sel_donors
    [max_dist furthest_donor] = max(min(tmp_dist_mat(sel_donors,:),[],1),[],2);
    sel_donors = [sel_donors furthest_donor];
end

pool_snps = sel_snps;
pool_donors = donors(sel_donors,:);
pool_dist_mat = tmp_dist_mat(sel_donors,sel_donors);
pool_snp_by_donor_mat = sel_snp_by_donor_mat(:,sel_donors);

writetable(pool_snps,sprintf('pool_%d/pool_%d_snps.txt',num_sel_donors,num_sel_donors));
writetable(pool_donors,sprintf('pool_%d/pool_%d_donors.txt',num_sel_donors,num_sel_donors));
dlmwrite(sprintf('pool_%d/pool_%d_dist_mat.txt',num_sel_donors,num_sel_donors),pool_dist_mat);
dlmwrite(sprintf('pool_%d/pool_%d_snp_by_donor_mat.txt',num_sel_donors,num_sel_donors),pool_snp_by_donor_mat);









