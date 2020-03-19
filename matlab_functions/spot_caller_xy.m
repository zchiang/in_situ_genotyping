function[mat] = spot_caller_xy(peak, all_stacks,spot_dim,corr_thresh,dims)

    num_channels = size(all_stacks,3);
    num_cycles = size(all_stacks,4);
    
    peak_pixels = reshape(squeeze(all_stacks(peak(1),peak(2),:,:)),[],1);
    spot_area = zeros(spot_dim(1)*2+1,spot_dim(2)*2+1,num_channels,num_cycles);

    x_bounds = [max(1, peak(1)-spot_dim(1)) min(size(all_stacks,1), peak(1)+spot_dim(1))];
    y_bounds = [max(1, peak(2)-spot_dim(2)) min(size(all_stacks,2), peak(2)+spot_dim(2))];

    x_offsets = [1-(peak(1)-spot_dim(1)-x_bounds(1)) size(spot_area,1)-(peak(1)+spot_dim(1)-x_bounds(2))];
    y_offsets = [1-(peak(2)-spot_dim(2)-y_bounds(1)) size(spot_area,2)-(peak(2)+spot_dim(2)-y_bounds(2))];

    spot_area(x_offsets(1):x_offsets(2),y_offsets(1):y_offsets(2),:,:) = all_stacks(x_bounds(1):x_bounds(2), y_bounds(1):y_bounds(2),:,:);
    area_pixels = size(spot_area,1)*size(spot_area,2);
    corr_mat = corrcoef(reshape(spot_area,area_pixels,num_channels*num_cycles)');
    correlations = reshape(corr_mat((area_pixels-1)/2+1,:),size(spot_area,1),size(spot_area,2));
    spot = correlations;

    regions = regionprops(correlations > corr_thresh, 'Centroid');
    pixel_count = sum(correlations(:) > corr_thresh);

    mat = zeros(num_channels, num_cycles+1);
    %max(all_stacks(peak(1),peak(2),:,1))
    mat(1:3,1) = [peak(1)+dims(2)-1 peak(2)+dims(1)-1 pixel_count];
    mat(4,1) = max(all_stacks(peak(1),peak(2),:,1));


    %spot_area.*(repmat(correlations,1,1,4,num_cycles) > corr_thresh)
        %mat(:,2:num_cycles+1) = squeeze(sum(sum(sum((all_stacks.*(repmat(vol,1,1,1,4,num_cycles) > 0)),1),2),3))
    mat(:,2:num_cycles+1) = squeeze(sum(sum((spot_area.*(repmat(correlations,1,1,4,num_cycles) > corr_thresh)),1),2));


    %else
        %mat = [];
        %spot = [];
        %status = 0;
    %end

    %disp(sprintf('Spot quantified')); toc

end