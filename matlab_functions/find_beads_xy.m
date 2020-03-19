function[centroids_new centroid_locs] = find_beads_xy(stack,bead_thresh,bead_size)

    xlen = size(stack,1);
    ylen = size(stack,2);
    
    conn = bwconncomp(stack>bead_thresh);
    regions = regionprops(conn,'Centroid','BoundingBox');
    centroids = cat(1, regions.Centroid);

    centroids_new = [];
    centroid_locs=zeros(xlen,ylen);

    for i = 1:length(centroids(:,1))
        
        x = round(centroids(i,2));
        y = round(centroids(i,1));
        
        xlim = [max(1,x-bead_size(1)), min(xlen,x+bead_size(1))];
        ylim = [max(1,y-bead_size(2)), min(ylen,y+bead_size(2))];
        
        area = stack(xlim(1):xlim(2),ylim(1):ylim(2));
        
        if x>bead_size(1) & x<xlen-bead_size(1) & y>bead_size(2) & y<ylen-bead_size(2)% & median(area(:)) < bead_thresh(2)
        
            centroids_new = [centroids_new; centroids(i,:)];
    	    centroid_locs(xlim(1):xlim(2),ylim(1):ylim(2)) = 1;

        end
    end
    
end
