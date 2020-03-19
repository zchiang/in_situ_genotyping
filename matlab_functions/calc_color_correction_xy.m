function[stackBr cc_offsets] = calc_color_correction_xy(stackA, stackB, centroids, centroid_locs, bead_size)

    metric = registration.metric.MeanSquares;
    optimizer = registration.optimizer.RegularStepGradientDescent;
    optimizer.MaximumIterations = 50;

    cc_offsets = zeros(2,1);
    cc_transform = zeros(3,3);

    stackBr = stackB;
    max_corr = corr(stackA(centroid_locs>0),stackB(centroid_locs>0));
    init_corr = max_corr;
    %disp(sprintf('%s: Intial correlation is %.05f',sec2time(toc),max_corr));

    %for i=1:size(centroids,1)
    for i=1:min(size(centroids,1),10)
               
        x = round(centroids(i,2)); y = round(centroids(i,1));
        xlen = length(stackA(:,1,1)); ylen = length(stackA(1,:,1));
        xbox = bead_size(2); ybox = bead_size(1);   
        xlim = [max(1,x-xbox), min(xlen,x+xbox)]; ylim = [max(1,y-ybox), min(ylen,y+ybox)];
        
        boxA = stackA(xlim(1):xlim(2),ylim(1):ylim(2));
        boxB = stackB(xlim(1):xlim(2),ylim(1):ylim(2));
        
        if min(boxB(:)) == max(boxB(:))
            continue
        end
        
        imgA = max(boxA,[],3);
        imgB = max(boxB,[],3);
        
        c = normxcorr2_general(imgB, imgA, round(size(imgA,1)*size(imgA,2)*.8));
        [max_c, imax] = max(abs(c(:)));
        [ypeak, xpeak] = ind2sub(size(c),imax(1));
        corr_offset = [(xpeak-size(imgB,2)) (ypeak-size(imgB,1))];
        xoffset = corr_offset(1);
        yoffset = corr_offset(2);    
                
        %disp(sprintf('%s: Offset is %d (x) by %d (y)',sec2time(toc),xoffset,yoffset));

        xbegin_A = max(1+xoffset, 1);
        xend_A   = min(size(stackA,2), size(stackA,2)+xoffset);
        ybegin_A = max(1+yoffset, 1);
        yend_A   = min(size(stackA,1), size(stackA,1)+yoffset);

        xbegin_B = max(1-xoffset, 1);
        xend_B   = min(size(stackA,2), size(stackA,2)-xoffset);
        ybegin_B = max(1-yoffset, 1);
        yend_B   = min(size(stackA,1), size(stackA,1)-yoffset);
     
        box_xbegin_A = max(1+xoffset, 1);
        box_xend_A   = min(size(boxA,2), size(boxA,2)+xoffset);
        box_ybegin_A = max(1+yoffset, 1);
        box_yend_A   = min(size(boxA,1), size(boxA,1)+yoffset);

        box_xbegin_B = max(1-xoffset, 1);
        box_xend_B   = min(size(boxA,2), size(boxA,2)-xoffset);
        box_ybegin_B = max(1-yoffset, 1);
        box_yend_B   = min(size(boxA,1), size(boxA,1)-yoffset);
        
        stackBc1 = zeros(xlen,ylen);
        stackBc1(ybegin_A:yend_A,xbegin_A:xend_A) = stackB(ybegin_B:yend_B,xbegin_B:xend_B);
        stackBc1(stackBc1 == 0) = mode(stackB(:));

        boxBc = zeros(xbox,ybox);
        boxBc(box_ybegin_A:box_yend_A,box_xbegin_A:box_xend_A) = boxB(box_ybegin_B:box_yend_B,box_xbegin_B:box_xend_B);
        boxBc(boxBc == 0) = mode(boxB(:));
        
        offset_corr = corr(stackA(centroid_locs>0),stackBc1(centroid_locs>0),'type','Spearman');
        if offset_corr > max_corr
            max_corr = offset_corr;
            stackBr = stackBc1;
            cc_offsets = [xoffset yoffset];
            cc_transform = zeros(3,3);
        end

        %disp(sprintf('%s: Offset correlation is %.05f',sec2time(toc),offset_corr));

        
        [imreg_new, ~, tform_new] = imregister2(boxBc, boxA, 'translation', optimizer, metric, 'PyramidLevels', 1);
        reg_img = imwarp(stackBc1, tform_new, 'outputView', imref2d(size(stackBc1)));;
        reg_corr = corr(stackA(centroid_locs>0),reg_img(centroid_locs>0),'type','Spearman');
        %figure; scatter(stackA(centroid_locs>0),reg_img(centroid_locs>0));
        
        %if reg_corr > max_corr
        if reg_corr > max_corr
            max_corr = reg_corr;
            stackBr =  reg_img;
            max_tform = tform_new;
            cc_offsets = [xoffset+max_tform.T(3,1) yoffset+max_tform.T(3,2)];
            cc_transform = max_tform.T;
        end
        
        %disp(sprintf('%s: Registration correlation is %.05f',sec2time(toc),reg_corr));
        
    end
    

   
    %disp(sprintf('%s: Final correlation is %.05f',sec2time(toc),max_corr));
    %disp(sprintf('%s: Improved correlation from %.03f to %.03f',sec2time(toc),init_corr,max_corr));
    
end
