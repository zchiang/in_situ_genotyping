function[stackBr cy_offset] = fov_offset_xy(stackA, stackB, wga_channel)

    xlen = length(stackA(:,1,1));
    ylen = length(stackA(1,:,1));
    num_channels = length(stackA(1,1,:));
    
    %imgA = max(stackA,[],3);
    %imgB = max(stackB,[],3);
    imgA = max(stackA(:,:,wga_channel),[],3);
    imgB = max(stackB(:,:,wga_channel),[],3);

    c = normxcorr2_general(imgB, imgA);
    %figure; contourf(c)
    if xlen>1000
        bounds = 500;
    else
        bounds = 0;
    end
    mask = zeros(size(c));
    mask(xlen-bounds:xlen+bounds,ylen-bounds:ylen+bounds) = 1;
    c = c.*mask;
   
    [max_c, imax] = max(abs(c(:)));
    [ypeak, xpeak] = ind2sub(size(c),imax(1));
    cy_offset = [(xpeak-size(imgB,2)) (ypeak-size(imgB,1))];

    xoffset = cy_offset(1);
    yoffset = cy_offset(2);

    xbegin_A = max(1+xoffset, 1);
    xend_A   = min(size(imgA,2), size(imgA,2)+xoffset);
    ybegin_A = max(1+yoffset, 1);
    yend_A   = min(size(imgA,1), size(imgA,1)+yoffset);

    xbegin_B = max(1-xoffset, 1);
    xend_B   = min(size(imgA,2), size(imgA,2)-xoffset);
    ybegin_B = max(1-yoffset, 1);
    yend_B   = min(size(imgA,1), size(imgA,1)-yoffset);

    stackBr = zeros(size(stackB));
    stackBr(ybegin_A:yend_A,xbegin_A:xend_A,:,:) = stackB(ybegin_B:yend_B,xbegin_B:xend_B,:,:);
    stackBr(stackBr == 0) = mode(stackB(:));
    
end