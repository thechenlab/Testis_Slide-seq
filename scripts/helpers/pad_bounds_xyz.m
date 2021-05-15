function[bounds_pad] = pad_bounds_xy(bounds,pad,zpad,xlen,ylen,zlen)

    for i=1:size(bounds,1)
        
        bounds_pad(i,1) = max(1,bounds(i,1)-pad);
        bounds_pad(i,2) = max(1,bounds(i,2)-pad);
        bounds_pad(i,3) = max(1,bounds(i,3)-zpad);
        bounds_pad(i,4) = min(bounds(i,4)+pad*2,xlen-bounds(i,1)+pad);
        bounds_pad(i,5) = min(bounds(i,5)+pad*2,ylen-bounds(i,2)+pad);
        bounds_pad(i,6) = min(bounds(i,6)+zpad*2,zlen-bounds(i,3)+zpad);
        
    end
end 