function[img] = read_3d_tif(filename, xlen, ylen, zlen)

    img = zeros(xlen, ylen, zlen);
    for z = 1:zlen
    	img(:,:,z) = imread(filename, z);
    end

end
