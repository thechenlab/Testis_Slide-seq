function[] = write_3d_tif(filename, img)

    imwrite(uint16(img(:,:,1)),filename);
    for z=2:size(img,3)
        imwrite(uint16(img(:,:,z)),filename,'WriteMode','append');
    end
	
end
