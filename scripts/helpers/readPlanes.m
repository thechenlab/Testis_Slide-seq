function[img] = readPlanes(reader, series, zlen, c, t)

    img = zeros(reader.getSizeX,reader.getSizeY,zlen,'uint16')
    for z=1:zlen
        img(:,:,z) = readPlane(reader,series,z,c,t);
    end

end