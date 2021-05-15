function[img] = readPlane(reader, series, z, c, t)

    reader.setSeries(series - 1);
    plane = reader.getIndex(z - 1, c -1, t - 1) + 1;
    img = bfGetPlane(reader, plane);

end