
vol=maskfat;
vol=rot90(vol,3);
vol=fliplr(vol);


mha_write_volume('testimg.mha', vol)
