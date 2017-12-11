function gain = CreateGain(flat_image,offset,spectra_cube,X,Y,PL,u_sn,flat_surf,u_fp)

model_image = forward_model(u_sn,flat_surf,u_fp,1,X,Y,PL,spectra_cube);

flat_image = double(demosaic(uint16(double(flat_image) - double(offset)),'gbrg'));
flat_image(lt(flat_image,1)) = 1;

gain = model_image./flat_image;