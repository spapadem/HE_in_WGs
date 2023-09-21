N = 1;

using VTUFileHandler

# uimag= VTUFile();
# ureal= VTUFile("wgfem_real.vtu");

set_uncompress_keywords(["xRamp"]) # uncrompress data field xramp
set_interpolation_keywords(["xRamp"]) # apply math operators to xramp
vtu = VTUFile("wgfem_imag.vtu"); # read the vtu