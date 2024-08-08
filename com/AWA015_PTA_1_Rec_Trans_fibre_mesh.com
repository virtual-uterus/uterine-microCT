# Directory to the mesh
set dir /home/mroe734/Documents/phd/mesh/rat

# Read elements and nodes of the surface mesh
gfx read node AWA015_PTA_1_Rec_Trans_surface_mesh.exnode 
gfx read elem AWA015_PTA_1_Rec_Trans_surface_mesh.exelem
gfx def faces egroup uterus

# Directory to the fibres
set dir /home/mroe734/Documents/phd/microCT/data/AWA015_PTA_1_Rec_Trans/downsampled/ST/binary
gfx read node Streamlines_L4_FB
gfx read elem Streamlines_L4_FB

# Create scene and spectrum
gfx create win
gfx create spectrum angle linear reverse range 0 90
gfx define font large "20 decorative normal normal"

# Create uterus mesh and colour it
gfx modify g_element /uterus surfaces material tissue
gfx modify material tissue alpha 0.5

# Create fibre lines and colour them
gfx modify g_element /streamlines lines data angle spectrum angle line_width 2

# Create colour bar
gfx create colour_bar spectrum angle number_format %.1e divisions 9 font large centre -1.1 0 0.5
gfx modify g_element /streamlines point glyph colour_bar spectrum angle LOCAL NORMALISED_WINDOW_FIT_LEFT scale_factors 0.8

# Display everything
gfx modify window 1 set slow_transparency
gfx modify window 1 image view_all
gfx modify window 1 layout front_side
