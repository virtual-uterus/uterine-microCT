# Directory to the mesh
set dir /home/mroe734/Documents/phd/mesh/rat

# Read elements and nodes of the volumetric mesh
gfx read node AWA015_PTA_2_Ova_Rec_Trans_volumetric_mesh_annotated.exnode 
gfx read elem AWA015_PTA_2_Ova_Rec_Trans_volumetric_mesh_annotated.exelem
gfx def faces egroup uterus

# Create scene and spectrum
gfx create win
gfx create spectrum thickness
gfx define font large "20 decorative normal normal"

# Create uterus mesh and colour it
gfx modify g_element /uterus surfaces data thickness spectrum thickness
gfx modify spectrum thickness autorange

# Create colour bar
gfx create colour_bar spectrum thickness number_format %.1e font large centre -1.1 0 0.5
gfx modify g_element /uterus point glyph colour_bar spectrum thickness LOCAL NORMALISED_WINDOW_FIT_LEFT scale_factors 0.8

gfx modify window 1 image view_all

# For full mesh views
#gfx modify window 1 image rotate 1 0 0 90
#gfx modify window 1 image rotate 0 0 1 -90

# For clipped view of the mesh
gfx modify window 1 view near_clipping_plane 1545
gfx modify window 1 image rotate 0 0 1 -90
gfx modify window 1 image rotate 1 0 0 90
