# Directory to the mesh
set dir /home/mroe734/Documents/phd/mesh/rat

# Read elements and nodes of the surface mesh
gfx read node AWA015_PTA_1_Rec_Trans_surface_mesh.exnode 
gfx read elem AWA015_PTA_1_Rec_Trans_surface_mesh.exelem
gfx def faces egroup uterus

# Create scene and spectrum
gfx create win

# Create uterus mesh and colour it
gfx modify g_element /uterus surfaces material tissue

gfx modify window 1 image view_all


