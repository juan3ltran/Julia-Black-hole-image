
module config
include("schwarzschild.jl")
include("thin_disk.jl")
include("image_plane.jl")
using .thindisk
using .imageplane
using .schwarzschild


export  spectrum,g,geodesics!,image_plane,photon_coords,thin_disk

export metric, geodesics,M,D,iota,screen_side,n_pixels,detector,R_min,R_max,acc_structure,filename
#metric
metric=schwarzschild.g
geodesics=schwarzschild.geodesics!

#BH parameters
M=1
# DETECTOR PARAMETERS 
D = 100*M
iota = pi/2.2
screen_side = 22*M
n_pixels = 1500
detector = image_plane(D, iota, screen_side,  n_pixels)


# ACCRETION STRUCTURE 
R_min = 6*M
R_max = 20*M
acc_structure = thindisk.thin_disk(R_min, R_max)

end
