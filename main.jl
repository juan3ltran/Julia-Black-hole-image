include("common.jl")

imagen=cm.image()
cm.create_photons(imagen,config.detector)
cm.create_image(imagen,config.acc_structure)
imagen.image_data
using Plots
pythonplot()
heatmap(imagen.image_data,size=(600,600),colorbar=false, c= :hot)