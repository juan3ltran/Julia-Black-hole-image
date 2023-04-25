
include("config.jl")
module cm
using ..config
using DifferentialEquations
export spectrum,g,geodesics!,image_plane,photon_coords,thin_disk,metric, geodesics,M,D,iota,screen_side,n_pixels,detector,R_min,R_max,acc_structure,filename
export initCond, Photon, initial_conditions,geo_integ, image, create_image, create_photons

function initCond(x,k,metric=config.metric,M=config.M)
    g_tt, g_rr, g_thth, g_phph = metric(x, M)   
    
    k_t = g_tt*k[1] 
    k_r = g_rr*k[2]
    k_th = g_thth*k[3]
    k_phi = g_phph*k[4]
    
    return [x[1], x[2], x[3], x[4], k_t, k_r, k_th, k_phi]
end

mutable struct Photon
    alpha
    beta
    i
    j
    xin
    kin
    fP
    iC
    function Photon(alpha, beta)
        i=1
        j=1
        xin=1
        kin=1
        fP=1
        iC=1
        new(alpha,beta,i,j,xin,kin,fP,iC)
    end
end

function initial_conditions(P::Photon)
    P.iC=initCond(P.xin,P.kin)    
end

function geo_integ(p::Photon,in_edge,out_edge)
    λ=(0,-200)
    prob=ODEProblem(config.geodesics!,p.iC,λ,config.M)
    sol= solve(prob,Vern7(),saveat = 0.105)
    indx=length(sol[2,:])
    p.fP = [0,0,0,0,0,0,0,0]
    for i in range(1,indx)
        if sol[2,i]<2*config.M +1e-5
            indx= i
            break
        elseif abs(sol[2,i]*cos(sol[3,i])) < 1e-1
            if sol[2,i] > in_edge && sol[2,i] < out_edge
                indx = i
                p.fP = sol[i]
                break
            
            end
        end
    end
end

mutable struct image
    photon_list
    image_data
    detector

    function image()
        new([],1,1)        
    end
end

function create_photons(x::image,detector::config.image_plane)
    x.detector=detector
    i=1
    for a in x.detector.alphaRange
        j = 1
        for b in x.detector.betaRange
            p = Photon(a, b)
            p.xin, p.kin = config.photon_coords(x.detector,a, b) 
            p.j, p.i = i, j
            initial_conditions(p)
            push!(x.photon_list,p)
            j += 1
        end
        i += 1
    end
    
end

function create_image(x::image,acc_structure::config.thin_disk)
    x.image_data=zeros(Float64,x.detector.numPixels,x.detector.numPixels)
    photon=1
    for p in x.photon_list
        geo_integ(p,acc_structure.Rmin,acc_structure.Rmax)
        x.image_data[p.i,p.j]=config.spectrum(acc_structure,p.fP[2])
        println("photon #$(photon)")
        photon+=1
    end    
end

end




