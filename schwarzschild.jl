
module schwarzschild
export geodesics!
export g

function geodesics!(du,u,M,λ)
    du[1] = u[5]*u[2]^2/(u[2]^2 - 2*M*u[2])
    du[2] = (1 - 2*M/u[2])*u[6]
    du[3] = u[7]/u[2]^2
    du[4] = u[8]/((u[2]*sin(u[3]))^2)    
    du[5] = 0.
    du[6] = -M*(u[6]/u[2])^2 + u[7]^2/u[2]^3 +u[8]^2/((u[2]^3)*sin(u[3])^2)-M*(u[5]/(u[2]-2*M))^2 
    du[7]= (cos(u[3])/sin(u[3])^3)*(u[8]/u[2])^2
    du[8]= 0.    
end

function g(x,M)
    gtt = -(1 - 2*M/x[2])
    grr = 1/(1- 2*M/x[2])
    gthth = x[2]^2
    gphph = (x[2]*sin(x[3]))^2
    
    return [gtt, grr, gthth, gphph]   
end


end