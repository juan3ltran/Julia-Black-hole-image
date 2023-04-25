module imageplane

export image_plane
export photon_coords

mutable struct image_plane
    D
    iota
    s_side
    numPixels
    alphaRange
    betaRange
    function image_plane(D::Any,iota::Any,s_side::Any,numPixels::Any)
        if numPixels%2!=0
            numPixels+=1
        end
        aplhaRange=Array(LinRange(-s_side,s_side,numPixels))
        betaRange=Array(LinRange(-s_side,s_side,numPixels))
        new(D,iota,s_side,numPixels,aplhaRange,betaRange)
    end
end

function photon_coords(ip::image_plane,alpha,beta,freq=1)
    r = sqrt(alpha^2 + beta^2 + ip.D^2)
    theta = acos((beta*sin(ip.iota) + ip.D*cos(ip.iota))/r)
    phi = atan(alpha/(ip.D*sin(ip.iota) - beta*cos(ip.iota)))


    xin = [0., r, theta, phi]
                    
    w0 =  freq    
    aux = alpha^2 + (-beta*cos(ip.iota) + ip.D*sin(ip.iota))^2 
    kr =  (ip.D/r)*w0   
    ktheta = (w0/sqrt(aux))*(-cos(ip.iota) 
                + (beta*sin(ip.iota) + ip.D*cos(ip.iota))*(ip.D/(r^2)))  
    kphi = - alpha*sin(ip.iota)*w0/aux     
    kt = sqrt(kr^2 + r^2 * ktheta^2 + r^2*(sin(theta))^2 *kphi^2)       
    
    kin = [kt, kr, ktheta, kphi]
    return xin, kin

end






end