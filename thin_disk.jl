module thindisk

export thin_disk
export spectrum

mutable struct thin_disk
    Rmin
    Rmax
end

function spectrum(x::thin_disk,r)
    m = (1-0)/(x.Rmin - x.Rmax)
    intensity = m * (r - x.Rmax)
    if r>x.Rmin && r < x.Rmax
        return intensity
    else
        return 0
    end
end

end