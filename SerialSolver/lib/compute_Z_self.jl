function compute_Z_self(R, Cd, w)
    Z_self = zeros(ComplexF64, length(R))
    for cont in range(1,length(R))
        for aux in 1:4
            if R[cont, aux] != 0 && Cd[cont, aux] != 0
                Z_self[cont] += 1 / (1 / R[cont, aux] + 1im * w * Cd[cont, aux])
            elseif R[cont, aux] != 0
                Z_self[cont] += R[cont, aux]
            elseif Cd[cont, aux] != 0
                Z_self[cont] += 1 / (1im * w * Cd[cont, aux])
            end
        end
    end
    return Z_self
end
