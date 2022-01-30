export VNum, ϕ, proj, sangle, cartesian, dim_count

struct VNum{T}

    mag::T
    angles::Vector{T}

    function VNum(mag::T, angles::Vector{T}) where T
        for i in eachindex(angles)
            if mag < 0
                angles[i] += pi
            end
            angles[i] = (angles[i] < 0 ? angles[i] + 2 * pi : angles[i]) % (2 * pi)
        end
        new{T}(abs(mag), angles)
    end

    VNum(mag::T, angles...) where T = VNum(mag, T[item for item in angles])


    function VNum(coordinates::Vector{T}) where T
        mag::T = 0
        for coord in coordinates
            mag += coord ^ 2
        end

        partialp = mag
        mag = sqrt(mag)

        addPI = coordinates[end] < 0
        deleteat!(coordinates, length(coordinates))

        for i in 1:length(coordinates)
            value = acos(coordinates[i] / sqrt(partialp))
            partialp -= coordinates[i] ^ 2
            coordinates[i] = value
        end

        if addPI
            coordinates[end] = 2 * pi - coordinates[end]
        end

        ϕ(mag, coordinates)
    end

    "Solution to adding vectors"
    function Base.:+(a::VNum, b::VNum)
        c = sqrt(abs(a) ^ 2 + abs(b) ^ 2 + 2 * proj(a, b))
        ϕ(c, broadcast((x, y) -> asin((abs(b) / c) * sin(pi - (y - x))) + x, a.angles, b.angles))
    end
    
    "Solution to subtracting vectors"
    function Base.:-(c::VNum{T}, a::VNum{T}) where T
        b = sqrt(abs(c) ^ 2 + abs(a) ^ 2 - 2 * proj(c, a))
        (b == 0) && return ϕ(T(0), zeros(T, dim_count(c)))
        ϕ(b, broadcast((x, y) -> (pi - asin(sin(x - y) * abs(c) / b) + y), c.angles, a.angles))
    end 

    Base.:*(a::VNum, b::VNum) = ϕ(abs(a) * abs(b), a.angles + b.angles)
    Base.:/(a::VNum, b::VNum) = ϕ(abs(a) / abs(b), a.angles - b.angles)
    

    Base.:^(a::VNum, b::Number) = ϕ(abs(a) ^ b, a.angles .* b)
    Base.:+(a::VNum, b::Number) = a + ϕ(b)
    Base.:-(a::VNum, b::Number) = a - ϕ(b)
    Base.:*(a::VNum, b::Number) = a * ϕ(b)
    Base.:/(a::VNum, b::Number) = a / ϕ(b)
    
    Base.:/(b::Number, a::VNum) = ϕ(b) / a
    Base.:*(b::Number, a::VNum) = ϕ(b) * a
    Base.:-(b::Number, a::VNum) = ϕ(b) - a
    Base.:+(b::Number, a::VNum) = a + ϕ(b)
    Base.:^(b::Number, a::VNum) = ϕ(b) ^ a
    
    Base.string(a::VNum) = "ϕ($(a.mag), $(a.angles))\t≈\tCartesian$(cartesian(a))" 
    Base.show(io::IO, a::VNum) = show(io, string(a))
    Base.print(io::IO, a::VNum) = print(io, string(a))
end

ϕ(n::T) where T = VNum(n, n < 0 ? T[pi] : T[0])
ϕ(mag, angles::Vector) = VNum(mag, angles)
ϕ(mag, angles...) = VNum(mag, angles...)
sangle(a::VNum) = a.angles[1]
Base.abs(a::VNum) = a.mag
proj(a::VNum, b::VNum) = abs(a) * abs(b) * cos(sangle(b) - sangle(a))
dim_count(a::VNum) = length(a.angles)

function cartesian(p::VNum{T}) where T
     coordinates = zeros(T, dim_count(p) + 1)
     coordinates[1] = p.mag * cos(p.angles[1])

     sin_product = p.mag * sin(p.angles[1])
     for i in 2:(length(coordinates) - 1)
         coordinates[i] = sin_product * cos(p.angles[i])
         sin_product *= sin_product * sin(p.angles[i])
     end

     coordinates[end] = sin_product
     coordinates
end