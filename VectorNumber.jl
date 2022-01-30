export VNum, ϕ, proj, sangle, cartesian, dim_count

struct VNum{T}

    mag::T
    angles::Vector{T}

    VNum(mag::T, angles::Vector{T}) where T = mag < 0 ? new{T}(-mag, (angles .+ pi) .% (2 * pi)) : new{T}(mag, (angles) .% (2 * pi))
    VNum(mag::T, angles...) where T = VNum(mag, collect(angles))


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


    Base.:+(a::VNum, b::VNum) = ϕ(sqrt(abs(a) ^ 2 + abs(b) ^ 2 + 2 * proj(a, b)), broadcast((x, y) -> (x * abs(a) + y * abs(b) / (abs(a) + abs(b)), a.angles, b.angles)))
    Base.:-(a::VNum, b::VNum) = ϕ(sqrt(abs(a) ^ 2 - abs(b) ^ 2 - 2 * proj(a, b)), broadcast((x, y) -> (x * abs(a) - y * abs(b) / (abs(a) - abs(b)), a.angles, b.angles)))
    Base.:*(a::VNum, b::VNum) = ϕ(abs(a) * abs(b), broadcast((x, y) -> x + y, a.angles, b.angles))
    Base.:/(a::VNum, b::VNum) = ϕ(abs(a) / abs(b), broadcast((x, y) -> x - y, a.angles, b.angles))
    
    Base.:^(a::VNum, b::Number) = ϕ(abs(a) ^ b, a.angles .^ b)
    Base.:^(a::VNum, b::Number) = ϕ(abs(a) * b, a.angles)
    Base.:^(b::Number, a::VNum) = ϕ(abs(a) * b, a.angles)
    Base.:+(a::VNum, b::Number) = a + ϕ(b)
    Base.:+(a::Number, b::VNum) = a + ϕ(b)
    Base.:-(a::VNum, b::Number) = a - ϕ(b)
    Base.:-(a::Number, b::VNum) = ϕ(b) - a
    Base.:*(a::VNum, b::Number) = a * ϕ(b)
    Base.:*(a::Number, b::VNum) = ϕ(b) * a
    Base.:/(a::VNum, b::Number) = a / ϕ(b)
    Base.:/(a::Number, b::VNum) = ϕ(b) / a
end

ϕ(n) = ϕ(n, n < 0 ? [pi] : [0])
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