export Vec, VecF64

struct Vec{N, T}
    coords::SVector{N, T}
    Vec{N, T}(x) where {N, T} = new(x)
end
const VecF64 = Vec{N, Float64} where {N}

Vec(x...) = Vec{length(x), Base.promote_typeof(x...)}(x)
Vec(x::AbstractArray{T}) where {T} = Vec{length(x), T}(x)

(::Type{VecF64})(x...) = Vec{length(x), Float64}(x)
(::Type{VecF64})(p::Vec) = Vec{length(p), Float64}(p.coords)

Base.zero(::Type{Vec{N, T}}) where {N, T} = Vec(ntuple(_ -> zero(T), N)...)
Base.zero(::Vec{N, T}) where {N, T} = Vec(ntuple(_ -> zero(T), N)...)

Base.eltype(::Vec{N, T}) where {N, T} = T

Base.length(::Vec{N}) where {N} = N

Base.getindex(x::Vec, i::Integer) = x.coords[i]

function Base.getproperty(p::Vec, sym::Symbol)
    return if sym == :x
        p.coords.x
    elseif sym == :y
        p.coords.y
    elseif sym == :z
        p.coords.z
    else
        getfield(p, sym)
    end
end

Base.isless(p::Vec{2}, q::Vec{2}) = p.y < q.y || (p.y == q.y && p.x > q.x)

Base.:(*)(p::Vec, x::Number) = Vec(p.coords * x)
Base.:(*)(x::Number, p::Vec) = p * x

LinearAlgebra.normalize(p::Vec) = Vec(normalize(p.coords))
