module TestUtils
using Serialization

using Test

function fieldvalues(obj)
    [getfield(obj, f) for f in fieldnames(typeof(obj))]
end

deep_compare(a, b; rtol::Real=sqrt(eps()), atol::Real=0.) = a == b # base case: use ==
deep_compare(a::Union{AbstractFloat, Complex}, b::Union{AbstractFloat, Complex}; rtol::Real=sqrt(eps()), atol::Real=0.) = isapprox(a, b) # for floats, use isapprox
deep_compare(a::Union{Array{AbstractFloat}, Array{Complex}}, b::Union{Array{AbstractFloat}, Array{Complex}}; rtol::Real=sqrt(eps()), atol::Real=0.) = isapprox(a, b) # for arrays of floats, use isapprox element-wise
deep_compare(a::AbstractArray, b::AbstractArray; rtol::Real=sqrt(eps()), atol::Real=0.) = all(deep_compare.(a, b; rtol, atol)) # for arrays of other types, recurse
deep_compare(a::Tuple, b::Tuple; rtol::Real=sqrt(eps()), atol::Real=0.) = all(deep_compare.(a, b; rtol, atol)) # for tuples, recurse
deep_compare(a::T, b::T; rtol::Real=sqrt(eps()), atol::Real=0.) where {T <: Any} = deep_compare(fieldvalues(a), fieldvalues(b); rtol, atol) # for composite types, recurse

function test_from_serialized(fn::Function, filename::String)
    argskwargs, out = open(filename, "r") do f
        deserialize(f)
    end

    @test deep_compare([fn(a...; kw...) for (a, kw) in argskwargs], out)
end

function asymmetric_floatvec_to_complexvec(vec_a::Vector{Float32}, vec_b::Vector{Float32})
    shortest = min(length(vec_a), length(vec_b))
    ComplexF32.(vec_a[1:shortest], vec_b[1:shortest])
end

end