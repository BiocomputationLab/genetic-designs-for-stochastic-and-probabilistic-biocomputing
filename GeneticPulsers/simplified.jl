using Distributions
using StatsBase

function point_process(rate, T)
    N = rand(Poisson(rate * T))
    return sort(rand(Uniform(0, T), N))
end

function match(a::Vector, b::Vector)
    z1 = falses(length(a))
    z2 = falses(length(b))
    z1[1] = any(x -> 0 < x < a[1], b)
    z2[1] = any(x -> 0 < x < b[1], a)
    for i in 2:length(a)
        z1[i] = any(x -> a[i - 1] < x < a[i], b)
    end
    for i in 2:length(b)
        z2[i] = any(x -> b[i - 1] < x < b[i], a)
    end
    return a[z1], b[z2]
end

function difference(a::Vector, b::Vector)
    z1 = falses(length(a))
    z2 = falses(length(b))
    z1[1] = !any(x -> 0 < x < a[1], b)
    z2[1] = !any(x -> 0 < x < b[1], a)
    for i in 2:length(a)
        z1[i] = !any(x -> a[i - 1] < x < a[i], b)
    end
    for i in 2:length(b)
        z2[i] = !any(x -> b[i - 1] < x < b[i], a)
    end
    return a[z1], b[z2]
end

intervals(p) = p .- [0; p[1:end-1]]

function match(r1::Real, r2::Real)
    p1 = point_process(r1, 1)
    p2 = point_process(r2, 1)
    p3, p4 = match(p1, p2)
    return p1, p2, p3, p4, sort(vcat(p3, p4))
end

function difference(r1::Real, r2::Real)
    p1 = point_process(r1, 1)
    p2 = point_process(r2, 1)
    p3, p4 = difference(p1, p2)
    return p1, p2, p3, p4, sort(vcat(p3, p4))
end

function match(r1::Real, r2::Real, N::Int)
    p = zeros(Int, 5, N)
    for i in 1:N
        p[:, i] .= length.(match(r1, r2))
    end
    return p
end

function difference(r1::Real, r2::Real, N::Int)
    p = zeros(Int, 5, N)
    for i in 1:N
        p[:, i] .= length.(difference(r1, r2))
    end
    return p
end

x1, x2 = point_process(1000, 1), point_process(2000, 1)
y1, y2 = match(x1, x2)
z1, z2 = match(y1, x2)
plt = plot(vline.([x1, x2, z1, z2], xlims=(0, 1))..., layout=(8, 1))
matchers(r1, r2) = r1 * r2 / (r1 + r2)

function stats(p)
    @show mean(intervals(p))
    @show mean(intervals(p))^2
    @show var(intervals(p))
    return nothing
end

