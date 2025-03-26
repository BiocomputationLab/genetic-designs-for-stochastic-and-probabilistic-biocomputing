include("models.jl")

using Distributions

function PulseCallback(idx, times, amplitude)
    on! = let idx=idx, amplitude=amplitude
        function (integrator)
            integrator.u[idx] = integrator.u[idx] + amplitude
            reset_aggregated_jumps!(integrator)
        end
    end
    return PresetTimeCallback(times, on!, save_positions=(true, true))
end

function Matcher(; name=gensym())
    @variables t
    @species X1(t)=0 X2(t)=0
    @named Y1 = DoubleRegulated(X1, first(@species Y2₊p(t)=0), -H, -H)
    @named Y2 = DoubleRegulated(X2, Y1.p, -H, -H)
    @named Z1 = DoubleRegulated(X1, Y1.p, H, H)
    @named Z2 = DoubleRegulated(X2, Y2.p, H, H)
    sys = ReactionSystem([
        Reaction(α₂, [X1], nothing),
        Reaction(α₂, [X2], nothing),
    ], t, [X1, X2], [α₂]; name=name)
    return compose(sys, [Y1, Y2, Z1, Z2])
end

function matcher(μ, w, p; name=gensym())
    system = Matcher(; name=name)
    system = complete(system)
    problem = DiscreteProblem(system, [], (0., 10000.0))
    jumpproblem = JumpProblem(system, problem, Direct(), save_positions=(false, false))
    npulses = rand(Poisson(μ))
    pulsets = sort(rand(Uniform(0, 10000), npulses))

    train1 = sort(rand(Uniform(0, 10000), npulses))
    train2 = train1[rand(npulses) .< p]
    train2 = train2 .+ rand(Exponential(w), length(train2))
    callback = CallbackSet(
        PulseCallback(1, train1, 100),
        PulseCallback(2, sort(train2), 100),
    )
    return solve(jumpproblem, SSAStepper(), saveat=10, callback=callback)
end

function signal_match(s1, s2)
    T = max(maximum(s1), maximum(s2))
    m = Matcher()
    dp = DiscreteProblem(m, [], (0., T), [])
    jp = JumpProblem(m, dp, Direct(), save_positions=(false, false))
    callback = CallbackSet(PulseCallback(m, 1, s1, 30, 200), PulseCallback(m, 2, s2, 30, 200))
    solution = solve(jp, SSAStepper(), saveat=1, callback=callback)
    x1 = pulse_detector(solution, 1, 200, 6, 0)
    x2 = pulse_detector(solution, 2, 200, 6, 0)
    z1 = pulse_detector(solution, 8, 200, 6, 0)
    z2 = pulse_detector(solution, 10, 200, 6, 0)

    signalplot = plot(solution, label="X1", idxs=[1])
    plot!(signalplot, solution, label="X2", idxs=[2])
    plot!(signalplot, z1, label="Z1")
    plot!(signalplot, z2, label="Z2")
    return solution, signalplot
end

x1 = sort(rand(Uniform(0., 10_000.), 10))
x2 = x1 .+ 200
x3 = x2 .+ randn(length(x1)) .* 50
x4 = x2 .+ randn(length(x1)) .* 100


function pulse_count(sol, idx)
    signal = pulse_detector(sol, idx, 200, 6, 0)
    signal = signal[2:end] .- signal[1:end-1]
    return count(isequal(1), signal)
end



# p1 = plot(sol, idxs=[1, 2, 8, 10])
# p2 = plot(range(0, sol.t[end], Int(floor(sol.t[end]))), pulse_detector(sol, 8, 200, 6, 0.0))
# p3 = plot(p1, p2, layout=(2, 1))
