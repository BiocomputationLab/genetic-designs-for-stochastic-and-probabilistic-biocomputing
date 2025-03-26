using DifferentialEquations
using Catalyst
using ModelingToolkit
using Plots
using StatsBase

const H = 4

b₁ = 0.5
a₁ = log(2) / 5
b₂ = 0.25
a₂ = log(2) / 30
k = b₁ * b₂ / a₁ / a₂ / 10

ps = @parameters β₁=b₁ β₂=b₂ α₁=a₁ α₂=a₂ K=k

function Regulated(by, n; name)
    by = ParentScope(by)
    @variables t
    vs = @species m(t)=0 p(t)=0
    rxs = [
        Reaction(β₁ / (1 + (K/by)^n), nothing, [m]),
        Reaction(β₁ / 100, nothing, [m]),
        Reaction(β₂, [m], [p, m]),
        Reaction(α₁, [m], nothing),
        Reaction(α₂, [p], nothing),
    ]
    return ReactionSystem(rxs, t, vs, ps; name=name)
end

function DoubleRegulated(by, and_by, n1, n2; name)
    by = ParentScope(by)
    and_by = ParentScope(and_by)
    @variables t
    vs = @species m(t)=0 p(t)=0
    rxs = [
        Reaction(β₁ / (1 + (K/by)^n1) / (1 + (K/and_by)^n2), nothing, [m]),
        Reaction(β₁ / 100, nothing, [m]),
	Reaction(β₂, [m], [m, p]),
        Reaction(α₁, [m], nothing),
        Reaction(α₂, [p], nothing),
    ]
    return ReactionSystem(rxs, t, vs, ps; name=name)
end

function IFL(X, xy, yz, xz; name)
    @named Y = Regulated(ParentScope(X), xy)
    @named Z = DoubleRegulated(ParentScope(X), Y.p, xz, yz)
    @variables t
    return compose(ReactionSystem(Reaction[], t, [X], []; name=name), [Y, Z])
end

IFL1(X; name) = IFL(X, H, -H, H; name=name)
IFL2(X; name) = IFL(X, -H, -H, -H; name=name)
IFL3(X; name) = IFL(X, H, H, -H; name=name)
IFL4(X; name) = IFL(X, -H, H, H; name=name)

function Repeaters(N; name=gensym())
    @variables t
    @species X(t)=0
    syss = [IFL4(X; name=:FFL1)]
    for i in 2:N
        push!(syss, IFL4(syss[i-1].Z.p; name=Symbol(:FFL, i)))
    end
    sys = ReactionSystem([Reaction(α₂, [X], nothing)], t, [X], [α₂]; name=name)
    return compose(sys, syss)
end

function pulses(times, value, idx)
    on! = let idx=idx, value=value
	function (integrator)
	    integrator.u[idx] = integrator.u[idx] + value
	    reset_aggregated_jumps!(integrator)
	end
    end
    return CallbackSet(PresetTimeCallback(times, on!, save_positions=(false, false)))
end

function repeater(a, N::Int=2)
    system = Repeaters(N)
    # params = reduce(vcat, [[s.Y.β₁ => a, s.Z.β₁ => a] for s in system.systems])
    system = complete(system)
    problem = DiscreteProblem(system, [], (0., 10000.0))
    jumpproblem = JumpProblem(system, problem, Direct(), save_positions=(false, false))
    callback = pulses(sort(rand(10) .* 10000.0), a, 1)
    return solve(jumpproblem, SSAStepper(); callback=callback, saveat=10)
end     

function does_pulse(solution, idx, threshold)
    return any(u -> u[idx] >= threshold, solution.u)
end

function pulse_detector(solution, idx, lag, threshold, influence)
    y = Float32.([solution(t)[idx] for t in range(solution.t[1], solution.t[end], Int(floor(solution.t[end])))])
    n = length(y)
    signals = zeros(n) # init signal results
    filteredY = copy(y) # init filtered series
    avgFilter = zeros(n) # init average filter
    stdFilter = zeros(n) # init std filter
    avgFilter[lag - 1] = mean(y[1:lag]) # init first value
    stdFilter[lag - 1] = std(y[1:lag]) # init first value

    for i in range(lag, stop=n-1)
        if abs(y[i] - avgFilter[i-1]) > threshold*stdFilter[i-1]
            if y[i] > avgFilter[i-1] && y[i] > 5
                signals[i] += 1 # postive signal
            end
            # Make influence lower
            filteredY[i] = influence*y[i] + (1-influence)*filteredY[i-1]
        else
            signals[i] = 0
            filteredY[i] = y[i]
        end
        avgFilter[i] = mean(filteredY[i-lag+1:i])
        stdFilter[i] = std(filteredY[i-lag+1:i])
    end
    return signals
end

function does_pulse(solution, idx)
    u = [solution(t)[idx] for t in range(solution.t[1], solution.t[end], 256)]
    ubar = mean(u)
    ustd = std(u)
    @show ubar, ustd
    @show maximum(getindex.(solution.u, idx))
    return maximum(getindex.(solution.u, idx)) >= ubar + 5 * ustd
end

function error_rate(a, N)
    @named system = Repeaters(1)
    params = reduce(vcat, [[s.Y.β₁ => a, s.Z.β₁ => a] for s in system.systems])
    system = complete(system)
    problem = DiscreteProblem(system, [], (0., 280.0), [])
    jumpproblem = JumpProblem(system, problem, Direct(), save_positions=(false, false))
    callback = pulses([100.], 50., Int(floor(2 * a)), 1)
    e = falses(N)
    Threads.@threads for i in 1:N
        solution = solve(jumpproblem, SSAStepper(); callback=callback, saveat=1)
        e[i] = !does_pulse(solution, 5, a)
    end
    return sum(e) / N
end

function Fanout(; name)
    @variables t
    @species X(t)=0
    @named Z1 = Regulated(X, 2)
    @named Z2 = Regulated(X, 4)
    @named Z3 = Regulated(X, 1)
    sys = ReactionSystem(Reaction[], t, [X], []; name=name)
    return compose(sys, [Z1, Z2, Z3])
end

function conventional(a, b::Int)
    @named system = Conventional()
    params = [s.β₁ => a for s in system.systems]
    problem = DiscreteProblem(system, [], (0., 3000.0), params)
    jumpproblem = JumpProblem(system, problem, Direct(), save_positions=(false, false))
    callback = pulses([0.], 9999., b, 1)
    return solve(jumpproblem, SSAStepper(); callback=callback, saveat=10)
end

function PulseFanout(N::Int; name)
    @variables t
    @species X(t)=0
    sys = ReactionSystem(Reaction[], t, [X], []; name=name)
    return compose(sys, [IFL3(X; name=Symbol(:Z, i)) for i in 1:N])
end

function pulsefanout(as, b::Int)
    @named system = PulseFanout(length(as))
    params = reduce(vcat, [s.Z.β₁ => a, s.Y.β₁ => a] for (s, a) in zip(system.systems, as))
    problem = DiscreteProblem(system, [], (0., 3000.0), params)
    jumpproblem = JumpProblem(system, problem, Direct(), save_positions=(false, false))
    callback = pulses(range(10, 3000, 10), 50., b, 1)
    return solve(jumpproblem, SSAStepper(); callback=callback, saveat=10)
end
