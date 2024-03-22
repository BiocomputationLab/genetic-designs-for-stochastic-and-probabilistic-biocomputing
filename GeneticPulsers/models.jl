using DifferentialEquations
using Catalyst
using ModelingToolkit
using BenchmarkTools

const H = 4

function Regulated(by, n; name)
    by = ParentScope(by)
    return @reaction_network $name begin
        @parameters β₁=2 β₂=1 α₁=1 α₂=1/10 K=β₁*β₂/α₁/α₂/10
        @species m(t)=0 p(t)=0
        β₁ / (1 + (K/$by)^$n), 0 --> m
        β₂, m --> p + m
	α₁, m --> 0
	α₂, p --> 0
    end
end

function DoubleRegulated(by, and_by, n, m; name)
    by = ParentScope(by)
    and_by = ParentScope(and_by)
    return @reaction_network $name begin
        @parameters β₁=2 β₂=1 α₁=1 α₂=1/10 K=β₁*β₂/α₁/α₂/10
        @species m(t)=0 p(t)=0
        β₁ / (1 + (K/$by)^$n) / (1 + (K/$and_by)^$m), 0 --> m
	β₂, m --> p + m
	α₁, m --> 0
	α₂, p --> 0
    end
end

function IFL(X, xy, yz, xz; name)
    @named Y = Regulated(ParentScope(X), xy)
    @named Z = DoubleRegulated(ParentScope(X), Y.p, xz, yz)
    return compose(@reaction_network, [Y, Z]; name=name)
end

IFL3(X; name) = IFL(X, H, H, -H; name=name)
IFL4(X; name) = IFL(X, -H, H, H; name=name)

function Repeaters(N; name)
    @variables t
    @species X(t)=0
    syss = [IFL3(X; name=:FFL1)]
    for i in 2:N
        push!(syss, IFL3(syss[i-1].Z.p; name=Symbol(:FFL, i)))
    end
    sys = ReactionSystem(Reaction[], t, [X], []; name=name)
    return compose(sys, syss)
end

function pulses(times, width, value, idx)
    on! = let idx=idx, value=value
	function (integrator)
	    integrator.u[idx] = value
	    reset_aggregated_jumps!(integrator)
	end
    end
    off! = let idx=idx
        function (integrator)
            integrator.u[idx] = 0
	    reset_aggregated_jumps!(integrator)
        end
    end
    return CallbackSet(
	PresetTimeCallback(times, on!, save_positions=(true, true)), 
	PresetTimeCallback(times .+ width, off!, save_positions=(true, true)),
    )
end

function repeater(a, N::Int=2)
    @named system = Repeaters(N)
    params = reduce(vcat, [[s.Y.β₁ => a, s.Z.β₁ => a] for s in system.systems])
    problem = DiscreteProblem(system, [], (0., 280.0), params)
    jumpproblem = JumpProblem(system, problem, Direct())
    callback = pulses([10.], 20., a, 1)
    return solve(jumpproblem, SSAStepper(); callback=callback)
end     

function does_pulse(solution, idx, threshold)
    return any(u -> u[idx] >= threshold, solution.u)
end

function error_rate(a, N)
    @named system = Repeaters(1)
    params = reduce(vcat, [[s.Y.β₁ => a, s.Z.β₁ => a] for s in system.systems])
    problem = DiscreteProblem(system, [], (0., 280.0), params)
    jumpproblem = JumpProblem(system, problem, Direct())
    callback = pulses([10.], 20., Int(floor(2 * a)), 1)
    e = falses(N)
    Threads.@threads for i in 1:N
        solution = solve(jumpproblem, SSAStepper(); callback=callback)
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
