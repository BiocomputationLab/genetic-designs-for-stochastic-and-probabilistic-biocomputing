include("models.jl")

# a population of repeaters, responding to a bunch of trains of the same arrival rate

using Distributions

b₁ = 0.5
a₁ = log(2) / 5
b₂ = 0.25
a₂ = log(2) / 30
k = b₁ * b₂ / a₁ / a₂ / 10

ps = @parameters β₁=b₁ β₂=b₂ α₁=a₁ α₂=a₂ K=k

function pulse_train(amplitude, n, T)
    times = sort(rand(Uniform(0.0, T), rand(Poisson(n))))
    return pulses(times, amplitude, 1)
end

function repeater_problem(T)
    system = complete(Repeaters(1))
    problem = DiscreteProblem(system, [], (0., T))
    return JumpProblem(system, problem, Direct(), save_positions=(false, false))
end

function Normal(strength, T)
    ps = @parameters β₁=strength β₂=b₂ α₁=a₁ α₂=a₂
    @variables t
    vs = @species m(t)=0 p(t)=0
    rxs = [
        Reaction(β₁, nothing, [m]),
        Reaction(β₁ / 100, nothing, [m]),
        Reaction(β₂, [m], [p, m]),
        Reaction(α₁, [m], nothing),
        Reaction(α₂, [p], nothing),
    ]
    return ReactionSystem(rxs, t, vs, ps; name=name)
end

function normal_problem(strength, T)
    system = complete(Normal(strength, T))
    problem = DiscreteProblem(system, [], (0., T))
    return JumpProblem(system, problem, Direct(), save_positions=(false, false))
end

function X(strength, T, N)
    population = zeros(Int(T) + 1)
    problem = normal_problem(T)
    for i in 1:N
        println("cell $i")
        @time solution = solve(problem, SSAStepper(); saveat=1)
        population .= population .+ getindex.(solution.u, 2)
    end
    return population
end

function Z(amplitude, n, T, N)
    population = zeros(Int(T) + 1)
    problem = repeater_problem(T)
    for i in 1:N
        println("cell $i")
        train = pulse_train(amplitude, n, T)
        @time solution = solve(problem, SSAStepper(); callback=train, saveat=1)
        population .= population .+ getindex.(solution.u, 6)
    end
    return population
end

ncells = 32
T = 10_000

pulsed_population = Z(40, 10, T, ncells)
plt = plot(pulsed_population ./ ncells)
normal_population = X(, 10, T, ncells)
plot!(plt, normal_population ./ ncells)
display(plt)
@show mean(pulsed_population) / std(pulsed_population)
@show mean(normal_population) / std(normal_population)

