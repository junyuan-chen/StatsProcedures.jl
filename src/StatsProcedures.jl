module StatsProcedures

using Base: sym_in
using Graphs: SimpleDiGraph, add_edge!, topological_sort_by_dfs
using MacroTools: @capture, isexpr, postwalk

import Base: ==, show, eltype, firstindex, lastindex, getindex, iterate, length

export StatsStep,
       AbstractStatsProcedure,
       SharedStatsStep,
       PooledStatsProcedure,
       StatsSpec,
       proceed,
       @specset

"""
    StatsStep{Alias, F<:Function, ById}

Specification for a step in an [`AbstractStatsProcedure`](@ref).
An instance of `StatsStep` is callable.

# Parameters
- `Alias::Symbol`: alias of the type for pretty-printing.
- `F<:Function`: type of the function to be called by `StatsStep` with `F.instance`.
- `ById::Bool`: whether arguments from multiple [`StatsSpec`](@ref)s should be grouped by `object-id` or `isequal`.

# Methods
    (step::StatsStep{A,F})(ntargs::NamedTuple; verbose::Union{Bool,IO}=false)

Call an instance of function of type `F` with arguments extracted from `ntargs`
via [`groupargs`](@ref) and [`combinedargs`](@ref).

A message with the name of the `StatsStep` is printed to `stdout`
if a keyword `verbose` takes the value `true`
or `ntargs` contains a key-value pair `verbose=true`.
Alternative `IO` stream can be specified by setting the value of `verbose`.
The value from `ntargs` supersedes the keyword argument
in case both are specified.

## Returns
- `NamedTuple`: named intermediate results.
"""
struct StatsStep{Alias, F<:Function, ById} end

_f(::StatsStep{A,F}) where {A,F} = F.instance
_byid(::StatsStep{A,F,I}) where {A,F,I} = I

"""
    required(s::StatsStep)

Return a tuple of `Symbol`s representing the names of arguments
used to form [`groupargs`](@ref) that do not have defaults.
See also [`default`](@ref) and [`transformed`](@ref).
"""
required(::StatsStep) = ()

"""
    default(s::StatsStep)

Return a `NamedTuple` of arguments with keys showing the names
and values representing the defaults to be used to form [`groupargs`](@ref).
See also [`required`](@ref) and [`transformed`](@ref).
"""
default(::StatsStep) = NamedTuple()

"""
    transformed(s::StatsStep, ntargs::NamedTuple)

Return a tuple of arguments transformed from fields in `ntargs`
to be used to form [`groupargs`](@ref).
See also [`required`](@ref) and [`default`](@ref).
"""
transformed(::StatsStep, ::NamedTuple) = ()

_get(nt::NamedTuple, args::Tuple{Vararg{Symbol}}) =
    map(s->getfield(nt, s), args)

# This method is useful for obtaining default values
_get(nt::NamedTuple, args::NamedTuple) =
    map(s->getfield(ifelse(sym_in(s, keys(nt)), nt, args), s), keys(args))

"""
    groupargs(s::StatsStep, ntargs::NamedTuple)

Return a tuple of arguments that allow classifying multiple `ntargs`s into groups.
Equality (defined by `isequal`) of the returned tuples across `ntargs`s imply that
it is possible to exectute step `s` for only once
to obtain results for these `ntargs`s.

This function is important for [`proceed`](@ref) to work properly.
However, in most cases, there is no need to define new methods
for concrete [`StatsStep`](@ref)s.
Instead, one should define methods for
[`required`](@ref), [`default`](@ref) or [`transformed`](@ref).
See also [`combinedargs`](@ref).
"""
groupargs(s::StatsStep, @nospecialize(ntargs::NamedTuple)) =
    (_get(ntargs, required(s))..., _get(ntargs, default(s))...,
    transformed(s, ntargs)...)

"""
    combinedargs(s::StatsStep, allntargs::Any)

Return a tuple of arguments obtained by combining a collection of arguments
across multiple specifications.

The element type of `allntargs` can be assumed to be `NamedTuple`.
This function allows combining arguments that differ
across specifications in the same group classified based on [`groupargs`](@ref)
into objects that are accepted by the call of `s`.
See also [`proceed`](@ref).
"""
combinedargs(::StatsStep, ::Any) = ()

"""
    copyargs(s::StatsStep)

Return a tuple of indices of arguments that are mutable objects
that may be modified in-place by `s`.
This allows making copys of mutable objects that are initially shared
across multiple specifications.
"""
copyargs(::StatsStep) = ()

_parse_verbose(verbose::Bool) = (verbose, stdout)
_parse_verbose(verbose::IO) = (true, verbose)

function (step::StatsStep{A,F})(@nospecialize(ntargs::NamedTuple);
        verbose::Union{Bool,IO}=false) where {A,F}
    haskey(ntargs, :verbose) && (verbose = ntargs.verbose)
    verbose, io = _parse_verbose(verbose)
    verbose && printstyled(io, "Running ", step, "\n", color=:green)
    ret = F.instance(groupargs(step, ntargs)..., combinedargs(step, (ntargs,))...)
    return merge(ntargs, ret)
end

(step::StatsStep)(; kwargs...) = step(NamedTuple(); kwargs...)

show(io::IO, ::StatsStep{A}) where A = print(io, A)

function show(io::IO, ::MIME"text/plain", s::StatsStep{A,F}) where {A,F}
    print(io, A, " (", typeof(s).name.name, " that calls ")
    fmod = F.name.mt.module
    fmod == Main || print(io, fmod, ".")
    print(io, F.name.mt.name, ")")
end

"""
    AbstractStatsProcedure{Alias, T<:NTuple{N,StatsStep} where N}

Supertype for all types specifying the procedure for statistical estimation or inference.

Fallback methods for indexing and iteration are defined for
all subtypes of `AbstractStatsProcedure`.

# Parameters
- `Alias::Symbol`: alias of the type for pretty-printing.
- `T<:NTuple{N,StatsStep}`: steps involved in the procedure.
"""
abstract type AbstractStatsProcedure{Alias, T<:NTuple{N,StatsStep} where N} end

result(::Type{<:AbstractStatsProcedure}, @nospecialize(ntargs::NamedTuple)) = ntargs

length(::AbstractStatsProcedure{A,T}) where {A,T} = length(T.parameters)
eltype(::Type{<:AbstractStatsProcedure}) = StatsStep
firstindex(::AbstractStatsProcedure{A,T}) where {A,T} = firstindex(T.parameters)
lastindex(::AbstractStatsProcedure{A,T}) where {A,T} = lastindex(T.parameters)

function getindex(::AbstractStatsProcedure{A,T}, i) where {A,T}
    fs = T.parameters[i]
    return fs isa Type && fs <: StatsStep ? fs.instance : [f.instance for f in fs]
end

getindex(::AbstractStatsProcedure{A,T}, i::Int) where {A,T} = T.parameters[i].instance

iterate(p::AbstractStatsProcedure, state=1) =
    state > length(p) ? nothing : (p[state], state+1)

"""
    prerequisites(p::AbstractStatsProcedure, s::StatsStep)

Return a tuple of [`StatsStep`](@ref)s that must be completed
prior to executing step `s` when conducting procedure `p`.
If `p` consists of more than one [`StatsStep`](@ref),
this method must be defined for each step in `p`
in order to allow [`pool`](@ref) to share a [`StatsStep`](@ref) across multiple procedures.
"""
prerequisites(::AbstractStatsProcedure{A,Tuple{T}}, ::T) where {A,T} = ()

show(io::IO, ::AbstractStatsProcedure{A}) where A = print(io, A)

function show(io::IO, ::MIME"text/plain", p::AbstractStatsProcedure{A,T}) where {A,T}
    nstep = length(p)
    print(io, A, " (", typeof(p).name.name, " with ", nstep, " step")
    if nstep > 0
        nstep > 1 ? print(io, "s):\n  ") : print(io, "):\n  ")
        for (i, step) in enumerate(p)
            print(io, step)
            i < nstep && print(io, " |> ")
        end
    else
        print(io, ")")
    end
end

"""
    SharedStatsStep{Alias, F<:Function, ById}

A [`StatsStep`](@ref) that is possibly shared by
multiple instances of procedures that are subtypes of [`AbstractStatsProcedure`](@ref).
See also [`PooledStatsProcedure`](@ref) and [`pool`](@ref).

# Fields
- `step::StatsStep{Alias, F, ById}`: the `step` that may be shared.
- `ids::Vector{Int}`: indices of procedures that share `step`.
"""
struct SharedStatsStep{Alias, F<:Function, ById}
    step::StatsStep{Alias, F, ById}
    ids::Vector{Int}
end

SharedStatsStep(step::StatsStep, ids) = SharedStatsStep(step, Int[ids...])

_sharedby(s::SharedStatsStep) = s.ids
_f(::SharedStatsStep{A,F}) where {A,F} = F.instance
_byid(::SharedStatsStep{A,F,I}) where {A,F,I} = I
groupargs(s::SharedStatsStep, @nospecialize(ntargs::NamedTuple)) = groupargs(s.step, ntargs)
combinedargs(s::SharedStatsStep, v::AbstractArray) = combinedargs(s.step, v)
copyargs(s::SharedStatsStep) = copyargs(s.step)

==(x::SharedStatsStep, y::SharedStatsStep) =
    x.step == y.step && x.ids == y.ids

show(io::IO, s::SharedStatsStep) = print(io, s.step)

function show(io::IO, ::MIME"text/plain", s::SharedStatsStep)
    nps = length(s.ids)
    print(io, s.step, " (StatsStep shared by ", nps, " procedure")
    nps > 1 ? print(io, "s)") : print(io, ")")
end

"""
    PooledStatsProcedure

A collection of procedures and shared steps.

An instance of `PooledStatsProcedure` is indexed and iterable among the shared steps
in a way that helps avoid repeating identical steps.
See also [`pool`](@ref).

# Fields
- `procs::Vector{AbstractStatsProcedure}`: a collection of [`AbstractStatsProcedure`](@ref)s.
- `steps::Vector{SharedStatsStep}`: sorted [`SharedStatsStep`](@ref)s.
"""
struct PooledStatsProcedure
    procs::Vector{AbstractStatsProcedure}
    steps::Vector{SharedStatsStep}
end

"""
    pool(ps::Vector{AbstractStatsProcedure})

Return a [`PooledStatsProcedure`](@ref) with identical [`StatsStep`](@ref)s
shared among several procedures in `ps` as [`SharedStatsStep`](@ref)s.
Steps are sorted based on dependencies specified with [`prerequisites`](@ref).
"""
function pool(ps::Vector{AbstractStatsProcedure})
    nps = length(ps)
    nps == 1 && return PooledStatsProcedure(ps, [SharedStatsStep(s, [1]) for s in ps[1]])
    nv = 0
    steps = SharedStatsStep[]
    lookup = Dict{StatsStep, Int}()
    for i in 1:nps
        for s in ps[i]
            iv = get(lookup, s, 0)
            if iszero(iv)
                nv += 1
                push!(steps, SharedStatsStep(s, [i]))
                lookup[s] = nv
            else
                push!(steps[iv].ids, i)
            end
        end
    end
    N = length(steps)
    if N < sum(length, ps)
        dag = SimpleDiGraph(N)
        for p in ps
            for s in p
                for pr in prerequisites(p, s)
                    add_edge!(dag, lookup[pr], lookup[s])
                end
            end
        end
        try
            order = topological_sort_by_dfs(dag)
            return PooledStatsProcedure(ps, steps[order])
        catch
            error("procedures cannot be pooled together as dependencies among steps cannot be resolved")
        end
    else
        return PooledStatsProcedure(ps, steps)
    end
end

length(p::PooledStatsProcedure) = length(p.steps)
eltype(::Type{PooledStatsProcedure}) = SharedStatsStep
firstindex(p::PooledStatsProcedure) = firstindex(p.steps)
lastindex(p::PooledStatsProcedure) = lastindex(p.steps)

getindex(p::PooledStatsProcedure, i) = getindex(p.steps, i)

iterate(p::PooledStatsProcedure, state=1) = iterate(p.steps, state)

==(x::PooledStatsProcedure, y::PooledStatsProcedure) =
    x.procs == y.procs && x.steps == y.steps

show(io::IO, p::PooledStatsProcedure) = print(io, typeof(p).name.name)

function show(io::IO, ::MIME"text/plain", p::PooledStatsProcedure)
    nstep = length(p.steps)
    print(io, typeof(p).name.name, " with ", nstep, " step")
    nstep > 1 ? print(io, "s ") : print(io, " ")
    nps = length(p.procs)
    print(io, "from ", nps, " procedure")
    nps > 1 ? print(io, "s:") : print(io, ":")
    for p in p.procs
        print(io, "\n  ", typeof(p).parameters[1])
    end
end

"""
    StatsSpec{T<:AbstractStatsProcedure}

Specification for a statistical procedure of type `T`.

An instance of `StatsSpec` is callable and
its fields provide all information necessary for conducting the procedure.
An optional name for the specification can be specified.

# Fields
- `name::String`: a name for the specification (takes `""` if not specified).
- `args::NamedTuple`: arguments for the [`StatsStep`](@ref)s in `T` (default values are merged into `args` if not found in `args`).

# Methods
    (sp::StatsSpec{T})(; verbose::Union{Bool,IO}=false, keep=nothing, keepall::Bool=false)

Execute the procedure of type `T` with the arguments specified in `args`.
By default, a dedicated result object for `T` is returned if it is available.
Otherwise, the last value returned by the last [`StatsStep`](@ref) is returned.

## Keywords
- `verbose::Union{Bool,IO}=false`: print the name of each step to `stdout` or the specified `IO` stream when it is called.
- `keep=nothing`: names (of type `Symbol`) of additional objects to be returned.
- `keepall::Bool=false`: return all objects returned by each step.
"""
struct StatsSpec{T<:AbstractStatsProcedure, A<:NamedTuple}
    name::String
    args::A
    StatsSpec(name, T::Type{<:AbstractStatsProcedure}, args::NamedTuple) =
            new{T, typeof(args)}(string(name), args)
end

"""
    ==(x::StatsSpec{T}, y::StatsSpec{T})

Test whether two instances of [`StatsSpec`](@ref)
with the same parameter `T` also have the same field `args`.

See also [`≊`](@ref).
"""
==(x::StatsSpec{T}, y::StatsSpec{T}) where T = x.args == y.args

"""
    ≊(x::NamedTuple, y::NamedTuple)

Test whether two instances of `NamedTuple` contain
the same set of key-value pairs while ignoring the order.

See https://discourse.julialang.org/t/check-equality-of-two-namedtuples-with-order-of-the-fields-ignored
"""
≊(x::NamedTuple{N1,T1}, y::NamedTuple{N2,T2}) where {N1,T1,N2,T2} =
    length(N1) === length(union(N1,N2)) &&
        all(k->getfield(x,k)==getfield(y,k), keys(x))

"""
    ≊(x::StatsSpec{T}, y::StatsSpec{T})

Test whether two instances of [`StatsSpec`](@ref)
with the same parameter `T` also have the field `args`
containing the same sets of key-value pairs
while ignoring the orders.
"""
≊(x::StatsSpec{T}, y::StatsSpec{T}) where T = x.args ≊ y.args

_procedure(::StatsSpec{T}) where T = T

function (sp::StatsSpec{T})(;
        verbose::Union{Bool,IO}=false, keep=nothing, keepall::Bool=false) where T
    args = verbose == false ? sp.args : merge(sp.args, (verbose=verbose,))
    ntall = result(T, foldl(|>, T(), init=args))
    if keepall
        return ntall
    else
        if keep === nothing
            if isempty(ntall)
                return nothing
            else
                return haskey(ntall, :result) ? ntall.result : ntall[end]
            end
        else
            # Cannot iterate Symbol
            if keep isa Symbol
                keep = (keep,)
            else
                eltype(keep)==Symbol ||
                    throw(ArgumentError("expect Symbol or collections of Symbols for the value of option `keep`"))
            end
            in(:result, keep) || (keep = (keep..., :result))
            names = ((n for n in keep if haskey(ntall, n))...,)
            return NamedTuple{names}(ntall)
        end
    end
end

show(io::IO, sp::StatsSpec) = print(io, sp.name=="" ? "unnamed" : sp.name)

# This allows customizing how arguments are printed
_show_args(::IO, ::StatsSpec) = nothing

function show(io::IO, ::MIME"text/plain", sp::StatsSpec{T}) where T
    print(io, sp.name=="" ? "unnamed" : sp.name, " (", typeof(sp).name.name,
        " for ", T.parameters[1], ")")
    _show_args(io, sp)
end

function _count!(objcount::AbstractDict, obj)
    count = get(objcount, obj, 0)
    objcount[obj] = count + 1
end

"""
    proceed(sps::AbstractVector{<:StatsSpec}; kwargs...)

Conduct the procedures for the [`StatsSpec`](@ref)s in `sps`
while trying to avoid repeating identical steps for the [`StatsSpec`](@ref)s.
See also [`@specset`](@ref).

# Keywords
- `verbose::Union{Bool,IO}=false`: print the name of each step to `stdout` or the specified `IO` stream when it is called.
- `keep=nothing`: names (of type `Symbol`) of additional objects to be returned.
- `keepall::Bool=false`: return all objects generated by procedures along with arguments from the [`StatsSpec`](@ref)s.
- `pause::Int=0`: break the iteration over [`StatsStep`](@ref)s after finishing the specified number of steps (for debugging).

# Returns
- `Vector`: results for each specification in the same order of `sps`.

By default, either a dedicated result object for the corresponding procedure
or the last value returned by the last [`StatsStep`](@ref)
becomes an element in the returned `Vector` for each [`StatsSpec`](@ref).
When either `keep` or `keepall` is specified,
a `NamedTuple` with additional objects is formed for each [`StatsSpec`](@ref).
"""
function proceed(sps::Vector{<:StatsSpec};
        verbose::Union{Bool,IO}=false, keep=nothing, keepall::Bool=false, pause::Int=0)
    nsps = length(sps)
    nsps == 0 && throw(ArgumentError("expect a nonempty vector"))
    verbose, io = _parse_verbose(verbose)

    gids = IdDict{AbstractStatsProcedure, Vector{Int}}()
    objcount = IdDict{Any, Int}()
    traces = Vector{NamedTuple}(undef, nsps)
    @inbounds for i in 1:nsps
        push!(get!(Vector{Int}, gids, _procedure(sps[i])()), i)
        args = sps[i].args
        foreach(x->_count!(objcount, x), args)
        traces[i] = args
    end

    steps = pool(collect(AbstractStatsProcedure, keys(gids)))
    tasks_byid = IdDict{Tuple, Vector{Int}}()
    tasks_byeq = Dict{Tuple, Vector{Int}}()
    ntask_total = 0
    step_count = 0
    paused = false
    @inbounds for step in steps
        tasks = _byid(step) ? tasks_byid : tasks_byeq
        verbose && print(io, "Running ", step, "...")
        # Group arguments by objectid or isequal
        for i in _sharedby(step)
            traceids = gids[steps.procs[i]]
            for j in traceids
                push!(get!(Vector{Int}, tasks, groupargs(step, traces[j])), j)
            end
        end
        ntask = length(tasks)
        nprocs = length(_sharedby(step))
        verbose && print(io, "Scheduled ", ntask, ntask > 1 ? " tasks" : " task", " for ",
            nprocs, nprocs > 1 ? " procedures" : " procedure", "...\n")

        for (gargs, ids) in tasks
            # Handle potential in-place operations on mutable objects
            nids = length(ids)
            icopy = copyargs(step)
            if length(icopy) > 0
                gargs = Any[gargs...]
                for i in copyargs(step)
                    a = gargs[i]
                    objcount[a] > nids && (gargs[i] = deepcopy(a))
                end
            end

            ret = _f(step)(gargs..., combinedargs(step, view(traces, ids))...)
            for id in ids
                foreach(x->_count!(objcount, x), ret)
                traces[id] = merge(traces[id], ret)
            end
        end
        ntask_total += ntask
        empty!(tasks)
        verbose && print(io, "  Finished ", ntask, ntask > 1 ? " tasks" : " task", " for ",
            nprocs, nprocs > 1 ? " procedures\n" : " procedure\n")
        step_count += 1
        step_count === pause && (paused = true) && break
    end

    nprocs = length(steps.procs)
    verbose && printstyled(io, "All steps finished (", ntask_total,
        ntask_total > 1 ? " tasks" : " task", " for ", nprocs,
        nprocs > 1 ? " procedures)\n" : " procedure)\n", bold=true, color=:green)
    if !paused
        @inbounds for i in 1:nsps
            traces[i] = result(_procedure(sps[i]), traces[i])
        end
    end

    if keepall
        return traces
    elseif keep === nothing
        return [haskey(r, :result) ? r.result : isempty(r) ? nothing : r[end] for r in traces]
    else
        # Cannot iterate Symbol
        if keep isa Symbol
            keep = (keep,)
        else
            eltype(keep) == Symbol || throw(ArgumentError(
                "expect Symbol or collections of Symbols for the value of option `keep`"))
        end
        in(:result, keep) || (keep = (keep..., :result))
        @inbounds for i in 1:nsps
            names = ((n for n in keep if haskey(traces[i], n))...,)
            traces[i] = NamedTuple{names}(traces[i])
        end
        return traces
    end
end

function _parse!(options::Expr, args)
    for arg in args
        # Assume a symbol means the kwarg takes value true
        if isa(arg, Symbol)
            if arg == :noproceed
                return true
            else
                key = Expr(:quote, arg)
                push!(options.args, Expr(:call, :(=>), key, true))
            end
        elseif isexpr(arg, :(=))
            if arg.args[1] == :noproceed
                return arg.args[2]
            else
                key = Expr(:quote, arg.args[1])
                push!(options.args, Expr(:call, :(=>), key, arg.args[2]))
            end
        else
            throw(ArgumentError("unexpected option $arg"))
        end
    end
    return false
end

"""
    _args_kwargs(exprs)

Return an expression of `Vector{Any}` and an expression of `Dict{Symbol,Any}`
where the latter collects any `Expr` in `exprs` with `head` being `:(=)`
and the former collects the rest.
"""
function _args_kwargs(exprs)
    args = :(Any[])
    kwargs = :(Dict{Symbol,Any}())
    for expr in exprs
        if expr isa Expr && expr.head==:(=)
            key = Expr(:quote, expr.args[1])
            push!(kwargs.args, Expr(:call, :(=>), key, expr.args[2]))
        else
            push!(args.args, expr)
        end
    end
    return args, kwargs
end

function _spec_walker1(x, parsers, formatters, args_set)
    @capture(x, spec_(formatter_(parser_(rawargs__))...)(;o__)) || return x
    # spec may be a GlobalRef
    name = spec isa Symbol ? spec : spec.name
    name == :StatsSpec || return x
    push!(parsers, parser)
    push!(formatters, formatter)
    return :(push!($args_set, $parser($(rawargs...))))
end

function _spec_walker2(x, parsers, formatters, args_set)
    @capture(x, spec_(formatter_(parser_(rawargs__))...)) || return x
    # spec may be a GlobalRef
    name = spec isa Symbol ? spec : spec.name
    name == :StatsSpec || return x
    push!(parsers, parser)
    push!(formatters, formatter)
    return :(push!($args_set, $parser($(rawargs...))))
end

"""
    @specset [option option=val ...] default_args... begin ... end
    @specset [option option=val ...] default_args... for v in (...) ... end
    @specset [option option=val ...] default_args... for v in (...), w in (...) ... end

Construct a vector of [`StatsSpec`](@ref)s with shared default values for arguments
and then conduct the procedures by calling [`proceed`](@ref).
The vector of [`StatsSpec`](@ref)s and a vector of result objects are returned,
unless alternative behavior is specified with the `option`s.

# Arguments
- `[option option=val ...]`: optional settings for `@specset` including keyword arguments for [`proceed`](@ref).
- `default_args...`: optional default values for arguments shared by all [`StatsSpec`](@ref)s.
- `code block`: a `begin/end` block or a `for` loop containing arguments for constructing [`StatsSpec`](@ref)s.

# Notes
`@specset` transforms `Expr`s that construct [`StatsSpec`](@ref)
to collect the sets of arguments from the code block
and infers how the arguments entered by users need to be processed
based on the names of functions called within [`StatsSpec`](@ref).
For end users, `Macro`s that generate `Expr`s for these function calls should be provided.

Optional default arguments are merged
with the arguments provided for each individual specification
and supersede any default value associated with each [`StatsStep`](@ref)
via [`default`](@ref).
These default arguments should be specified in the same pattern as
how arguments are specified for each specification inside the code block,
as `@specset` processes these arguments by calling
the same functions found in the code block.

Options for the behavior of `@specset` can be provided in a bracket `[...]`
as the first argument with each option separated by white space.
For options that take a Boolean value,
specifying the name of the option is enough for setting the value to be true.
    
The following options are available for altering the behavior of `@specset`:
- `noproceed::Bool=false`: do not call [`proceed`](@ref) and only return the vector of [`StatsSpec`](@ref)s.
- `verbose::Union{Bool,IO}=false`: print the name of each step to `stdout` or the specified `IO` stream when it is called.
- `keep=nothing`: names (of type `Symbol`) of additional objects to be returned.
- `keepall::Bool=false`: return all objects generated by procedures along with arguments from the [`StatsSpec`](@ref)s.
- `pause::Int=0`: break the iteration over [`StatsStep`](@ref)s after finishing the specified number of steps (for debugging).
"""
macro specset(args...)
    nargs = length(args)
    nargs == 0 && throw(ArgumentError("no argument is found for @specset"))
    options = :(Dict{Symbol,Any}())
    noproceed = false
    default_args = nothing
    if nargs > 1
        if isexpr(args[1], :vect, :hcat, :vcat)
            noproceed = _parse!(options, args[1].args)
            nargs > 2 && (default_args = _args_kwargs(args[2:end-1]))
        else
            default_args = _args_kwargs(args[1:end-1])
        end
    end
    specs = macroexpand(__module__, args[end])
    isexpr(specs, :block, :for) ||
        throw(ArgumentError("last argument to @specset must be begin/end block or for loop"))

    # parser and formatter may be GlobalRef
    parsers, formatters, args_set = [], [], Dict{Symbol,Any}[]
    walked = postwalk(x->_spec_walker1(x, parsers, formatters, args_set), specs)
    walked = postwalk(x->_spec_walker2(x, parsers, formatters, args_set), walked)
    nparser = length(unique!(parsers))
    nparser == 1 ||
        throw(ArgumentError("found $nparser parsers from arguments while expecting one"))
    nformatter = length(unique!(formatters))
    nformatter == 1 ||
        throw(ArgumentError("found $nformatter formatters from arguments while expecting one"))

    parser, formatter = parsers[1], formatters[1]
    if default_args === nothing
        defaults = :(Dict{Symbol,Any}())
    else
        defaults = esc(:($parser($(default_args[1]), $(default_args[2]))))
    end

    blk = quote
        $(esc(walked))
        local nspec = length($args_set)
        local sps = Vector{StatsSpec}(undef, nspec)
        for i in 1:nspec
            sps[i] = StatsSpec($(esc(formatter))(merge($defaults, $args_set[i]))...)
        end
        # The same args_set will be reused if exectuted again without recalling the macro
        empty!($args_set)
    end

    if noproceed
        return quote
            $blk
            sps
        end
    else
        return quote
            $blk
            sps, proceed(sps; $options...)
        end
    end
end

end # module
