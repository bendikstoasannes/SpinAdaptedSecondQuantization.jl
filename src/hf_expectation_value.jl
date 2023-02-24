export hf_expectation_value

function hf_expectation_value(ex::Expression{T}) where {T<:Number}
    terms = Term{T}[]

    for t in ex.terms
        append!(terms, hf_expectation_value(t).terms)
    end

    Expression(terms)
end

function hf_expectation_value(t::Term)
    op_ex = hf_expectation_value(t.operators, t.constraints)

    Expression([fuse(noop_part(t), t2) for t2 in op_ex.terms])
end

function hf_expectation_value(ops::Vector{Operator}, constraints::Constraints)
    if isempty(ops)
        return Expression(1)
    elseif length(ops) == 1
        return hf_expectation_value(only(ops))
    end

    firstop = first(ops)
    rest = ops[2:end]

    if firstop isa SingletExcitationOperator
        if get(constraints, firstop.p, GeneralOrbital) == VirtualOrbital
            return zero(Expression{Int64})
        else
            if get(constraints, firstop.q, GeneralOrbital) == OccupiedOrbital
                hf_expectation_value(firstop) * hf_expectation_value(rest, constraints)
            elseif get(constraints, firstop.q, GeneralOrbital) == VirtualOrbital
                hf_expectation_value(commutator(Expression(firstop), Expression(rest)) * constrain(constraints))
            else
                assume_occ = hf_expectation_value(firstop) *
                             hf_expectation_value(rest, constraints) *
                             constrain(firstop.q => OccupiedOrbital)

                assume_vir = hf_expectation_value(commutator(
                    Expression(firstop) *
                    constrain(firstop.q => VirtualOrbital),
                    Expression(rest)
                ))

                assume_occ + assume_vir
            end
        end
    else
        throw("hf_expectation_value not implemented for $(typeof(firstop))!")
    end
end

function hf_expectation_value(o::SingletExcitationOperator)
    2 * δ(o.p, o.q) * constrain(o.p => OccupiedOrbital)
end
