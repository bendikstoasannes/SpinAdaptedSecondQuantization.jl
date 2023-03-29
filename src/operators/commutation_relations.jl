# Implement commutation relations here

# To implement a commutation relation you implement the function
# `reductive_commutator(a, b)` between the two operator types.
# The function should return (∓1, [a, b]±), where the first integer should
# be the sign change of commutation (1 for commutator, -1 for anticommutator)

# When implementing a commutation relation between two different operator types
# only one order is required, as a generic function takes care of the other
# i.e. implement only [a, b]± and not also [b, a]±

function reductive_commutator(
    a::SingletExcitationOperator,
    b::SingletExcitationOperator
)
    p = a.p
    q = a.q
    r = b.p
    s = b.q

    (1, δ(q, r) * E(p, s) - δ(p, s) * E(r, q))
end

function reductive_commutator(a::FermionOperator, b::FermionOperator)
    (-1, δ(a.p, b.p) * (a.spin == b.spin) * (a.dag != b.dag))
end

function reductive_commutator(e::SingletExcitationOperator, a::FermionOperator)
    p = e.p
    q = e.q
    r = a.p

    (1, if a.dag
        δ(q, r) * fermiondag(p, a.spin)
    else
        -δ(p, r) * fermion(q, a.spin)
    end)
end

function reductive_commutator(Q1::SingletPairOperator, a1::FermionOperator)
    p = Q1.p
    q = Q1.q
    r = a1.p
    σ = a1.spin

    if σ == α
        τ = β
    else
        τ = α
    end

    (1, !(Q1.dagger == a1.dag) * if a1.dag
        δ(p, r) * fermion(q, τ) + δ(q, r) * fermion(p, τ)
    else
        -δ(p, r) * fermiondag(q, τ) - δ(q, r) * fermiondag(p, τ)
    end)

end
