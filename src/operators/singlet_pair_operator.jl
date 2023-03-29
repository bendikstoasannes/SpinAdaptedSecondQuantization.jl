export Qdag, Q

"""
    SingletPairOperator

Pair singlet creation and annihilation operators
"""

struct SingletPairOperator <: Operator
    p::Int
    q::Int
    dagger::Bool
end

function Base.print(io::IO, constraints::Constraints, Q::SingletPairOperator)
    dagger = Q.dagger ? '†' : ""
    print(io, 'Q', dagger, '_')
    print_mo_index(io, constraints, Q.p, Q.q)
end

function exchange_indices(Q::SingletPairOperator, mapping)
    SingletPairOperator(
        exchange_index(Q.p, mapping),
        exchange_index(Q.q, mapping),
        Q.dagger
    )
end

function get_all_indices(Q::SingletPairOperator)
    (Q.p, Q.q)
end

function Base.:(==)(Q1::SingletPairOperator, Q2::SingletPairOperator)
    (Q1.dagger == Q2.dagger) && (((Q1.p, Q1.q) == (Q2.p, Q2.q)) || (((Q1.p, Q1.q) == (Q2.q, Q2.p))))
end

function Base.isless(Q1::SingletPairOperator, Q2::SingletPairOperator)
    (Q1.dagger, Q1.p, Q1.q) < (Q2.dagger, Q2.p, Q2.q)
end

"""
Pair creation operator
Qdag(p, q) = a†(pα) * a†(qβ) - a†(pβ) * a†(qα)

Pair annihilation operator
Q(p, q) = a(pβ) * a(qα) - a(pα) * a(pβ)
"""

Qdag(p, q) = Expression(SingletPairOperator(p, q, true))
Q(p, q) = Expression(SingletPairOperator(p, q, false))

function convert_to_elementary_operators(Q::SingletPairOperator)
    if Q.dagger
        Expression(
            [(fermiondag(Q.p, α) * fermiondag(Q.q, β))[1], - (fermiondag(Q.p, β) * fermiondag(Q.q, α))[1]]
        )

    else
        Expression(
            [(fermion(Q.p, β) * fermion(Q.q, α))[1], - (fermion(Q.p, α) * fermion(Q.q, β))[1]]
        )
    end
end

function act_on_ket(Q::SingletPairOperator)
    Expression(Q) * if Q.dagger
        virtual(Q.p) * virtual(Q.q)
    else
        occupied(Q.p) * occupied(Q.q)
    end
end
