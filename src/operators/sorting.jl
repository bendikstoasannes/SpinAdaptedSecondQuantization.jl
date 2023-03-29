# Implement ordering of new operator types here:

function Base.isless(::FermionOperator, ::SingletExcitationOperator)
    false
end

function Base.isless(::SingletPairOperator, ::SingletExcitationOperator)
    false
end

function Base.isless(::SingletPairOperator, ::FermionOperator)
    true
end
