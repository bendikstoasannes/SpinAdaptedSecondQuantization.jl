using Test

using SpinAdaptedSecondQuantization

@testset "indices" begin
    println()

    p = gen("p")
    i = occ("i")
    a = vir("a")
    @show p i a

    @test (space(p) <: GeneralOrbital) == true
    @test (space(p) <: OccupiedOrbital) == false
    @test (space(i) <: GeneralOrbital) == true
    @test (space(i) <: OccupiedOrbital) == true

    @test (space(i) == GeneralOrbital) == false
    @test (space(i) == OccupiedOrbital) == true

    @test isdisjoint(p, i) == false
    @test isdisjoint(a, i) == true

    @test isocc(i) == true
    @test isvir(i) == false
    @test isocc(a) == false
    @test isvir(a) == true

    println()
end

@testset "kronecker delta" begin
    println()

    p = gen("p")
    i = occ("i")
    a = vir("a")

    dpa = SpinAdaptedSecondQuantization.KroneckerDelta(p, a)
    dia = SpinAdaptedSecondQuantization.KroneckerDelta(i, a)
    dip = SpinAdaptedSecondQuantization.KroneckerDelta(i, p)

    @show dpa dia dip

    @test dpa != 0
    @test dia == 0
    @test dip != 0

    println()
end
