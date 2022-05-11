using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Test, Zygote
using ChainRulesTestUtils
using DifferentiableStateSpaceModels.Examples

@testset "SGU First Order" begin
    m = @include_example_module(Examples.sgu)
    p_d = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
           ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442)
    p_f = nothing

    c = SolverCache(m, Val(1), p_d)
    settings = PerturbationSolverSettings(; print_level = 0)
    sol = generate_perturbation(m, p_d, p_f; cache = c, settings)
    generate_perturbation_derivatives!(m, p_d, p_f, c)

    # TODO check that these are actually the correct results, currently only a regression test
    # Also, does not check the derivative details in the cache
    @test sol.retcode == :Success
    @test sol.A ≈ [0.42 1.7519543732312002e-17 -4.1734466435163804e-17 9.83532499523268e-19;
                   0.9006855710385308 0.9738235409647132 -1.685401161087487 0.027873826122536263;
                   0.6721065922179434 -0.006588281444814813 0.5003113265680597 -0.0001885768865857994;
                   0.016707717342759755 0.01806442668489436 -0.031264191538169755 0.000517059474573019]
    @test sol.C ≈
          [0.9006855710385274 0.9738235409647122 -1.6854011610874828 0.02787382612253622;
           1.2604305985051094 -0.03920634002915657 0.50643088763379 -0.001122206086526858;
           1.2903225806451621 -1.1817797818825905e-16 0.4129032258064517 -1.1862432475956954e-18;
           1.8774193548387106 4.1726226271227295e-17 0.600774193548387 6.519722854879265e-19;
           6.72106592217944 -0.06588281444815657 -3.996886734319368 -0.0018857688658582676;
           0.6721065922179433 -0.006588281444815545 0.5003113265680625 -0.00018857688658582318;
           -0.4905618464963626 0.207427359629203 -0.7024101653193269 0.005937209270617432;
           -0.643512113444817 0.04451868414563279 1.1217837540922557 0.0012742617208145878;
           -0.6059154121590582 0.017609608142098897 1.1338146985036983 -0.018751472641015628;
           0.0006683086937103991 0.0007225770673958052 -0.0012505676615268426 2.06823789829215e-5;
           0.01670771734276043 0.018064426684894292 -0.03126419153816984 0.0005170594745730165;
           1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    @test sol.Q == I
    @test sol.y ≈ [0.7442, 0.11060245636938473, 0.007390615600776165, 0.3964158265111139,
                   -1.0794906932984638, 1.2230943996955819, 1.7243861964374876, -3.5143212879378156,
                   0.0, 0.0, -3.2188758248681983]
    @test sol.x ≈ [0.0, 0.7442, 1.2230943996955819, -3.2188758248681983]
    @test sol.B ≈ [-0.0129; 0.0; 0.0; 0.0]
    @test sol.D === nothing
    @test sol.g_x ≈
          [0.9006855710385274 0.9738235409647122 -1.6854011610874828 0.02787382612253622;
           1.2604305985051094 -0.03920634002915657 0.50643088763379 -0.001122206086526858;
           1.2903225806451621 -1.1817797818825905e-16 0.4129032258064517 -1.1862432475956954e-18;
           1.8774193548387106 4.1726226271227295e-17 0.600774193548387 6.519722854879265e-19;
           6.72106592217944 -0.06588281444815657 -3.996886734319368 -0.0018857688658582676;
           0.6721065922179433 -0.006588281444815545 0.5003113265680625 -0.00018857688658582318;
           -0.4905618464963626 0.207427359629203 -0.7024101653193269 0.005937209270617432;
           -0.643512113444817 0.04451868414563279 1.1217837540922557 0.0012742617208145878;
           -0.6059154121590582 0.017609608142098897 1.1338146985036983 -0.018751472641015628;
           0.0006683086937103991 0.0007225770673958052 -0.0012505676615268426 2.06823789829215e-5;
           0.01670771734276043 0.018064426684894292 -0.03126419153816984 0.0005170594745730165]
    @test sol.η ≈ [-1; 0; 0; 0]
    @test sol.Γ ≈ [0.0129]
end

# No D
function test_first_order(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.B) + sum(sol.C) +
           sum(sol.x_ergodic.Σ.mat)
end
# Gradients.  Can't put in a testset until #117 fixed
#@testset "SGU 1st order Gradients" begin
const m_sgu = @include_example_module(Examples.sgu)
p_d = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
       ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442)
p_f = nothing

test_first_order(p_d, p_f, m_sgu)
gradient((args...) -> test_first_order(args..., p_f, m_sgu), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order(args..., p_f, m_sgu), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)  # note the rtol is not the default, but it is good enough

function test_second_order_no_D(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f, Val(2))
    return sum(sol.y) + sum(sol.x) + sum(sol.A_0) + +sum(sol.A_1) + sum(sol.A_2) +
           sum(sol.B) + sum(sol.C_0) + sum(sol.C_1) + sum(sol.C_2) +
           sum(sol.g_xx) + sum(sol.g_σσ) + sum(sol.g_x)
end
test_second_order_no_D(p_d, p_f, m_sgu)
gradient((args...) -> test_second_order_no_D(args..., p_f, m_sgu), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_second_order_no_D(args..., p_f, m_sgu), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-8)  # note the rtol is not the default, but it is good enough
# end

@testset "SGU Second Order" begin
    m = @include_example_module(Examples.sgu)
    p_d = (γ = 2.0, ω = 1.455, ρ = 0.42, σe = 0.0129, δ = 0.1, ψ = 0.000742, α = 0.32,
           ϕ = 0.028, β = 1.0 / (1.0 + 0.04), r_w = 0.04, d_bar = 0.7442)
    p_f = nothing

    c = SolverCache(m, Val(2), p_d)
    sol = generate_perturbation(m, p_d, p_f, Val(2); cache = c)
    generate_perturbation_derivatives!(m, p_d, p_f, c)  # Solves and fills the cache

    # TODO check that these are actually the correct results, currently only a regression test
    # Also, does not check the derivative details in the cache
    @test sol.retcode == :Success
    @test sol.Q == I
    @test sol.y ≈ [0.7442, 0.11060245636938473, 0.007390615600776165, 0.3964158265111139,
                   -1.0794906932984638, 1.2230943996955819, 1.7243861964374876, -3.5143212879378156,
                   0.0, 0.0, -3.2188758248681983]
    @test sol.x ≈ [0.0, 0.7442, 1.2230943996955819, -3.2188758248681983]
    @test sol.B ≈ [-0.0129; 0.0; 0.0; 0.0]
    @test sol.D === nothing
    @test sol.g_x ≈
          [0.9006855710385274 0.9738235409647122 -1.6854011610874828 0.02787382612253622;
           1.2604305985051094 -0.03920634002915657 0.50643088763379 -0.001122206086526858;
           1.2903225806451621 -1.1817797818825905e-16 0.4129032258064517 -1.1862432475956954e-18;
           1.8774193548387106 4.1726226271227295e-17 0.600774193548387 6.519722854879265e-19;
           6.72106592217944 -0.06588281444815657 -3.996886734319368 -0.0018857688658582676;
           0.6721065922179433 -0.006588281444815545 0.5003113265680625 -0.00018857688658582318;
           -0.4905618464963626 0.207427359629203 -0.7024101653193269 0.005937209270617432;
           -0.643512113444817 0.04451868414563279 1.1217837540922557 0.0012742617208145878;
           -0.6059154121590582 0.017609608142098897 1.1338146985036983 -0.018751472641015628;
           0.0006683086937103991 0.0007225770673958052 -0.0012505676615268426 2.06823789829215e-5;
           0.01670771734276043 0.018064426684894292 -0.03126419153816984 0.0005170594745730165]
    @test sol.η ≈ [-1; 0; 0; 0]
    @test sol.Γ ≈ [0.0129]

    @test sol.A_0 ≈
          [-3.9402619560802135e-20, 0.0006391584993811839, 0.00015405791020935275,
           1.185639016352093e-5]
    @test sol.A_1 ≈
          [0.42 1.7519543732312002e-17 -4.1734466435163804e-17 9.83532499523268e-19;
           0.9006855710385308 0.9738235409647132 -1.685401161087487 0.027873826122536263;
           0.6721065922179434 -0.006588281444814813 0.5003113265680597 -0.0001885768865857994;
           0.016707717342759755 0.01806442668489436 -0.031264191538169755 0.000517059474573019]
    @test sol.A_2[:, :, 1] ≈ [0.0 0.0 0.0 0.0;
                              -0.6058913712160838 -0.0031020337841137775 -0.40284701288926217 -8.878975162048575e-5;
                              -0.03550021179812225 0.0004999494367584898 -0.14268843831504846 1.4310091185952039e-5;
                              -0.0038546588776544625 0.007926733837883587 -0.021291238555819474 0.0002268875123904971]
    @test sol.A_2[:, :, 2] ≈ [0.0 0.0 0.0 0.0;
                              -0.0031020337841140833 -0.012114331571761804 0.026596726056904527 0.01838062634333291;
                              0.0004999494367584289 -0.0030632500289434696 0.007620095670675804 -0.00021437736130562276;
                              0.007926733837883589 0.008407899373508049 -0.014447148738243992 0.0005880527713927997] # add more slices
    @test sol.C_0 ≈ [0.0006391584993811839, 0.0001036018842634584, -6.50134529592093e-21,
                     -4.420914801226232e-21, 0.0015405791020935281, 0.0001540579102093528,
                     -0.0005481221988433036, -0.00042997911595382926, -0.0004299791159538292,
                     4.7425560654083843e-7, 1.185639016352093e-5, 0.0, 0.0, 0.0, 0.0]
    @test sol.C_1 ≈
          [0.9006855710385274 0.9738235409647122 -1.6854011610874828 0.02787382612253622;
           1.2604305985051094 -0.03920634002915657 0.50643088763379 -0.001122206086526858;
           1.2903225806451621 -1.1817797818825905e-16 0.4129032258064517 -1.1862432475956954e-18;
           1.8774193548387106 4.1726226271227295e-17 0.600774193548387 6.519722854879265e-19;
           6.72106592217944 -0.06588281444815657 -3.996886734319368 -0.0018857688658582676;
           0.6721065922179433 -0.006588281444815545 0.5003113265680625 -0.00018857688658582318;
           -0.4905618464963626 0.207427359629203 -0.7024101653193269 0.005937209270617432;
           -0.643512113444817 0.04451868414563279 1.1217837540922557 0.0012742617208145878;
           -0.6059154121590582 0.017609608142098897 1.1338146985036983 -0.018751472641015628;
           0.0006683086937103991 0.0007225770673958052 -0.0012505676615268426 2.06823789829215e-5;
           0.01670771734276043 0.018064426684894292 -0.03126419153816984 0.0005170594745730165;
           1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]
    @test sol.C_2[:, :, 1] ≈
          [-0.6058913712160838 -0.003102033784113778 -0.40284701288926217 -8.878975162048576e-5;
           0.3641856225490242 0.02778598688564625 0.041906068536917314 0.0007953204400113982;
           -9.250270498653839e-15 -1.6521083580989144e-16 2.6335403233441523e-15 -2.259643246988632e-18;
           -6.2196715024389655e-15 -2.032077242658297e-16 2.0907564007853216e-15 -2.070244672048874e-18;
           -20.682729326608023 0.20426072693774494 13.686097933010494 0.005846570499509324;
           -0.035500211798123044 0.0004999494367584962 -0.14268843831504868 1.431009118595195e-5;
           -0.2699169932252056 -0.04172126179447118 0.03163930042536679 -0.0011941908856709196;
           1.5804491334669595 -0.039703298061926455 -0.6000140721981689 -0.0011364305545266838;
           1.545156746130972 -0.014443488006995794 -0.6113076361456838 0.017661920088352694;
           -0.0001486033987300181 0.00032310566021555423 -0.0008620966077361734 9.248278166631295e-6;
           -0.0038546588776544816 0.00792673383788359 -0.02129123855581945 0.00022688751239049705;
           0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0] # TODO, add more slices
    @test sol.g_xx[:, :, 1] ≈
          [-1.2117827424321677 -0.006204067568227556 -0.8056940257785243 -0.00017757950324097152;
           0.7283712450980484 0.0555719737712925 0.08381213707383463 0.0015906408800227963;
           -1.8500540997307678e-14 -3.304216716197829e-16 5.2670806466883046e-15 -4.519286493977264e-18;
           -1.2439343004877931e-14 -4.064154485316594e-16 4.181512801570643e-15 -4.140489344097748e-18;
           -41.365458653216045 0.4085214538754899 27.372195866020988 0.011693140999018648;
           -0.07100042359624609 0.0009998988735169924 -0.28537687663009736 2.86201823719039e-5;
           -0.5398339864504113 -0.08344252358894236 0.06327860085073359 -0.002388381771341839;
           3.160898266933919 -0.07940659612385291 -1.2000281443963379 -0.0022728611090533676;
           3.090313492261944 -0.028886976013991587 -1.2226152722913677 0.03532384017670539;
           -0.0002972067974600362 0.0006462113204311085 -0.0017241932154723467 1.849655633326259e-5;
           -0.007709317755308963 0.01585346767576718 -0.0425824771116389 0.0004537750247809941] # TODO add more slices
    @test sol.g_σσ ≈ [0.0012783169987623678, 0.0002072037685269168, -1.300269059184186e-20,
                      -8.841829602452464e-21, 0.0030811582041870562, 0.0003081158204187056,
                      -0.0010962443976866073, -0.0008599582319076585, -0.0008599582319076584,
                      9.485112130816769e-7, 2.371278032704186e-5]
end