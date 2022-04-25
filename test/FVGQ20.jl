using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Test, Zygote, Statistics
using ChainRulesTestUtils

# need the m as const, so can't put in a testset immediately.
# @testset "FVGQ20 First Order" begin

isdefined(Main, :FVGQ20) || include(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                             "test/generated_models/FVGQ20.jl"))
const m_fvgq = PerturbationModel(Main.FVGQ20)
p_d = (β = 0.998, h = 0.97, ϑ = 1.17, κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36),
       σ_μ = exp(-5.43), σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
p_f = (δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5))

c = SolverCache(m_fvgq, Val(1), p_d)
sol = generate_perturbation(m_fvgq, p_d, p_f; cache = c)
generate_perturbation_derivatives!(m_fvgq, p_d, p_f, c)

@test sol.retcode == :Success
@test sol.A ≈
      [0.9597066320040756 -0.0016717402006548444 -0.0028476359898478745 2.150192831740316e-18 -0.006077464817194368 0.00022438953543368282 -0.002315915210730716 0.0002480013604186878 0.027777206268578707 -0.00383188379414908 -0.209970722216548 -0.9904138263557944 -0.008023324298765148 -0.004801556142159064;
       -0.06116209317503432 0.9081857564695927 -0.07117720022389305 1.4035175039973e-17 -0.1589622523374888 0.005869135737739553 -0.04337176006047146 0.001967204841407287 -0.0020170583816015585 -0.05634485080808087 -0.07734742427383892 -0.2888481452521955 -0.20985817937725523 -0.06009042851246476;
       0.047149679377472 0.05989700727193322 0.5175441967636624 -8.785633176355972e-18 -0.19497474630453365 0.007198773511681691 0.06605149935298901 -0.001731691775730205 0.0014797381625939485 0.05674820197564877 0.0010549545833333184 -0.06512643814483171 -0.2574010161679107 0.05252057306976145;
       0.525972081627061 0.7067647675166118 -0.8250514428090694 -2.0410009966095642e-16 -3.139152046485734 0.1159024182868347 1.0968606775309409 -0.017539660702373366 0.013562717425228943 0.8004565581415306 -0.024947971816113274 -0.8294252329988305 -4.144233763530546 0.6080461766818812;
       0.049282152560429886 0.05468765564912119 0.14582797660460228 5.256507747973834e-17 0.6831385774943626 -0.02522254799519321 0.03467394774591567 -0.00174994781502999 0.0014404882645437057 0.020696434081023268 0.014595652596892201 -0.001347799955404436 0.9018632790315961 0.05276005651146111;
       1.2382273640942494 1.294467918405045 -0.30675891690633167 -2.1373930578377843e-17 -1.0073328208733665 0.037192308059684084 0.5258480374286546 -0.04336430305608645 0.03512566747437948 0.13191493606079913 0.0883642693379878 -1.3357085443757348 -1.3298567974907014 1.3041462486643989;
       0.009599564765915364 0.012194890997879147 -0.0228957392867415 -1.3387600064770917e-17 -0.03969640365736839 0.0014656534991444518 0.8641998463024083 -0.0003525684071545177 0.00030127123908204754 0.011553801582813085 0.0002147862930288743 -0.013259590928356465 -0.05240624662057069 0.01069306619655468;
       -0.0611620931749605 0.9081857564697143 -0.0711772002238785 -8.141584778182808e-17 -0.15896225233760544 0.005869135737619167 -0.043371760060328696 0.9693452541489845 -0.002017058381581333 -0.05634485080799174 -10.48861951082978 -10.706368869464443 -0.2098581793773417 -0.060090428512097956;
       -1.2274387057378921e-18 1.129716330863846e-18 -2.6452513161575937e-18 4.303366036318858e-20 1.932605834302807e-18 -4.9112240214915196e-20 -5.5502118613030306e-18 3.2874122628061265e-19 0.12 -4.230985837363517e-18 3.0274528284264972e-18 3.2292323574232476e-18 1.7423640484797458e-18 -1.6502917887199003e-17;
       6.5743501365246386e-18 1.69567809349733e-17 1.8453482603434553e-18 1.4104454171393302e-17 1.0721975244983874e-16 -5.853396590979275e-18 2.3016001541456748e-17 -4.972737819827135e-19 -3.546634649666224e-18 0.9300000000000004 3.521544976745966e-16 4.400525620152486e-16 2.113466586628283e-16 2.707831559404897e-16;
       2.2568972501382555e-17 1.1941207210767652e-16 -7.886281870178793e-17 -8.262052872093175e-31 -1.7082021990020085e-16 6.306950502338531e-18 1.0385973714711463e-16 -9.933762539483766e-18 7.997775589584522e-19 6.113195184305982e-17 8.986093821921679e-17 9.333325263923877e-17 -6.301399549422952e-17 -6.637664205478931e-17;
       -2.3295865035336392e-17 -3.648847602877501e-18 2.5686152263757588e-17 -9.1576245100504e-31 5.249861165160181e-17 -1.9383311031746743e-18 -2.9498949743095986e-17 2.4517184677316216e-19 -7.10836567443952e-19 -2.014313976204385e-17 3.4655363692211897e-18 1.5444545238035702e-17 1.0724413225403257e-17 8.242097119011588e-17;
       -1.8469211863662103e-16 -1.6043964703540862e-16 6.123339946925812e-18 1.6457798003133015e-30 1.012214486991828e-16 -3.737254681960731e-18 -4.189986148448415e-17 6.012143234755371e-18 -5.662286392385015e-18 -3.254482025735823e-17 -1.1428499435437241e-17 1.891055157616275e-16 4.069714234131922e-17 -6.806721471364512e-17;
       -7.065490002090046e-18 -4.1123352461888506e-18 -3.786545338319949e-17 8.380886195309565e-18 -2.7207396233759008e-17 -2.88508375354103e-19 -9.860474762111539e-18 3.1037019916073786e-19 -6.578099119381701e-19 -1.0956403414738988e-17 -2.7651598835201586e-17 3.323414050181414e-17 -1.8308120469745517e-16 0.9499999999999996]
@test sol.C ≈
      [0.046682850885188414 0.05930396759586929 0.5124199967950123 0.0 -0.19304430327244906 0.007127498526411165 0.06539752411292642 -0.0017145463126460162 0.0014650872902314323 0.056186338595056756 0.0010445094854684119 -0.0644816219419623 -0.2548524912562035 0.052000567405129226;
       0.04848049555975208 0.05379806906178239 0.1434558395931211 0.0 0.6720261809473577 -0.024812260881409483 0.03410991773797231 -0.0017214819741272426 0.0014170563030341846 0.020359771812883194 0.014358229788893229 -0.0013258757233242705 0.8871929577850087 0.05190182556133141;
       0.4733058027331458 0.6359954782844232 -0.7424379526200766 -0.8998686798130806 -2.8248246078053247 0.10429695613070526 0.9870305698304098 -0.015783391320695707 0.012204664624805324 0.7203057862292772 -0.02244989846562421 -0.7463737894434745 -3.729266165627129 0.5471617102883314;
       0.7924415036295586 0.8284343678742737 -0.19631975872429247 0.0 -0.6446734730418766 -0.616178253872422 0.33653254769120866 -0.027752313116301518 0.0224797460923251 0.08442300123162487 0.056551418981011355 -0.854827568836762 -0.851082564294199 0.8346283115774529;
       0.7095449970937021 0.7149458978760908 -0.06799364123668006 0.0 -0.10822645419616717 0.003995890475364564 1.0031833138000825 -0.025718053394227534 0.020408148007625304 -0.051400061887996056 0.06913294700293195 -0.7124691120650971 -0.14287798709487276 0.7379111193434185;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.9966057734548978 0.0 0.0 0.0]
@test sol.Q ≈
      [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
@test sol.y ≈ [0.8154780241332105, 0.2783051143056207, 1.01, 1.1112732584577985,
               1.0165356601604938, 1.562547340626902, 1.0022208173087115, 1.0186852828361062,
               1.0, 1.0, 1.2953143777555123, 1.0044580087526465, 0.03489877588693624, 1.0,
               0.8982498823496541, 12.044079264159267, 13.382310293510294, 0.004448101265822784,
               0.009950330853168092, 0.01640043478966387, 0.0, 0.0, 0.0, 0.0034]
@test sol.x ≈ [0.8154780241332105, 0.2783051143056207, 1.01, 1.1112732584577985,
               1.0165356601604938, 1.562547340626902, 1.0022208173087115, 8.53122233338444, 1.0,
               1.0, 1.0034057865562385, 1.0028039236612292, 1.0, 0.4687642021880705]
@test sol.B ≈
      [0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0;
       0.2209099779593782 0.0 0.0 0.0 0.0 0.0; 0.0 0.09442022319630236 0.0 0.0 0.0 0.0;
       0.0 0.0 0.004383095802668776 0.0 0.0 0.0;
       0.0 0.0 0.0 0.018873433135151486 0.0 0.0;
       0.0 0.0 0.0 0.0 0.002879899158088243 0.0;
       0.0 0.0 0.0 0.0 0.0 0.049787068367863944]
@test sol.η ≈
      [0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 1.0 0.0 0.0 0.0 0.0 0.0;
       0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 1.0 0.0 0.0 0.0; 0.0 0.0 0.0 1.0 0.0 0.0;
       0.0 0.0 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0]
@test sol.Γ ≈
      [0.2209099779593782 0.0 0.0 0.0 0.0 0.0; 0.0 0.09442022319630236 0.0 0.0 0.0 0.0;
       0.0 0.0 0.004383095802668776 0.0 0.0 0.0;
       0.0 0.0 0.0 0.018873433135151486 0.0 0.0;
       0.0 0.0 0.0 0.0 0.002879899158088243 0.0;
       0.0 0.0 0.0 0.0 0.0 0.049787068367863944]
# end

function test_first_order_smaller(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.A)
end

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order_smaller(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)

# Checked following with no errors: x, y, Γ, B
# Errors on: A, C, g_x

function test_first_order_temp(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.g_x)
end

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order_temp(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)

# The bigger test, not required since fails for smaller           
function test_first_order(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.B) + sum(sol.C) +
           sum(cov(sol.D)) + sum(sol.x_ergodic.Σ.mat)
end
test_first_order(p_d, p_f, m_fvgq)
gradient((args...) -> test_first_order(args..., p_f, m_fvgq), p_d)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)
###############
# checking indivdual functions with simpler setup
# Only one parameter necessary for failures. β or h for example.  σ_m doesn't fail since only in Γ
p_d = (; h = 0.97)
p_f = (β = 0.998, δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5),
       ϑ = 1.17,
       κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36), σ_μ = exp(-5.43),
       σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
c = SolverCache(m_fvgq, Val(1), p_d)
sol = generate_perturbation(m_fvgq, p_d, p_f; cache = c)
generate_perturbation_derivatives!(m_fvgq, p_d, p_f, c)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order_smaller(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)
###############
# @testset "FVGQ20 Second Order" begin
isdefined(Main, :FVGQ20) || include(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                             "test/generated_models/FVGQ20.jl"))
const m_fvgq_2 = PerturbationModel(Main.FVGQ20)
p_d = (β = 0.998, h = 0.97, ϑ = 1.17, κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36),
       σ_μ = exp(-5.43), σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
p_f = (δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5))

c = SolverCache(m_fvgq_2, Val(2), p_d)
sol = generate_perturbation(m_fvgq_2, p_d, p_f, Val(2); cache = c)
generate_perturbation_derivatives!(m_fvgq_2, p_d, p_f, c)

# Only test whether it can run
@test sol.retcode == :Success
# end

# some of these are instead covered in the sequential repo.

# @testset "FVGQ20 Kalman likelhood derivative in 1st order" begin
#     path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
#     file_prefix = "FVGQ20"
#     A = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/data/FVGQ20_A.csv"))
#     B = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/data/FVGQ20_B.csv"))
#     C = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/data/FVGQ20_C.csv"))
#     D_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                              "test/data/FVGQ20_D.csv"))
#     D = MvNormal(Diagonal(abs2.(D_raw)))
#     observables_raw = Matrix(DataFrame(CSV.File(joinpath(path, "FVGQ20_observables.csv");
#                                                 header = false)))
#     noise_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                  "test/data/FVGQ20_noise.csv"))
#     observables = [observables_raw[i, :] for i in 1:size(observables_raw, 1)]
#     noise = [noise_raw[i, :] for i in 1:size(noise_raw, 1)]
#     u0_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                               "test/data/FVGQ20_ergodic.csv"))
#     u0 = MvNormal(Symmetric(u0_raw))

#     minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables)

#     res = gradient(minimal_likelihood_test_kalman_first, A, B, C, D, u0, noise, observables)

#     # Some tests
#     @test finite_difference_gradient(A -> minimal_likelihood_test_kalman_first(A, B, C, D,
#                                                                                u0, noise,
#                                                                                observables),
#                                      A) ≈ res[1] rtol = 1e-3
#     @test finite_difference_gradient(B -> minimal_likelihood_test_kalman_first(A, B, C, D,
#                                                                                u0, noise,
#                                                                                observables),
#                                      B) ≈ res[2] rtol = 1e-3
#     @test finite_difference_gradient(C -> minimal_likelihood_test_kalman_first(A, B, C, D,
#                                                                                u0, noise,
#                                                                                observables),
#                                      C) ≈ res[3] rtol = 1e-3

#     # missing test for the D = MvNormal.  No finite_difference_gradient support yet.

#     # @test finite_difference_gradient(u0 -> minimal_likelihood_test_kalman_first(A, B, C, D, u0, noise, observables), u0) ≈ res[5] rtol=1E-7

#     observables_grad = finite_difference_gradient(observables_mat -> minimal_likelihood_test_kalman_first(A,
#                                                                                                           B,
#                                                                                                           C,
#                                                                                                           D,
#                                                                                                           u0,
#                                                                                                           noise,
#                                                                                                           [observables_mat[i,
#                                                                                                                            :]
#                                                                                                            for i in
#                                                                                                                1:size(observables_mat,
#                                                                                                                       1)]),
#                                                   observables_raw)
#     @test [observables_grad[i, :] for i in 1:size(observables_raw, 1)] ≈ res[7] rtol = 1e-5
# end

# function minimal_likelihood_test_joint_first(A, B, C, D, u0, noise, observables)
#     problem = LinearStateSpaceProblem(A, B, C, u0, (0, length(observables)); noise = noise,
#                                       obs_noise = D, observables = observables)
#     return solve(problem, NoiseConditionalFilter()).loglikelihood
# end

# @testset "FVGQ20 joint likelhood derivative in 1st order" begin
#     path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
#     file_prefix = "FVGQ20"
#     A = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/data/FVGQ20_A.csv"))
#     B = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/data/FVGQ20_B.csv"))
#     C = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels), "test/data/FVGQ20_C.csv"))
#     D_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                              "test/data/FVGQ20_D.csv"))
#     D = MvNormal(Diagonal(map(abs2, vec(D_raw))))
#     observables_raw = Matrix(DataFrame(CSV.File(joinpath(path, "FVGQ20_observables.csv");
#                                                 header = false)))
#     noise_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                  "test/data/FVGQ20_noise.csv"))
#     observables = [observables_raw[i, :] for i in 1:size(observables_raw, 1)]
#     noise = [noise_raw[i, :] for i in 1:size(noise_raw, 1)]
#     u0 = zeros(size(A, 1))

#     res = gradient(minimal_likelihood_test_joint_first, A, B, C, D, u0, noise, observables)

#     # Some tests
#     @test finite_difference_gradient(A -> minimal_likelihood_test_joint_first(A, B, C, D,
#                                                                               u0, noise,
#                                                                               observables),
#                                      A) ≈ res[1] rtol = 1e-5
#     @test finite_difference_gradient(B -> minimal_likelihood_test_joint_first(A, B, C, D,
#                                                                               u0, noise,
#                                                                               observables),
#                                      B) ≈ res[2] rtol = 1e-5
#     @test finite_difference_gradient(C -> minimal_likelihood_test_joint_first(A, B, C, D,
#                                                                               u0, noise,
#                                                                               observables),
#                                      C) ≈ res[3] rtol = 1e-5

#     # missing test for the D = MvNormal.  No finite_difference_gradient support yet.

#     @test finite_difference_gradient(u0 -> minimal_likelihood_test_joint_first(A, B, C, D,
#                                                                                u0, noise,
#                                                                                observables),
#                                      u0) ≈ res[5] rtol = 1e-7

#     noise_grad = finite_difference_gradient(noise_mat -> minimal_likelihood_test_joint_first(A,
#                                                                                              B,
#                                                                                              C,
#                                                                                              D,
#                                                                                              u0,
#                                                                                              [noise_mat[i,
#                                                                                                         :]
#                                                                                               for i in
#                                                                                                   1:size(noise_mat,
#                                                                                                          1)],
#                                                                                              observables),
#                                             noise_raw)
#     @test [noise_grad[i, :] for i in 1:size(noise_raw, 1)] ≈ res[6] rtol = 1E-7

#     observables_grad = finite_difference_gradient(observables_mat -> minimal_likelihood_test_joint_first(A,
#                                                                                                          B,
#                                                                                                          C,
#                                                                                                          D,
#                                                                                                          u0,
#                                                                                                          noise,
#                                                                                                          [observables_mat[i,
#                                                                                                                           :]
#                                                                                                           for i in
#                                                                                                               1:size(observables_mat,
#                                                                                                                      1)]),
#                                                   observables_raw)
#     @test [observables_grad[i, :] for i in 1:size(observables_raw, 1)] ≈ res[7] rtol = 1E-7
# end

# @testset "Gradients, generate_perturbation + simulation, 2nd order" begin
#     m = @include_example_module(Examples.rbc_observables)
#     p_f = (ρ=0.2, δ=0.02, σ=0.01, Ω_1=0.1)
#     p_d = (α=0.5, β=0.95)
#     p_d_input = [0.5, 0.95]

#     T = 9
#     ϵ_mat = [0.22, 0.01, 0.14, 0.03, 0.15, 0.21, 0.22, 0.05, 0.18]
#     x0 = zeros(m.n_x)

#     function sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; kwargs...)
#         p_d = (α=p_d_input[1], β=p_d_input[2])
#         sol = generate_perturbation(m, p_d, p_f, Val(2); kwargs...)
#         ϵ = map(i -> ϵ_mat[i:i], 1:T)
#         simul = solve(
#             DifferentiableStateSpaceModels.dssm_evolution,
#             DifferentiableStateSpaceModels.dssm_volatility,
#             [x0; x0],
#             (0, T),
#             sol;
#             h = DifferentiableStateSpaceModels.dssm_observation,
#             noise = ϵ,
#         )
#         return sum(sum(simul.z))
#     end

#     settings = PerturbationSolverSettings()
#     # @inferred sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings)
#     res_zygote = gradient(
#         (p_d_input, ϵ_mat) -> sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
#         p_d_input,
#         ϵ_mat,
#     )
#     res_finite_p = finite_difference_gradient(
#         p_d_input -> sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
#         p_d_input,
#     )
#     res_finite_ϵ = finite_difference_gradient(
#         ϵ_mat -> sum_test_joint_second(p_d_input, ϵ_mat, x0, T, p_f, m; settings),
#         ϵ_mat,
#     )
#     @test res_zygote[2] ≈ res_finite_ϵ
#     @test isapprox(res_zygote[1][1], res_finite_p[1]; rtol = 1e-5)
#     @test isapprox(res_zygote[1][2], res_finite_p[2]; rtol = 1e-5)
# end

# @testset "FVGQ20 joint likelhood derivative in 2nd order" begin
#     path = joinpath(pkgdir(DifferentiableStateSpaceModels), "test", "data")
#     file_prefix = "FVGQ20"
#     A_0_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_A_0.csv"))
#     A_0 = vec(A_0_raw)
#     A_1 = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                            "test/data/$(file_prefix)_A_1.csv"))
#     A_2_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_A_2.csv"))
#     A_2 = reshape(A_2_raw, length(A_0), length(A_0), length(A_0))
#     B = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                          "test/data/$(file_prefix)_B.csv"))
#     C_0_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_C_0.csv"))
#     C_0 = vec(C_0_raw)
#     C_1 = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                            "test/data/$(file_prefix)_C_1.csv"))
#     C_2_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                "test/data/$(file_prefix)_C_2.csv"))
#     C_2 = reshape(C_2_raw, length(C_0), length(A_0), length(A_0))
#     D_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                              "test/data/$(file_prefix)_D.csv"))
#     D = MvNormal(Diagonal(map(abs2, vec(D_raw))))
#     observables_raw = Matrix(DataFrame(CSV.File(joinpath(path,
#                                                          "$(file_prefix)_observables.csv");
#                                                 header = false)))
#     noise_raw = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
#                                  "test/data/$(file_prefix)_noise.csv"))
#     observables = [observables_raw[i, :] for i in 1:size(observables_raw, 1)]
#     noise = [noise_raw[i, :] for i in 1:size(noise_raw, 1)]
#     u0 = zeros(length(A_0))

#     res = gradient(minimal_likelihood_test_joint_second, A_0, A_1, A_2, B, C_0, C_1, C_2, D,
#                    u0, noise, observables)

#     # Some tests
#     @test finite_difference_gradient(A_0 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      A_0) ≈ res[1] rtol = 1E-5
#     @test finite_difference_gradient(A_1 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      A_1) ≈ res[2] rtol = 1E-5
#     # I didn't add the tests for A_2 and C_2 because of their dimensions
#     @test finite_difference_gradient(B -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                A_2, B, C_0,
#                                                                                C_1, C_2, D,
#                                                                                u0, noise,
#                                                                                observables),
#                                      B) ≈ res[4] rtol = 1E-5
#     @test finite_difference_gradient(C_0 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      C_0) ≈ res[5] rtol = 1E-5
#     @test finite_difference_gradient(C_1 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                  A_2, B,
#                                                                                  C_0, C_1,
#                                                                                  C_2, D, u0,
#                                                                                  noise,
#                                                                                  observables),
#                                      C_1) ≈ res[6] rtol = 1E-5
#     # missing test for the D = MvNormal.  No finite_difference_gradient support yet.

#     @test finite_difference_gradient(u0 -> minimal_likelihood_test_joint_second(A_0, A_1,
#                                                                                 A_2, B, C_0,
#                                                                                 C_1, C_2, D,
#                                                                                 u0, noise,
#                                                                                 observables),
#                                      u0) ≈ res[9] rtol = 1E-7

#     noise_grad = finite_difference_gradient(noise_mat -> minimal_likelihood_test_joint_second(A_0,
#                                                                                               A_1,
#                                                                                               A_2,
#                                                                                               B,
#                                                                                               C_0,
#                                                                                               C_1,
#                                                                                               C_2,
#                                                                                               D,
#                                                                                               u0,
#                                                                                               [noise_mat[i,
#                                                                                                          :]
#                                                                                                for i in
#                                                                                                    1:size(noise_mat,
#                                                                                                           1)],
#                                                                                               observables),
#                                             noise_raw)
#     @test [noise_grad[i, :] for i in 1:size(noise_raw, 1)] ≈ res[10] rtol = 1E-7

#     observables_grad = finite_difference_gradient(observables_mat -> minimal_likelihood_test_joint_second(A_0,
#                                                                                                           A_1,
#                                                                                                           A_2,
#                                                                                                           B,
#                                                                                                           C_0,
#                                                                                                           C_1,
#                                                                                                           C_2,
#                                                                                                           D,
#                                                                                                           u0,
#                                                                                                           noise,
#                                                                                                           [observables_mat[i,
#                                                                                                                            :]
#                                                                                                            for i in
#                                                                                                                1:size(observables_mat,
#                                                                                                                       1)]),
#                                                   observables_raw)
#     @test [observables_grad[i, :] for i in 1:size(observables_raw, 1)] ≈ res[11] rtol = 1E-7

#     # inference
#     @inferred solve(A_0, A_1, A_2, B, C_0, C_1, C_2, D, u0, (0, length(observables)),
#                     DifferentiableStateSpaceModels.QTILikelihood(); noise, observables)
# end
