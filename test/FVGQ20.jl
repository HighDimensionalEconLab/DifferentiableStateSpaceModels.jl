using DifferentiableStateSpaceModels, Symbolics, LinearAlgebra, Test, Zygote, Statistics,
      DifferenceEquations
using DelimitedFiles
using ChainRulesTestUtils
using FiniteDiff
using FiniteDiff: finite_difference_derivative, finite_difference_gradient,
                  finite_difference_jacobian, finite_difference_hessian

isdefined(Main, :FVGQ20) || include(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                             "test/generated_models/FVGQ20.jl"))
# Simple test to reproduce the issue
const m_fvgq = PerturbationModel(Main.FVGQ20)

function test_first_order_smaller(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
    return sum(sol.A)
end
function test_second_order_smaller(p_d, p_f, m)
    sol = generate_perturbation(m, p_d, p_f, Val(2)) # manually passing in order
    return sum(sol.A_2)
end
p_d = (; β = 0.998)
p_f = (h = 0.97, δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5),
       ϑ = 1.17,
       κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36), σ_μ = exp(-5.43),
       σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order_smaller(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)

# To examine intermediate values set breakpoints/etc. in the following code
c = SolverCache(m_fvgq, Val(1), p_d)
sol = generate_perturbation(m_fvgq, p_d, p_f; cache = c)
generate_perturbation_derivatives!(m_fvgq, p_d, p_f, c)

# # can also test with finite differences
# eps_fd = sqrt(eps())
# @test (test_first_order_smaller(merge(p_d, (; β = p_d.β + eps_fd)), p_f, m_fvgq) -
#        test_first_order_smaller(p_d, p_f, m_fvgq)) / eps_fd ≈
#       gradient((args...) -> test_first_order_smaller(args...,
#                                                      p_f,
#                                                      m_fvgq),
#                p_d)[1].β

#################################
# More general tests 
# need the m as const, so can't put in a testset immediately.
# @testset "FVGQ20 First Order" begin

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

# Add the tests for the derivatives
@test c.A_1_p[1] ≈
      [0.2578397517574386 0.0444286244840111 0.09037436257911871 9.01880518352488e-16 0.17955992282532945 -0.008576000033306115 0.09115455823125601 -0.016492185183436803 -0.05197630515950965 0.15336698752952785 0.47986617418153316 1.9853761629344475 0.24509013646551328 0.21780106313554073;
       -1.690428215283972 -1.727993514596424 -2.503289668248622 -2.000395244139009e-15 -5.800907015837849 0.16326941232901984 -1.3326909159610765 -0.0004651100413596634 -0.051424783984073706 -1.7954719561774803 -2.8371412639943285 -10.625643102398861 -7.447940539609109 -1.5581185249758627;
       -0.32488901434007605 -1.0144702771071206 -0.15097657751092217 2.366573014523068e-16 -0.1863751958485695 -0.05556125475707354 -0.2928678334741597 0.06068345932636875 -0.012963866515159048 -0.1732087404333219 0.03374881180243761 0.05622965396901809 0.011868755387788725 -0.46213699034220845;
       4.805831807403937 -1.515902819067753 -6.953882717298972 2.1163597350584652e-15 -29.71845176082458 0.09190852252258637 9.900475571550848 0.3587800689368677 0.12954644121354308 9.465165622330245 0.04618106844629505 -11.204382366929776 -35.081053340134815 3.903243085674419;
       -0.4406328730358733 -0.7206129692637564 -0.2548248345387716 2.374632942256155e-16 -0.1414713488102187 0.22400497955386403 -0.1576552923163599 0.06670871544360904 -0.015898344521860203 -0.09460911569866375 -0.018306301096526867 0.025270895916070926 -1.0904377254105966 -0.5036391020763764;
       -0.8498121542866701 -2.832753359644709 -4.5912763684965725 -4.087435026869448e-17 -13.808285076340246 0.187215618214435 2.787190738698333 1.31467310903037 -0.10279167480120857 0.23302368548646907 0.19395102479904183 -9.996853809079527 -16.896847231244745 -0.9960889487695421;
       -0.06614664566334728 -0.20654378262712333 -0.030738479096815537 6.405165158580094e-17 -0.03794555523915411 -0.01131214189710873 -0.05962720791392553 0.012355010802197645 -0.0026394129901466407 -0.03526489562138482 0.006871179379173727 0.011448226417917654 0.0024164509184641495 -0.09409001350871062;
       -1.6904282152839754 -1.7279935145964853 -2.5032896682486223 -5.490120073719792e-16 -5.80090701583785 0.1632694123290185 -1.3326909159610796 -0.00046511004136090734 -0.051424783984073845 -1.795471956177483 -383.19788118822277 -391.21466794920997 -7.447940539609161 -1.5581185249758693;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0;
       -1.3440453420237384e-16 -1.8511964174521392e-16 -9.022077324479055e-16 3.4219266044796424e-17 -8.220742782265363e-18 6.583923181015661e-18 -1.1624593461616652e-16 3.936720771652378e-18 1.2242784242979848e-17 -1.650006860974369e-16 5.483485181572879e-16 1.5452981490826457e-15 -4.173683383300853e-16 -5.8288474763594705e-15]
@test c.C_1_p[1] ≈
      [-0.32167229142581294 -1.0044260169377357 -0.1494817599118069 2.9138849139939464e-16 -0.18452989687977803 -0.0550111433238351 -0.2899681519546116 0.0600826329964037 -0.012835511401147397 -0.1714938024092278 0.03341466515093669 0.055672924721801166 0.0117512429581998 -0.4575613765764384;
       -0.38488758828694664 -0.654985107958978 -0.10693635428679672 1.5716350879575626e-16 0.534202847132583 0.1954991774110054 -0.12091248985975504 0.06389865428712704 -0.014219835145783177 -0.07266956654152869 -0.003621514560764386 0.023531289551776947 -0.18372902442882716 -0.44344072879610463;
       0.75701990663689 -6.158003465079388 -0.6613687154065674 6.782864988311531 -5.450259443008548 -0.7034448542323498 1.4692699954134363 0.44182409746240603 0.024580703859349122 3.0880173298987685 0.21077562165810093 -4.456593395318587 -3.4585630576106317 -0.6118889205952927;
       -6.6235107418973005 -8.168693563696895 -1.4321536869786737 1.6661082798335216e-15 -3.8910703488945924 4.84716209691505 -0.7981450612087537 1.0542823129769907 -0.2382503183767247 -0.49856648436411555 -0.30974017645420626 0.16048358779904434 -4.284110256625778 -7.040784173132267;
       -6.6881346123105825 -6.9354292691104025 -1.4700409699517927 3.691593383210926e-16 -3.319322746976298 0.08789405053571213 -0.8342961905240056 0.9810292060332197 -0.241939645190218 -1.0600746873911364 -0.13096312111504452 1.0554836496486608 -4.2389265152927225 -6.8378463416056565;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_first_order_smaller(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-7)

test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> test_second_order_smaller(args..., p_f, m_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-6) # 1e-7 is a little too tight

# # The bigger test, not required since failures often occur for smaller ones.
# function test_first_order(p_d, p_f, m)
#     sol = generate_perturbation(m, p_d, p_f)#, Val(1); cache = c) # manually passing in order
#     return sum(sol.y) + sum(sol.x) + sum(sol.A) + sum(sol.B) + sum(sol.C) +
#            sum(cov(sol.D)) + sum(sol.x_ergodic.Σ.mat)
# end
# test_first_order(p_d, p_f, m_fvgq)
# gradient((args...) -> test_first_order(args..., p_f, m_fvgq), p_d)

# test_rrule(Zygote.ZygoteRuleConfig(),
#            (args...) -> test_first_order(args..., p_f, m_fvgq), p_d;
#            rrule_f = rrule_via_ad,
#            check_inferred = false, rtol = 1e-6)
###############
# checking indivdual functions with simpler setup  Only one parameter typically necessary for failures
p_d = (; β = 0.998)
p_f = (h = 0.97, δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5),
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

#weird behavior using this.  Probably the general issue with inference and the model
# Run this separately from above code until inference fixed

# function make_univariate_at_symbol(f, X, sym, args...)
#     return x -> f(merge(X, [sym => x]), args...)
# end
#FiniteDiff.derivative(make_univariate_at_symbol(test_first_order_smaller, p_d, :β, p_f, m_fvgq), p_d.β)

# Verifying it with FD and not test_rrule

eps_fd = sqrt(eps())
@test (test_first_order_smaller(merge(p_d, (; β = p_d.β + eps_fd)), p_f, m_fvgq) -
       test_first_order_smaller(merge(p_d, (; β = p_d.β - eps_fd)), p_f, m_fvgq)) /
      (2 * eps_fd) ≈
      gradient((args...) -> test_first_order_smaller(args...,
                                                     p_f,
                                                     m_fvgq),
               p_d)[1].β

# Verifying the issue is with the A_p and the g_x_p and not the rrule using it
h_x_p_fd = (generate_perturbation(m_fvgq, merge(p_d, (; β = p_d.β + eps_fd)), p_f).A -
            generate_perturbation(m_fvgq, merge(p_d, (; β = p_d.β - eps_fd)), p_f).A) ./
           (2 * eps_fd)
@test c.h_x_p[1] ≈ h_x_p_fd rtol = 1e-5
g_x_p_fd = (generate_perturbation(m_fvgq, merge(p_d, (; β = p_d.β + eps_fd)), p_f).g_x -
            generate_perturbation(m_fvgq, p_d, p_f).g_x) ./ eps_fd
@test c.g_x_p[1] ≈ g_x_p_fd rtol = 1e-5 # poor approximation, but could just be finite differences

######
# Checking gradient calculations
using FiniteDiff, ForwardDiff
p_d = (; β = 0.998)
p_f = (h = 0.97, δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5),
       ϑ = 1.17,
       κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36), σ_μ = exp(-5.43),
       σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
c = SolverCache(m_fvgq, Val(1), p_d)
sol = generate_perturbation(m_fvgq, p_d, p_f; cache = c)
generate_perturbation_derivatives!(m_fvgq, p_d, p_f, c)
const m = m_fvgq
function H(X, m)
    y_p = X[1:(m.n_y)]
    y = X[(m.n_y + 1):(2 * m.n_y)]
    y_ss = X[(2 * m.n_y + 1):(3 * m.n_y)]
    x_p = X[(3 * m.n_y + 1):(3 * m.n_y + m.n_x)]
    x = X[(3 * m.n_y + m.n_x + 1):(3 * m.n_y + 2 * m.n_x)]
    x_ss = X[(3 * m.n_y + 2 * m.n_x + 1):(3 * m.n_y + 3 * m.n_x)]
    p = X[(3 * m.n_y + 3 * m.n_x + 1):end]

    out = similar(X, m.n_x + m.n_y) # output 
    m.mod.m.H!(out, y_p, y, y_ss, x_p, x, x_ss, p)
    return out
end
p = DifferentiableStateSpaceModels.order_vector_by_symbols(merge(p_d, p_f),
                                                           m.mod.m.p_symbols)
X = vcat(sol.y, sol.y, sol.y, sol.x, sol.x, sol.x, p)
H_val = H(X, m)

out = similar(sol.y, m.n_x + m.n_y)
m.mod.m.H!(out, sol.y, sol.y, sol.y, sol.x, sol.x, sol.x, p)
H_val ≈ out

# calculate H_yp with forwarddiff
H_yp(X, m) = ForwardDiff.jacobian(X -> H(X, m), X)[1:(m.n_y + m.n_x), 1:(m.n_y)]
H_y(X, m) = ForwardDiff.jacobian(X -> H(X, m), X)[1:(m.n_y + m.n_x),
                                                  (m.n_y + 1):(2 * m.n_y)]
H_xp(X, m) = ForwardDiff.jacobian(X -> H(X, m), X)[1:(m.n_y + m.n_x),
                                                   (3 * m.n_y + 1):(3 * m.n_y + m.n_x)]
H_x(X, m) = ForwardDiff.jacobian(X -> H(X, m), X)[1:(m.n_y + m.n_x),
                                                  (3 * m.n_y + m.n_x + 1):(3 * m.n_y + 2 * m.n_x)]

H_yp_FD = H_yp(X, m)
out = zeros(m.n_y + m.n_x, m.n_y)
m.mod.m.H_yp!(out, sol.y, sol.x, p)
out ≈ H_yp_FD

## Finite Difference Utilies for testing
# using finite differences on the gradient calculations
using FiniteDiff
using FiniteDiff: finite_difference_derivative, finite_difference_gradient,
                  finite_difference_jacobian, finite_difference_hessian
modify_at_index(X, ind, val) = [X[1:(ind - 1)]; val; X[(ind + 1):end]]
# generates closure to make it univariate modifying only a single index
function make_univariate_at_index(f, X, ind, args...)
    return x -> f(modify_at_index(X, ind, x), args...)
end

# Derivatives of the H
p_deriv_index = (3 * m.n_y + 3 * m.n_x + 1) # first one is β, verfifying below
@test m_fvgq.mod.m.p_symbols[1] == :β

make_univariate_at_index(H_yp, X, p_deriv_index, m)(X[p_deriv_index])

@test c.H_yp_p[1] ≈
      finite_difference_derivative(make_univariate_at_index(H_yp, X,
                                                            p_deriv_index, m),
                                   X[p_deriv_index])
@test c.H_y_p[1] ≈
      finite_difference_derivative(make_univariate_at_index(H_y, X,
                                                            p_deriv_index, m),
                                   X[p_deriv_index])
@test c.H_xp_p[1] ≈
      finite_difference_derivative(make_univariate_at_index(H_xp, X,
                                                            p_deriv_index, m),
                                   X[p_deriv_index])
@test c.H_x_p[1] ≈
      finite_difference_derivative(make_univariate_at_index(H_x, X,
                                                            p_deriv_index, m),
                                   X[p_deriv_index])
# for the Ψ
function H_for_hessian(YX, ind, m, y_ss, x_ss, p)
    y_p = YX[1:(m.n_y)]
    y = YX[(m.n_y + 1):(2 * m.n_y)]
    x_p = YX[(2 * m.n_y + 1):(2 * m.n_y + m.n_x)]
    x = YX[(2 * m.n_y + m.n_x + 1):(2 * m.n_y + 2 * m.n_x)]

    out = similar(YX, m.n_x + m.n_y) # output 
    m.mod.m.H!(out, y_p, y, y_ss, x_p, x, x_ss, p)
    return out[ind]
end
YX = vcat(sol.y, sol.y, sol.x, sol.x)
H_for_hessian(YX, 2, m, sol.y, sol.x, p)
H_hessian(YX, ind, p) = ForwardDiff.hessian(YX -> H_for_hessian(YX, ind, m, sol.y, sol.x, p),
                                            YX)
Ψ_fd(p) = [H_hessian(YX, ind, p) for ind in 1:(m.n_x + m.n_y)]
@test all(c.Ψ .≈ Ψ_fd(p))

# For the Ψ_p
# Calculate it with finite differences over the forward-difference hessian.  Probably could nest AD with some work
Ψ_p_fd = finite_difference_derivative(make_univariate_at_index(Ψ_fd, p, 1), p[1])

# Needs to be zeroed out!
out = [zero(c.Ψ[1]) for _ in 1:(m.n_x + m.n_y)]  # only the single parameter is necessary.  Stacking hessians
m.mod.m.Ψ_p!(out, Val(:β), sol.y, sol.x, p)

# Compare to finite differences
@test all(out .≈ Ψ_p_fd)

# Steady state and gradients
@test c.H_p[1] ≈ finite_difference_derivative(make_univariate_at_index(H, X,
                                                                       p_deriv_index, m),
                                              X[p_deriv_index])

function ȳ(p)
    out = zeros(m.n_y)
    m.mod.m.ȳ!(out, p)
    return out
end
function x̄(p)
    out = zeros(m.n_x)
    m.mod.m.x̄!(out, p)
    return out
end

out = zero(c.y)
m.mod.m.ȳ_p!(out, Val(:β), p)
@test out ≈ finite_difference_derivative(make_univariate_at_index(ȳ, p, 1), p[1])

out = zero(c.x)
m.mod.m.x̄_p!(out, Val(:β), p)
@test out ≈ finite_difference_derivative(make_univariate_at_index(x̄, p, 1), p[1]) rtol = 1e-7

# Only test whether it can run
@test sol.retcode == :Success
# end

# @testset "FVGQ20 Kalman likelhood derivative in 1st order" begin

function joint_likelihood_1(p_d, p_f, m, observables, noise)
    sol = generate_perturbation(m, p_d, p_f)
    problem = LinearStateSpaceProblem(sol.A, sol.B, zeros(m.n_x), (0, size(observables, 2));
                                      sol.C,
                                      observables_noise = sol.D,
                                      noise, observables)
    return solve(problem, DirectIteration()).logpdf
end

# CRTU has problems with generating random MvNormal, so just testing diagonals
function kalman_likelihood(p_d, p_f, m, observables)
    sol = generate_perturbation(m, p_d, p_f)
    problem = LinearStateSpaceProblem(sol.A, sol.B, sol.x_ergodic,
                                      (0, size(observables, 2)); sol.C,
                                      observables_noise = sol.D,
                                      u0_prior = sol.x_ergodic,
                                      noise = nothing, observables)
    return solve(problem).logpdf
end
observables_fvgq = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
                                    "test/data/FVGQ20_observables.csv"), ',')' |> collect

noise_fvgq = readdlm(joinpath(pkgdir(DifferentiableStateSpaceModels),
                              "test/data/FVGQ20_noise.csv"),
                     ',')' |>
             collect

# m_fvgq already loaded
p_d = (β = 0.998, h = 0.97, ϑ = 1.17, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36),
       σ_μ = exp(-5.43), σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
p_f = (δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5), κ = 9.51)
# NOTE: Moved κ = 9.51 over to fixed.  Lower precision for likelihood.... 1e-4 etc.

kalman_likelihood(p_d, p_f, m_fvgq, observables_fvgq)
joint_likelihood_1(p_d, p_f, m_fvgq, observables_fvgq, noise_fvgq)
test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> kalman_likelihood(args..., p_f, m_fvgq, observables_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-6)
test_rrule(Zygote.ZygoteRuleConfig(),
           (args...) -> joint_likelihood_1(args..., p_f, m_fvgq, observables_fvgq,
                                           noise_fvgq), p_d;
           rrule_f = rrule_via_ad,
           check_inferred = false, rtol = 1e-6)

@test true

###############
# @testset "FVGQ20 Second Order" begin
const m_fvgq_2 = PerturbationModel(Main.FVGQ20)
p_d = (β = 0.998, h = 0.97, ϑ = 1.17, κ = 9.51, α = 0.21, θp = 0.82, χ = 0.63,
       γR = 0.77, γy = 0.19, γΠ = 1.29, Πbar = 1.01, ρd = 0.12, ρφ = 0.93, ρg = 0.95,
       g_bar = 0.3, σ_A = exp(-3.97), σ_d = exp(-1.51), σ_φ = exp(-2.36),
       σ_μ = exp(-5.43), σ_m = exp(-5.85), σ_g = exp(-3.0), Λμ = 3.4e-3, ΛA = 2.8e-3)
p_f = (δ = 0.025, ε = 10, ϕ = 0, γ2 = 0.001, Ω_ii = sqrt(1e-5))

c = SolverCache(m_fvgq_2, Val(2), p_d)
sol = generate_perturbation(m_fvgq_2, p_d, p_f, Val(2); cache = c)
generate_perturbation_derivatives!(m_fvgq_2, p_d, p_f, c)

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
