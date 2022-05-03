module SW07
using LinearAlgebra, SymbolicUtils, LaTeXStrings
const max_order = 1
const n_y = 22
const n_x = 20
const n_p = 34
const n_ϵ = 2
const n_z = 42
const η = [0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0; 0.0 0.0]
const Q = LinearAlgebra.UniformScaling{Bool}(true)
const has_Ω = false
# Display definitions
const x_symbols = [:y_f_m, :y_m, :k_f, :k, :c_f_m, :c_m, :i_f_m, :i_m, :π_m, :w_m, :r_m, :ε_a, :b, :ε_g, :ε_i, :ε_r, :ε_pm, :η_p_aux, :ε_wm, :η_w_aux]
const y_symbols = [:k_s_f, :k_s, :r_f, :r, :r_k_f, :r_k, :z_f, :z, :w_f, :w, :l_f, :l, :q_f, :q, :y_f, :y, :i_f, :i, :c_f, :c, :π, :μ_pm]
const u_symbols = [:k_s_f, :k_s, :r_f, :r, :r_k_f, :r_k, :z_f, :z, :w_f, :w, :l_f, :l, :q_f, :q, :y_f, :y, :i_f, :i, :c_f, :c, :π, :μ_pm, :y_f_m, :y_m, :k_f, :k, :c_f_m, :c_m, :i_f_m, :i_m, :π_m, :w_m, :r_m, :ε_a, :b, :ε_g, :ε_i, :ε_r, :ε_pm, :η_p_aux, :ε_wm, :η_w_aux]
const p_symbols = [:ε_w, :ρ_ga, :ε_p, :l_bar, :Π_bar, :B, :μ_w, :μ_p, :α, :ψ, :φ, :δ, :σ_c, :λ, :ϕ_p, :ι_w, :ξ_w, :ι_p, :ξ_p, :σ_l, :ϕ_w, :r_π, :r_Δy, :r_y, :ρ, :ρ_a, :ρ_b, :ρ_g, :ρ_i, :ρ_r, :ρ_p, :ρ_w, :γ_bar, :gy_ss]
const H_latex = L"\begin{equation}
\left[
\begin{array}{c}
\alpha \mathrm{r\_k\_f}\left( t \right) - \mathrm{\varepsilon\_a}\left( t \right) + \left( 1 - \alpha \right) \mathrm{w\_f}\left( t \right) \\
\frac{\left( 1 - \psi \right) \mathrm{r\_k\_f}\left( t \right)}{\psi} - \mathrm{z\_f}\left( t \right) \\
 - \mathrm{k\_s\_f}\left( t \right) - \mathrm{r\_k\_f}\left( t \right) + \mathrm{l\_f}\left( t \right) + \mathrm{w\_f}\left( t \right) \\
 - \mathrm{k\_s\_f}\left( t \right) + \mathrm{k\_f}\left( t \right) + \mathrm{z\_f}\left( t \right) \\
\frac{\frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \mathrm{i\_f}\left( 1 + t \right)}{1 + \frac{1}{100} B} + \frac{\mathrm{q\_f}\left( t \right)}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{2} \varphi} + \mathrm{i\_f\_m}\left( t \right)}{1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B}} - \mathrm{i\_f}\left( t \right) + \mathrm{\varepsilon\_i}\left( t \right) \\
\frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}} \left( 1 - \delta \right) \mathrm{q\_f}\left( 1 + t \right)}{1 + \frac{1}{100} B} + \frac{\sigma_{c} \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right) b\left( t \right)}{1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}}} + \frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}} \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right) \mathrm{r\_k\_f}\left( 1 + t \right)}{1 + \frac{1}{100} B} - \mathrm{q\_f}\left( t \right) - \mathrm{r\_f}\left( t \right) \\
\frac{\lambda \mathrm{c\_f\_m}\left( t \right)}{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right)} + \frac{ - \left( 1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}} \right) \mathrm{r\_f}\left( t \right)}{\sigma_{c} \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right)} + \frac{\left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( -1 + \sigma_{c} \right) \left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right) \left(  - \mathrm{l\_f}\left( 1 + t \right) + \mathrm{l\_f}\left( t \right) \right)}{\alpha \sigma_{c} \phi_{w} \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right) \left( 1 - gy_{ss} - \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \right)} - \mathrm{c\_f}\left( t \right) + \frac{\mathrm{c\_f}\left( 1 + t \right)}{1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}}} + b\left( t \right) \\
\left( 1 - gy_{ss} - \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \right) \mathrm{c\_f}\left( t \right) - \mathrm{y\_f}\left( t \right) + \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right) \mathrm{z\_f}\left( t \right) + \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \mathrm{i\_f}\left( t \right) + \mathrm{\varepsilon\_g}\left( t \right) \\
 - \mathrm{y\_f}\left( t \right) + \phi_{p} \left( \alpha \mathrm{k\_s\_f}\left( t \right) + \left( 1 - \alpha \right) \mathrm{l\_f}\left( t \right) + \mathrm{\varepsilon\_a}\left( t \right) \right) \\
\sigma_{l} \mathrm{l\_f}\left( t \right) + \frac{ - \lambda \mathrm{c\_f\_m}\left( t \right)}{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}} \right)} + \frac{\mathrm{c\_f}\left( t \right)}{1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}}} - \mathrm{w\_f}\left( t \right) \\
\left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \mathrm{i\_f}\left( t \right) + \frac{\left( 1 - \delta \right) \mathrm{k\_f}\left( t \right)}{1 + \frac{1}{100} \gamma_{bar}} - \mathrm{k\_f}\left( 1 + t \right) + \left( 1 + \frac{1}{100} \gamma_{bar} \right)^{2} \varphi \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \mathrm{\varepsilon\_i}\left( t \right) \\
 - \mathrm{\varepsilon\_a}\left( t \right) + \alpha \mathrm{r\_k}\left( t \right) + \left( 1 - \alpha \right) w\left( t \right) - \mathrm{\mu\_pm}\left( t \right) \\
\frac{\left( 1 - \psi \right) \mathrm{r\_k}\left( t \right)}{\psi} - z\left( t \right) \\
 - \mathrm{k\_s}\left( t \right) - \mathrm{r\_k}\left( t \right) + l\left( t \right) + w\left( t \right) \\
 - \mathrm{k\_s}\left( t \right) + k\left( t \right) + z\left( t \right) \\
\frac{\frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right) i\left( 1 + t \right)}{1 + \frac{1}{100} B} + \frac{q\left( t \right)}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{2} \varphi} + \mathrm{i\_m}\left( t \right)}{1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B}} - i\left( t \right) + \mathrm{\varepsilon\_i}\left( t \right) \\
\frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}} \left( 1 - \delta \right) q\left( 1 + t \right)}{1 + \frac{1}{100} B} + \frac{\sigma_{c} \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right) b\left( t \right)}{1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}}} + \frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}} \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right) \mathrm{r\_k}\left( 1 + t \right)}{1 + \frac{1}{100} B} - q\left( t \right) - r\left( t \right) + \pi\left( 1 + t \right) \\
\frac{\lambda \mathrm{c\_m}\left( t \right)}{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right)} + \frac{ - \left( 1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}} \right) \left(  - \pi\left( 1 + t \right) + r\left( t \right) \right)}{\sigma_{c} \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right)} + \frac{\left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( -1 + \sigma_{c} \right) \left( 1 - \alpha \right) \left(  - l\left( 1 + t \right) + l\left( t \right) \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\alpha \sigma_{c} \phi_{w} \left( 1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}} \right) \left( 1 - gy_{ss} - \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \right)} - c\left( t \right) + \frac{c\left( 1 + t \right)}{1 + \frac{\lambda}{1 + \frac{1}{100} \gamma_{bar}}} + b\left( t \right) \\
 - y\left( t \right) + \left( 1 - gy_{ss} - \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \right) c\left( t \right) + \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right) z\left( t \right) + \left( \frac{\left( 1 - \alpha \right) \left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)}{\left( \frac{\left( 1 - \alpha \right)^{1 - \alpha} \alpha^{\alpha}}{\left( -1 + \delta + \frac{1 + \frac{1}{100} B}{\left( 1 + \frac{1}{100} \gamma_{bar} \right)^{ - \sigma_{c}}} \right)^{\alpha} \phi_{p}} \right)^{\frac{1}{1 - \alpha}} \alpha} \right)^{-1 + \alpha} \phi_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) i\left( t \right) + \mathrm{\varepsilon\_g}\left( t \right) \\
 - y\left( t \right) + \phi_{p} \left( \alpha \mathrm{k\_s}\left( t \right) + \left( 1 - \alpha \right) l\left( t \right) + \mathrm{\varepsilon\_a}\left( t \right) \right) \\
\frac{\iota_{p} \mathrm{\pi\_m}\left( t \right) + \frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \pi\left( 1 + t \right)}{1 + \frac{1}{100} B} + \frac{\left( 1 - \xi_{p} \right) \left( 1 + \frac{ - \iota_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right)}{1 + \frac{1}{100} B} \right) \mathrm{\mu\_pm}\left( t \right)}{\xi_{p} \left( 1 + \varepsilon_{p} \left( -1 + \phi_{p} \right) \right)}}{1 + \frac{\iota_{p} \left( 1 + \frac{1}{100} \gamma_{bar} \right)}{1 + \frac{1}{100} B}} - \pi\left( t \right) + \mathrm{\varepsilon\_pm}\left( t \right) \\
\frac{\iota_{w} \mathrm{\pi\_m}\left( t \right)}{1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B}} + \frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right) w\left( 1 + t \right)}{\left( 1 + \frac{1}{100} B \right) \left( 1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B} \right)} + \frac{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \pi\left( 1 + t \right)}{\left( 1 + \frac{1}{100} B \right) \left( 1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B} \right)} + \frac{ - \left( 1 + \frac{\iota_{w} \left( 1 + \frac{1}{100} \gamma_{bar} \right)}{1 + \frac{1}{100} B} \right) \pi\left( t \right)}{1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B}} + \frac{\left( 1 - \xi_{w} \right) \left( 1 + \frac{ - \xi_{w} \left( 1 + \frac{1}{100} \gamma_{bar} \right)}{1 + \frac{1}{100} B} \right) \left( \sigma_{l} l\left( t \right) + \frac{ - \lambda \mathrm{c\_m}\left( t \right)}{\left( 1 + \frac{1}{100} \gamma_{bar} \right) \left( 1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}} \right)} + \frac{c\left( t \right)}{1 + \frac{ - \lambda}{1 + \frac{1}{100} \gamma_{bar}}} - w\left( t \right) \right)}{\xi_{w} \left( 1 + \varepsilon_{w} \left( -1 + \phi_{w} \right) \right) \left( 1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B} \right)} - w\left( t \right) + \frac{\mathrm{w\_m}\left( t \right)}{1 + \frac{1 + \frac{1}{100} \gamma_{bar}}{1 + \frac{1}{100} B}} + \mathrm{\varepsilon\_wm}\left( t \right) \\
r_{{\Delta}y} \left(  - \mathrm{y\_f}\left( t \right) - \mathrm{y\_m}\left( t \right) + y\left( t \right) + \mathrm{y\_f\_m}\left( t \right) \right) + \rho \mathrm{r\_m}\left( t \right) + r_{y} \left( 1 - \rho \right) \left(  - \mathrm{y\_f}\left( t \right) + y\left( t \right) \right) + r_{\pi} \left( 1 - \rho \right) \pi\left( t \right) + \mathrm{\varepsilon\_r}\left( t \right) \\
\left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) i\left( t \right) + \frac{\left( 1 - \delta \right) k\left( t \right)}{1 + \frac{1}{100} \gamma_{bar}} - k\left( 1 + t \right) + \left( 1 + \frac{1}{100} \gamma_{bar} \right)^{2} \varphi \left( 1 + \frac{-1 + \delta}{1 + \frac{1}{100} \gamma_{bar}} \right) \mathrm{\varepsilon\_i}\left( t \right) \\
 - \mathrm{\varepsilon\_a}\left( 1 + t \right) + \rho_{a} \mathrm{\varepsilon\_a}\left( t \right) \\
 - b\left( 1 + t \right) + \rho_{b} b\left( t \right) \\
 - \mathrm{\varepsilon\_g}\left( 1 + t \right) + \rho_{g} \mathrm{\varepsilon\_g}\left( t \right) \\
 - \mathrm{\varepsilon\_i}\left( 1 + t \right) + \rho_{i} \mathrm{\varepsilon\_i}\left( t \right) \\
\rho_{r} \mathrm{\varepsilon\_r}\left( t \right) - \mathrm{\varepsilon\_r}\left( 1 + t \right) \\
\rho_{p} \mathrm{\varepsilon\_pm}\left( t \right) - \mathrm{\varepsilon\_pm}\left( 1 + t \right) - \mu_{p} \mathrm{\eta\_p\_aux}\left( t \right) + \mathrm{\eta\_p\_aux}\left( 1 + t \right) \\
 - \mathrm{\eta\_p\_aux}\left( 1 + t \right) \\
\rho_{w} \mathrm{\varepsilon\_wm}\left( t \right) - \mathrm{\varepsilon\_wm}\left( 1 + t \right) - \mu_{w} \mathrm{\eta\_w\_aux}\left( t \right) + \mathrm{\eta\_w\_aux}\left( 1 + t \right) \\
 - \mathrm{\eta\_w\_aux}\left( 1 + t \right) \\
 - \mathrm{y\_f\_m}\left( 1 + t \right) + \mathrm{y\_f}\left( t \right) \\
 - \mathrm{y\_m}\left( 1 + t \right) + y\left( t \right) \\
 - \mathrm{c\_f\_m}\left( 1 + t \right) + \mathrm{c\_f}\left( t \right) \\
 - \mathrm{c\_m}\left( 1 + t \right) + c\left( t \right) \\
 - \mathrm{i\_f\_m}\left( 1 + t \right) + \mathrm{i\_f}\left( t \right) \\
 - \mathrm{i\_m}\left( 1 + t \right) + i\left( t \right) \\
 - \mathrm{\pi\_m}\left( 1 + t \right) + \pi\left( t \right) \\
 - \mathrm{w\_m}\left( 1 + t \right) + w\left( t \right) \\
 - \mathrm{r\_m}\left( 1 + t \right) + r\left( t \right) \\
\end{array}
\right]
\end{equation}
"
const steady_states_latex = L"\begin{align}
\mathrm{y\_f\_m}\left( \infty \right) =& 0 \\
\mathrm{y\_m}\left( \infty \right) =& 0 \\
\mathrm{k\_f}\left( \infty \right) =& 0 \\
k\left( \infty \right) =& 0 \\
\mathrm{c\_f\_m}\left( \infty \right) =& 0 \\
\mathrm{c\_m}\left( \infty \right) =& 0 \\
\mathrm{i\_f\_m}\left( \infty \right) =& 0 \\
\mathrm{i\_m}\left( \infty \right) =& 0 \\
\mathrm{\pi\_m}\left( \infty \right) =& 0 \\
\mathrm{w\_m}\left( \infty \right) =& 0 \\
\mathrm{r\_m}\left( \infty \right) =& 0 \\
\mathrm{\varepsilon\_a}\left( \infty \right) =& 0 \\
b\left( \infty \right) =& 0 \\
\mathrm{\varepsilon\_g}\left( \infty \right) =& 0 \\
\mathrm{\varepsilon\_i}\left( \infty \right) =& 0 \\
\mathrm{\varepsilon\_r}\left( \infty \right) =& 0 \\
\mathrm{\varepsilon\_pm}\left( \infty \right) =& 0 \\
\mathrm{\eta\_p\_aux}\left( \infty \right) =& 0 \\
\mathrm{\varepsilon\_wm}\left( \infty \right) =& 0 \\
\mathrm{\eta\_w\_aux}\left( \infty \right) =& 0 \\
\mathrm{k\_s\_f}\left( \infty \right) =& 0 \\
\mathrm{k\_s}\left( \infty \right) =& 0 \\
\mathrm{r\_f}\left( \infty \right) =& 0 \\
r\left( \infty \right) =& 0 \\
\mathrm{r\_k\_f}\left( \infty \right) =& 0 \\
\mathrm{r\_k}\left( \infty \right) =& 0 \\
\mathrm{z\_f}\left( \infty \right) =& 0 \\
z\left( \infty \right) =& 0 \\
\mathrm{w\_f}\left( \infty \right) =& 0 \\
w\left( \infty \right) =& 0 \\
\mathrm{l\_f}\left( \infty \right) =& 0 \\
l\left( \infty \right) =& 0 \\
\mathrm{q\_f}\left( \infty \right) =& 0 \\
q\left( \infty \right) =& 0 \\
\mathrm{y\_f}\left( \infty \right) =& 0 \\
y\left( \infty \right) =& 0 \\
\mathrm{i\_f}\left( \infty \right) =& 0 \\
i\left( \infty \right) =& 0 \\
\mathrm{c\_f}\left( \infty \right) =& 0 \\
c\left( \infty \right) =& 0 \\
\pi\left( \infty \right) =& 0 \\
\mathrm{\mu\_pm}\left( \infty \right) =& 0
\end{align}
"
const steady_states_iv_latex = L"$$"
# Function definitions
include("SW07/zero_order_ip.jl")
include("SW07/first_order_ip.jl")
end
