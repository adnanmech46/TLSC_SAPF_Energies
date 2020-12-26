# TLSC SPAF - Nonlinear Model - Closed-loop operation

This notebook presents the implementation of the PI controllers based current control to the nonlinear model of the **Three-leg Split-Capacitor Shunt Active Power Filter presented by Urrea-Quintero et al. 2020 (Energies)**

The TLSC SAPF model, after $dq0$ transformation, is given by: 

$$L\left[ \begin{bmatrix} -\omega_{i_{S}}^{q}\\
\omega_{i_{S}}^{d}\\
0
\end{bmatrix} + \dfrac{di_{S}^{dq0}}{dt} \right]  =  -R_{L}i_{S}^{dq0} - \dfrac{v_{DC}}{2}u^{dq0} - \begin{bmatrix}
0 \\
0 \\
\sqrt{3}
\end{bmatrix} \varepsilon_{v} + v_{pcc}^{dq0}$$
$$C\dfrac{dv_{DC}}{dt} = i_{S}^{dq0}u^{dq0^{T}} - \dfrac{v_{DC}}{R_{o}}$$
$$2C\dfrac{d\varepsilon_{v}}{dt} = \sqrt{3}i_{S}^0 - \dfrac{2\varepsilon_{v}}{R_{o}}$$
