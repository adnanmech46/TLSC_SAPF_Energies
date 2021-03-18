# TLSC SPAF - Nonlinear Model - Closed-loop operation

This repository presents the implementation of the PI controllers to reproduce the results presented in the paper **Three-leg Split-Capacitor Shunt Active Power Filter presented by Urrea-Quintero et al. 2020 (Energies)**

---

<p align="center">
  <img align="middle" src="./Img/SAPFRobustControl.png" alt="Graphical abstract" height="450"/>
</p>

---

The TLSC SAPF model, after **dq0** transformation, is given by: 

<p align="center">
  <img align="middle" src="./Img/dq0SAPFModel.png" alt="SAPF model" height="150"/>
</p>

<!---
<img src="https://latex.codecogs.com/svg.latex?L\left[&space;\begin{bmatrix}&space;-\omega_{i_{S}}^{q}\\&space;\omega_{i_{S}}^{d}\\&space;0&space;\end{bmatrix}&space;&plus;&space;\dfrac{di_{S}^{dq0}}{dt}&space;\right]&space;=&space;-R_{L}i_{S}^{dq0}&space;-&space;\dfrac{v_{DC}}{2}u^{dq0}&space;-&space;\begin{bmatrix}&space;0&space;\\&space;0&space;\\&space;\sqrt{3}&space;\end{bmatrix}&space;\varepsilon_{v}&space;&plus;&space;v_{pcc}^{dq0}&space;\\&space;C\dfrac{dv_{DC}}{dt}&space;=&space;i_{S}^{dq0}u^{dq0^{T}}&space;-&space;\dfrac{v_{DC}}{R_{o}}\\&space;2C\dfrac{d\varepsilon_{v}}{dt}&space;=&space;\sqrt{3}i_{S}^0&space;-&space;\dfrac{2\varepsilon_{v}}{R_{o}}" title="L\left[ \begin{bmatrix} -\omega_{i_{S}}^{q}\\ \omega_{i_{S}}^{d}\\ 0 \end{bmatrix} + \dfrac{di_{S}^{dq0}}{dt} \right] = -R_{L}i_{S}^{dq0} - \dfrac{v_{DC}}{2}u^{dq0} - \begin{bmatrix} 0 \\ 0 \\ \sqrt{3} \end{bmatrix} \varepsilon_{v} + v_{pcc}^{dq0} \\ C\dfrac{dv_{DC}}{dt} = i_{S}^{dq0}u^{dq0^{T}} - \dfrac{v_{DC}}{R_{o}}\\ 2C\dfrac{d\varepsilon_{v}}{dt} = \sqrt{3}i_{S}^0 - \dfrac{2\varepsilon_{v}}{R_{o}}" />
--->

Three inner PI controllers are implemented to regulate SAPF active, in quadrature, and zero current components. An outer PI controller is implemented to indirectly regulate DC-link voltage through the direct current componen.  

---

<p align="center">
  <img align="middle" src="./Img/ControllersPerfComparison.png" alt="Controllers performance comparison" height="450"/>
</p>

---
