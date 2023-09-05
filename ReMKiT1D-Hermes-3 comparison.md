## Equation overview

### Electrons

$$\begin{align}
    &\frac{\partial n_e}{\partial t} + \frac{\partial \Gamma_e}{\partial x} = S,\\
    &m_e\frac{\partial \Gamma_e}{\partial t} + \frac{\partial}{\partial x} \left(n_ekT_e + m_e\Gamma_e u_e\right) + n_e e E = R_{ei} + S^{\Gamma}_e, \\
    &\frac{\partial W_e}{\partial t} + \frac{\partial}{\partial x} \left[\left(W_e + n_ekT_e\right)u_e + q_e\right] + \Gamma_e e E = Q_e, 
\end{align}$$

where $\Gamma_e=n_eu_e$ and $W_e = 3n_ekT_e/2 + n_e m_e u_e^2/2$.

$$q_e=-\kappa_e\nabla kT_e$$

$$R_{ei} = - 0.71n_e\nabla kT_e$$

$Q_e = Q_{h} + Q_{ei} + Q_{en}$, where $Q_{h}$ is the external heating.

For transport coefficients and sources see below.

### Ions 

$$\begin{align}
    &\frac{\partial n_i}{\partial t} + \frac{\partial \Gamma_i}{\partial x} = S,\\
    &m_i\frac{\partial \Gamma_i}{\partial t} + \frac{\partial}{\partial x} \left(n_ikT_e + m_e\Gamma_e u_e\right) - n_ie E = - R_{ei} + R_{CX} + S^{\Gamma}_i,  \\
    &\frac{\partial W_i}{\partial t} + \frac{\partial}{\partial x} \left[\left(W_i + n_ikT_i\right)u_i + q_i\right] + \Gamma_i e E = Q_i
\end{align}$$

where $R_{CX} = -\Gamma_in_nK_{CX}$ with $K_{CX}$ calculated using AMJUEL rate H.2-3.1.8, scaling the temperature dependence by $1/2$ to account for tracking deuterium instead of hydrogen.

$$q_i=-\kappa_i\nabla kT_i$$

$$Q_i = Q_{h} - Q_{ei} + Q_{in}$$

### Neutrals

$$\begin{align}
    &\frac{\partial n_n}{\partial t} + \frac{\partial \Gamma_n}{\partial x} = - S,\\
    &m_i\frac{\partial \Gamma_n}{\partial t} + \frac{\partial}{\partial x} \left(n_nkT_n + m_i\Gamma_n u_n\right) = - R_{CX} + S^{\Gamma}_n,  \\
    &\frac{\partial W_n}{\partial t} + \frac{\partial}{\partial x} \left[\left(W_n + n_nkT_n\right)u_n + q_n\right] = Q_n
\end{align}$$

where now $W_n=3nkT_n/2$ with $Q_n = -Q_{in}$. 

$$q_n = - \kappa_n \nabla kT_n$$

#### <span style="color:red">Hermes-3</span>
In 1D, the Hermes-3 neutral momentum equation is handled by this [general component](https://github.com/bendudson/hermes-3/blob/master/src/evolve_momentum.cxx) which also solves the ion momentum. There is an additional [component](https://github.com/bendudson/hermes-3/blob/master/src/neutral_mixed.cxx) which adds parallel diffusion and viscosity for 1D simulations. This one has a [documentation section](https://hermes3.readthedocs.io/en/latest/components.html#neutral-parallel-diffusion).

The following equation has been 

$$m \frac{\partial \Gamma_{n}} {\partial t} = -\frac{\partial}{\partial x} (u_{n} m \Gamma_{n}) 
 - \frac{\partial P_{n}}{\partial x}
 + \frac{\partial}{\partial x} (m \Gamma_n D_{n} \frac{1}{P_{n}} \frac{\partial P_{n}}{\partial x})
 + \frac{\partial}{\partial x} (\eta_{n} u_{n})
$$

Where the parallel projection of radial diffusivity is:
 $$D_{n} = dneut \frac{u_{th,n}^{2}}{\nu_{tot}} $$

And the viscosity is:
$$\eta_n = \frac{2}{5} \kappa_n = \frac{2}{5}n_n D_n$$

Currently the diffusion has no limiter in 1D.

### Electric field

Solve the same as in SOL-KiT, so acts mainly to balance electron pressure.


## Transport coefficients

The plasma transport coefficients are taken from [Makarov et al](https://pubs.aip.org/aip/pop/article/28/6/062308/973257/Equations-and-improved-coefficients-for-parallel). 

For $\kappa_{e,i}$ see the first term in equation (A27) in Makarov et al (see related equations A(15) and A(30)). To calculate (A30) for ions and electrons, the following Python code snippet was used 

```
ionZ = 1
sqrt2 = np.sqrt(2)

delta = 1 + (65*sqrt2/32 + 433*sqrt2/288 - 23*sqrt2/16)*ionZ + (5629/1152 - 529/128)*ionZ**2 #A30 in Makarov assuming single ion species and 0 mass ratio

elCondConst = 125*(1+433*sqrt2*ionZ/360)/(32*delta)
ionCondConst = 125/32
```

The Coulomb logarithms for e-e and i-i collisions were taken from the NRL Formulary (2013). 

For the e-i heat exchange term the standard form was used (see for example eq (19) in Makarov) with the NRL e-i Coulomb log.

For the neutrals the following heat conductivity is used 

$$\kappa_n = 2.4 n_n k T_n /(m_nn_i K_{CX})$$

## Sources and sinks

The particle source is given by 

$$S = n_e n_n K_{ion} - n_en_iK_{rec},$$

where the ionization and recombination rates are taken from AMJUEL (H.4.-2.1.5. and H.4.-2.1.8.). 

The momentum sources are given by the following:

$$ S^{\Gamma}_e = m_e(K_{ion} n_e \Gamma_n - K_{rec}n_e\Gamma_i),$$
where then $S_i^{\Gamma} = m_i S_e^{\Gamma}/m_e$ and $S_n^{\Gamma} = - S_i^{\Gamma}$.

For the electron energy source term due to e-n collisions 

$$Q_{en} = - n_e n_n K^E_{ion} - n_e n_i K^E_{rec} + n_en_i\epsilon_{ion}K_{rec}$$

where the first two terms are energy rates from AMJUEL (H.10.-2.1.5. and H.10.-2.1.8.) and the last term is the recombination energy source including the ionization potential $\epsilon_{ion}$.

For the ions, the energy source due to i-n interactions is

$$Q_{in} = (W_nn_i- W_in_n)K_{CX} + W_nn_eK_{ion} - W_in_e K_{rec}$$

## Boundary conditions

The upstream boundary condition is reflective for all species, with the heating term $Q_h$ specified over some length $L_h$ for the ions and electrons. 

At the sheath, the electron sheath heat transmission coefficients are 4.5, 2.5, and 0.25 for the electrons, ions, and neutrals respectively, with the neutral boundary condition being

$$ q_{n,sh} = \gamma_n n_n k T_n c_{s,n} = \gamma_n n_n k T_n \sqrt{kT_n/m_i} .$$

At the sheath recycling is set to 100%, and the ions and electron obey the Bohm condition. 