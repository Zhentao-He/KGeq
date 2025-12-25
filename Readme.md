$$
\begin{equation}
\Phi=0+\epsilon\Phi^{(1)}+\epsilon^2\Phi^{(2)}
+\mathcal{O}(\epsilon^3)
\end{equation}
$$

$$
\begin{equation}
\Box \Phi=\lambda^2\Phi
\end{equation}
$$

$$
\begin{equation}
\Box \Phi^{(1)}=0
\end{equation}
$$

In the hyperboloidal coordinates $(T,R,\vartheta,\phi)$[2010.00162]

$$
\begin{equation}
\Phi^{(1)}=\sum_m R\psi_m(T,R,\vartheta)
\exp(\mathrm{i}m\phi)
\end{equation}
$$

The master equation of $\psi_m(T,R,\vartheta)$ reads

$$
\begin{equation}
\left[ C_{TT}\partial _{T}^{2}+C_{TR}\partial _T\partial _R+C_{T\phi}\partial _T\partial _{\phi}+C_T\partial _T
+C_{RR}\partial _{R}^{2}+C_{R\phi}\partial _R\partial _{\phi}+C_R\partial _R+C_{\phi}\partial _{\phi}+C-_s\tilde{\Delta} \right] \psi_m=0
\end{equation}
$$

with $\partial _{\phi}\to\mathrm{i}m$.

Defining an auxiliary variable

$$
\begin{equation}
P=[C_{TT}\partial _{T}+C_{TR}\partial _R+C_{T\phi}\partial _{\phi}+C_T] \psi_m
\end{equation}
$$

to reduce the master equation as

$$
\begin{equation}
\begin{aligned}
    \partial_TP + \left[ 
C_{RR}\partial _{R}^{2}+C_{R\phi}\partial _R\partial _{\phi}+C_R\partial _R+C_{\phi}\partial _{\phi}+C-_s\tilde{\Delta} \right] \psi_m=0
\end{aligned}
\end{equation}
$$

# Pseudo-spectral method

[MatlCheb](https://github.com/Zhentao-He/MatlCheb) is required

求导矩阵
$R$ Chebyshev
$\vartheta$ spin-weighted harmonics

$$
\begin{equation}

\end{equation}
$$

exp filter

$$
\begin{equation}

\end{equation}
$$

# time-evolution

RK4
ode45

$$
\begin{equation}

\end{equation}
$$
