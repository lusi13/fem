# Solving the Heat Equation Using the Finite Element Method

To solve the Heat Equation with non-homogeneous, time-varying Dirichlet boundary conditions, we apply the **Finite Element Method**. 

First, we discretize space using a triangular mesh and Lagrange basis functions to approximate the horseshoe-shaped domain. Then, we solve the linear system and compute the time integral, discretizing time as well.

To achieve this, we store the coefficients derived from space discretization into the stiffness \( H \) and mass matrices \( M \), thereby obtaining a system of ODEs. Using the Crank-Nicholson scheme:

$$
\left(\frac{M}{\Delta t} + \theta H\right) u_{k+1} = \left(\frac{M}{\Delta t} - (1 - \theta)H\right) u_k,
$$

which is an unconditionally stable method for ODEs.

## Problem Definition

Given the following functions:

$$
\begin{aligned}
&f(x,t) : \Omega \times T \rightarrow \mathbb{R}, \\
&d_1(x,t) : \Gamma_{D,1} \times T \rightarrow \mathbb{R}, \\
&d_2(x,t) : \Gamma_{D,2} \times T \rightarrow \mathbb{R}, \\
&q(x,t) : \Gamma_N \times T \rightarrow \mathbb{R}, \\
&u_0(x) : \overline{\Omega} \rightarrow \mathbb{R}.
\end{aligned}
$$


we aim to find \( u(x,t) \in \Omega \times T \) such that:

$$
\begin{aligned}
    &\frac{\partial u(x, t)}{\partial t} - \frac{\partial^2 u(x, t)}{\partial x^2} - \frac{\partial^2 u(x, t)}{\partial y^2} = f(x, t) && \text{in } \Omega \times T, \\
    &u(x, t) = d_1(x, t) && \text{on } \Gamma_{D,1} \times T, \\
    &u(x, t) = d_2(x, t) && \text{on } \Gamma_{D,2} \times T, \\
    &\frac{\partial u(x, t)}{\partial n} = q(x, t) && \text{on } \Gamma_N \times T, \\
    &u(x, 0) = u_0(x) && \text{in } \Omega.
\end{aligned}
$$


*Figure: Geometry of the mesh. The Dirichlet time-varying boundary conditions affect only the edge "E1".*

---


