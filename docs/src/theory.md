## Frechet derivative of Matrix functions

Let $\mathcal{H}:=\mathbb{R}^{N\times N}_{\rm sym}$ be the vector space of $N\times N$ real symmetric matrices. For $H\in\mathcal{H}$ and $h_1,...h_n\in\mathcal{H}$, and $f$ is a $n$ times continuously differentiable function on a subset of $\mathbb{R}$ containing the spectrum of $H+t_1h_1+\cdots + t_nh_n$, the $n$ order Frechet derivative of $f(H)$ is
```math
\begin{align}
    \label{fd}
    \nonumber
    &{{\rm d}}^{n}f(H)h_1\cdots h_n =\\ &\frac{1}{2\pi i}\oint_\mathcal{C} f(z) \sum_{\mathcal{P}_n}(z-H)^{-1}h_{\mathcal{P}_n(1)}(z-H)^{-1}\cdots(z-H)^{-1}h_{\mathcal{P}_n(n)}(z-H)^{-1}\; dz,
\end{align}
```
where $\mathcal{C}$ is a contour in the complex plane enclosing all the eigenvalues of $H$, and $\mathcal{P}_n$ is an arbitrary permutation of $\{1,\cdots,n\}$. This can be proved by induction.

For $n = 1$, we have 
```math
\begin{align*}
    {\rm d} f(H)h&=\lim_{t\to 0} \frac{f(H+th)-f(H)}{t}\\&=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)\lim_{t\to 0}\frac{(z-H-th)^{-1}-(z-H)^{-1}}{t}\;dz.
\end{align*}
```
Note that
```math
\begin{align*}
    (z-H-th)^{-1}  &= (I-t(z-H)^{-1}h)^{-1}(z-H)^{-1} \\&= (I+t(z-H)^{-1}h+O(t^2))(z-H)^{-1},
\end{align*}
```
we have 
```math
\begin{align*}
    {\rm d} f(H)h=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)(z-H)^{-1}h(z-H)^{-1}\;dz,
\end{align*}
```
which satisfies the formula.

Assume the $n-1$ order derivative satisfies the formula. Then we have the $n$ order derivative
```math
\begin{align*}
    {{\rm d}}^{n}f(H)h_1\cdots h_n &=\lim_{t\to 0} \frac{{{\rm d}}^{n-1}f(H+th_n)h_1\cdots h_{n-1}-{{\rm d}}^{n-1}f(H)h_1\cdots h_{n-1}}{t}\\
    &=\frac{1}{2\pi i}\oint_\mathcal{C} f(z)\lim_{t\to 0}\sum_{\mathcal{P}_{n-1}}\frac{1}{t}\Big((z-H-th_n)^{-1}h_{\mathcal{P}_{n-1}(1)}\cdots h_{\mathcal{P}_{n-1}(n-1)}(z-H-th_n)^{-1}\\
    &\qquad -(z-H)^{-1}h_{\mathcal{P}_{n-1}(1)}\cdots h_{\mathcal{P}_{n-1}(n-1)}(z-H)^{-1}\Big)\;dz.
\end{align*}
```
Similarly, we have
```math 
\begin{align*}
   (z&-H-th_n)^{-1}h_{\mathcal{P}_{n-1}(1)}\cdots h_{\mathcal{P}(n-1)}(z-H-th_n)^{-1}\\
    &=(I+t(z-H)^{-1}h_n)(z-H)^{-1}h_{\mathcal{P}_{n-1}(1)}\cdots h_{\mathcal{P}_{n-1}(n-1)}(z-H)^{-1}(I+th_n(z-H)^{-1}) + O(t^2)\\
    &=(z-H)^{-1}h_{\mathcal{P}_{n-1}(1)}\cdots h_{\mathcal{P}_{n-1}(n-1)}(z-H)^{-1} \\
    &\quad+ t\Big((z-H)^{-1}h_n(z-H)^{-1}h_{\mathcal{P}_{n-1}(1)}\cdots h_{\mathcal{P}_{n-1}(n-1)}(z-H)^{-1}\\
    &\qquad\quad+(z-H)^{-1}h_{\mathcal{P}_{n-1}(1)}(z-H)^{-1}h_n\cdots h_{\mathcal{P}_{n-1}(n-1)}(z-H)^{-1}\\
    &\qquad\quad+(z-H)^{-1}h_{\mathcal{P}_{n-1}(1)}(z-H)^{-1}h_{\mathcal{P}_{n-1}(2)}\cdots h_n(z-H)^{-1}\Big) + O(t^2).
\end{align*}
```
Therefore, we can obtain the formula.

## Divided difference form
Let $(\lambda_i,\phi_i), i = 1,\dots,n$ be the eigenpairs of $H$, then we have
```math 
\begin{align*}
    (z-&H)^{-1}h_{\mathcal{P}_n(1)}(z-H)^{-1}\cdots(z-H)^{-1}h_{\mathcal{P}_n(n)}(z-H)^{-1}\\
    &=\sum_{i_1,\cdots,i_{n+1}=1}^n\phi_{i_1}(h_{\mathcal{P}_n(1)})_{i_1,i_2}\cdots (h_{\mathcal{P}_n(n)})_{i_n,i_{n+1}}\phi_{i_{n+1}}^*(z-\lambda_{i_1})^{-1}\cdots (z-\lambda_{i_{n+1}})^{-1},
\end{align*}
```
where $`(h_{\mathcal{P}_n(k)})_{i,j}=\phi_i^*h_{\mathcal{P}_n(k)}\phi_j`$. We let
```math 
    (z-\lambda_{i_1})^{-1}\cdots (z-\lambda_{i_{n+1}})^{-1} = \sum_{k=1}^{n+1}C_k (z-\lambda_{i_k})^{-1},
```
then
```math 
\sum_{k=1}^{n+1}C_k \prod_{\ell\neq k}(z-\lambda_{i_k})=1.
```
Let $z=\lambda_{i_k}$, we can obtain
```math 
C_k = \frac{1}{\prod_{\ell\neq k}(\lambda_{i_k}-\lambda_{i_\ell})}.
```
Therefore, we have
```math 
\begin{align*}
    \frac{1}{2\pi i}\oint_\mathcal{C} f(z)(z-\lambda_{i_1})^{-1}\cdots (z-\lambda_{i_{n+1}})^{-1}\;dz = \sum_{k=1}^{n+1}\frac{f(\lambda_{i_k})}{\prod_{\ell\neq k}(\lambda_{i_k}-\lambda_{i_\ell})}=f[\lambda_{i_1},\cdots,\lambda_{i_{n+1}}].
\end{align*}
```
Finally, we obtain
```math 
\begin{equation}
\label{dd1}
    {{\rm d}}^{n}f(H)h_1\cdots h_n =\sum_{i_1,\cdots,i_{n+1}=1}^n\phi_{i_1}\Bigg(\sum_{\mathcal{P}_n}(h_{\mathcal{P}_n(1)})_{i_1,i_2}\cdots (h_{\mathcal{P}_n(n)})_{i_n,i_{n+1}}\Bigg)f[\lambda_{i_1},\cdots,\lambda_{i_{n+1}}]\phi_{i_{n+1}}^*,
\end{equation}
```
which can be also written as 
```math 
\begin{equation}
\label{dd2}
(\Phi^*[{{\rm d}}^{n}f(H)h_1\cdots h_n]\Phi)_{k\ell}=\sum_{i_2,\cdots,i_n=1}^n\Bigg(\sum_{\mathcal{P}_n}(h_{\mathcal{P}_n(1)})_{k,i_2}\cdots (h_{\mathcal{P}_n(n)})_{i_n,\ell}\Bigg)f[k,\lambda_{i_2},\cdots,\lambda_{i_{n}},\lambda_\ell],
\end{equation}
```
where $\Phi = (\phi_1,\cdots,\phi_n)$.