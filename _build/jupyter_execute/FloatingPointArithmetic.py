#!/usr/bin/env python
# coding: utf-8

# (CH-Fl-Po-Ar)=
# # Floating Point Arithmetic
# 
# 
# In this chapter, we will introduce some basics on the real number system for modern computers and discuss the arithmetic operations of the number system.
# 
# (SEC-Re-Re-Nu)=
# ## Representation of Real Numbers
# 
# Any nonzero real number $x\in \mathbb{R}$ can be accurately represented with an infinite sequence of digits. This can be understood as the consequence that rational numbers are dense on any interval. Therefore, with the binary representation, we can write 
# 
# $$
# x = \pm (0.d_1 d_2 d_3\dots d_{t-1} d_t d_{t+1}\dots) \times 2^e,
# $$
# 
# where $e$ is an integer exponent and $d_1=1$, the other binary digits $d_i\in \{0, 1\}$. The mantissa part
# 
# $$
# 0.d_1 d_2 d_3\dots =  \frac{d_1}{2} + \frac{d_2}{2^2} + \frac{d_3}{2^3} + \cdots.
# $$

# ```{note}
# In order to guarantee the uniqueness of the above representation, we need further assumption that there exists an infinite subset $S\subset \mathbb{N}$ that $d_j\neq 1$ for $j\in S$. For example, under binary representation
# $0.111\dots = (0.1)\times 2^1,$
# then we will take the latter representation.
# ```
# 

# (SEC-Fl-Po-Nu)=
# ## Floating Point Numbers
# The floating point numbers generally refer to a set of real numbers with **finite** mantissa length. More precisely, we consider the set of real numbers $\mathbb{F} = \mathbb{F}(t, e_{\min}, e_{\max})\subset \mathbb{R}$ that
# 
# $$\mathbb{F} := \{x\in \mathbb{R} \mid x = \pm (0.d_1 d_2 d_3\dots d_{t-1} d_t) \times 2^e, d_1 =1,  e_{\min}\le e\le e_{\max}\}\cup \{0\}.$$
# 
# It can be seen that there are only finite numbers in $\mathbb{F}$ with the smallest positive element $x_{\min} = 2^{e_{\min}-1}$ and the largest element $x_{\max} = ( 1- 2^{-t} )\times 2^{e_{\max} }$. Therefore
# 
# $$\mathbb{F}\subset \overline{\mathbb{F}}:= \{ x\in \mathbb{R}\mid x_{\min} \le |x| \le x_{\max}\}\cup \{ 0\}.$$

# ```{note}
# The elements in $\mathbb{F}$ are called normalized. If we  allow $d_1 = 0$ in the definition of $\mathbb{F}$, then the numbers in the set are called denormalized.
# ```

# ````{prf:theorem} Distribution of Floating Numbers
# :label: THM-Di-Fl-Nu
# 
# For any $e_{\min} \le e\le e_{\max}$, the distribution of the floating point number system $\mathbb{F}$ on interval $[2^{e-1}, 2^e]$ is equidistant with distances of length $2^{e-t}$.
# 
# ````

# 
# ````{prf:proof}
# For any  $x \in \mathbb{F}\cap [2^{e-1}, 2^e]$, it can be represented by 
# 
# $$x = (0.d_1d_2\dots d_t)\times 2^{e}$$
# 
# where $d_1 = 1$. The mantissa is equidistantly distributed with distance $2^{-t}$, therefore the floating point numbers are equidistantly distributed with distances of length $2^{e-t}$. 
# ````

# To understand the approximation to real numbers by the floating point number system $\mathbb{F}$, it is important to consider the maximal relative distance between the numbers in $\overline{\mathbb{F}}$ and their respective closest element in $\mathbb{F}$, which is the following quantity:
# 
# $$\max_{x\in \overline{\mathbb{F}}, x\neq 0}\min_{z\in\mathbb{F}} \frac{|z - x|}{|x|}.$$
# 
# The following holds:
# 
# ````{prf:theorem} Machine Precision
# :label: THM-Ma-Pr
# 
# $$\max_{x\in \overline{\mathbb{F}}, x\neq 0}\min_{z\in\mathbb{F}} \frac{|z - x|}{|x|}\le  2^{-t}.$$
# 
# The number $\mathrm{u} := 2^{-t}$ is also called rounding unit or machine precision.
# ````
# 
# ````{prf:proof}
# Without loss of generality, we only need to consider the positive numbers in $\overline{\mathbb{F}}$, then one can represent any nonzero $x\in [x_{\min}, x_{\max}]$ by
# 
# $$x = (0.d_1 d_2 \dots d_t\dots)\times 2^e\in [2^{e-1}, 2^e].$$
# 
# Since the floating point numbers are equidistantly distributed on $[2^{e-1}, 2^e]$ from {prf:ref}``THM-Di-Fl-Nu``, one can find $z^{\ast}\in\mathbb{F}$ such that 
# 
# $$|z^{\ast} - x| \le \frac{1}{2} 2^{e-t},$$
# 
# therefore 
# 
# $$\frac{|z^{\ast}- x|}{|x|}\le \frac{1}{2} 2^{e-t} \frac{1}{2^{e-1}} = 2^{-t}.$$
# ````
# 
# ```{note}
# On modern computers, the following two floating point number systems 
# 
# $$\mathbb{F}_{32} := \mathbb{F}(24, -125, 128),\quad\mathbb{F}_{64} := \mathbb{F}(53, -1021, 1024)$$ 
# 
# are supported, they are often called single precision and double precision, respectively.  
# ```

# (SEC:ROUNDING)=
# ## Rounding
# The rounding operation $\textrm{fl}$ is to map any real numbers of $\overline{\mathbb{F}}$ into the floating point number system $\mathbb{F}$ with smallest error. Such rounding operation can be written out explicitly, let $x = \pm (0.d_1 d_2\dots d_t d_{t+1}\dots )\times 2^e$, then
# 
# $$
#     \textrm{fl}(x) =\begin{cases}
#         \pm (0.d_1 d_2\dots d_t) \times 2^e & \text{if } d_{t+1} = 0,\\
#         \pm (0.d_1 d_2\dots d_t  + 2^{-t}) \times 2^e & \text{if } d_{t+1} = 1.
#     \end{cases}
# $$
# 
# It is clear that rounding $\textrm{fl}$ is monotone and idempotent, which means
# 
# - $x\le y \Rightarrow \textrm{fl}(x) \le \textrm{fl}(y)$. 
# - $\textrm{fl}(z) = z$ if $z\in \mathbb{F}$.
# 
# ````{prf:theorem}
# :label: THM-REL-ERR
# For any $x\in \overline{\mathbb{F}}$, $|\textrm{fl}(x) - x| = \min_{z\in \mathbb{F}} |z - x|$. If $x\neq 0$, then
# 
# $$\frac{|\textrm{fl}(x) - x|}{|x|}\le \mathrm{u} = 2^{-t}.$$
# 
# ````
# 
# ````{prf:proof}
# The special case that $x = 0$ is trivial, we only consider $x\in [x_{\min}, x_{\max}]$, it can be seen that 
# 
# $$
#     |\textrm{fl}(x) - x| = | (0.d_1d_2\dots \tilde{d}_t) - (0.d_1d_2\dots d_t d_{t+1} \dots )|\times 2^e \le 2^{-(t+1)} \times 2^e,
# $$
# 
# where $\tilde{d}_t$ is the rounding bit, therefore 
# 
# $$
#     \frac{|\textrm{fl}(x) - x|}{|x|} \le \frac{2^{e-(t+1)}}{2^{e-1}} = 2^{-t}.
# $$
# ````
# 
# ```{prf:corollary}
# :label: COR-REL-ERR
# For any $x\in\overline{\mathbb{F}}$, $\textrm{fl}(x) = x(1+\delta)$ with $|\delta|\le \mathrm{u}$.
# ```

# (SEC-AR-OP)=
# ## Arithmetic Operations
# Let $\mathbb{F} = \mathbb{F}(t, e_{\min}, e_{\max})$ be a given floating point number system and we consider the basic binary operation $\circ\in \{ +, -, *, /\}$ on $\mathbb{F}$, in order to represent the outcome in $\mathbb{F}$, a straightforward realization is to define the binary operation $\boxed{\circ}$ as the following (for the case of division, we assume $y\neq 0$):
# 
# $$x \,\boxed{\circ}\, y := \textrm{fl}(x\circ y),$$
# 
# then for any $x, y\in\mathbb{F}$, if $x\circ y\in \overline{\mathbb{F}}$, then $x \,\boxed{\circ}\, y = (x\circ y)(1 + \delta)$ with $|\delta| \le \mathrm{u}$ from the {prf:ref}``COR-REL-ERR``.
# 
# ```{prf:remark} Cancellation Error
# If $x, y\in\mathbb{R}$, then the relative error from the following binary operation $\textrm{fl}(x)\,\boxed{+}\,\textrm{fl}(y)$ can be estimated by
# 
# $$
# \begin{aligned}
# \frac{|\textrm{fl}(x)\,\boxed{+}\,\textrm{fl}(y)- (x +  y) |}{|x + y|} &\le \frac{|\textrm{fl}(x)\,\boxed{+}\,\textrm{fl}(y) - (\textrm{fl}(x) +  \textrm{fl}(y)) |}{|x + y|} + \frac{| (\textrm{fl}(x) +  \textrm{fl}(y)) - (x+y)   |}{|x+y|}\\  
# &\le\mathrm{u} + (\mathrm{u} +\mathrm{u}^2) \frac{|x|+|y|}{|x+y|}.
# \end{aligned}
# $$
# 
# When $x$ and $y$ are close in magnitude but with opposite signs, the cancellation error will be significant.
# ```
# 
# (SSEC-ER-AC-MU)=
# ### Error Accumulation: Multiplication
# For complicated computations on modern computers, the errors from arithmetic operations will accumulate towards the final result (we do not consider techniques such as fused multiply-add (FMA) here). To quantify the accumulation effect, we will need the following lemma.
# 
# 
# ```{prf:lemma} 
# :label: LEM-ACC
# For real numbers $a_1, a_2, \dots, a_n$ with $|a_k|\le \delta$ for $k=1,\dots, n$, then for $n\delta < 1$, the following holds
# 
# $$\prod_{k=1}^n (1 + a_k) = 1 + b_n,$$
# 
# where $|b_n| \le \frac{n\delta}{1 - n\delta}$.
# ```
# 
# ```{prf:proof}
# The proof is quite easy with induction. When $n=1$, $|b_1| = |a_1|\le \delta \le \frac{\delta}{1 - \delta}$. Suppose the claim holds for $n = m$, then for $n = m +1$, we could see that
# 
# $$
#     \prod_{k=1}^{m+1} (1 + a_k) = (1 + b_m) (1 + a_{m+1}) = 1 + b_{m+1}, 
# $$
# 
# which implies that $b_{m+1} = b_m + a_{m+1} + a_{m+1}b_m$, with the given bounds on $a_{m+1}$ and $b_m$, we can estimate
# 
# $$
#     |b_{m+1}| = | b_m + a_{m+1} + a_{m+1}b_m | \le \frac{(m+1)\delta}{1 - m\delta} \le \frac{(m+1)\delta}{1 - (m+1)\delta}.
# $$
# ```
# 
# 
# In the following, we consider the naive floating point product $P_n$ of $n$ real numbers $\{x_j\}_{j=1}^n \subset \mathbb{R}$ with assumption that $(2n-1)\textrm{u} < 1$ by the following iteration
# 
# $$\newcommand{\fl}[1]{\textrm{fl}(#1)}
#     \begin{cases}
#         P_k =\fl{x_1} & k=1,\\
#         P_k = P_{k-1}~\boxed{\ast}~ \fl{x_k} & k\ge 2.
#     \end{cases}
# $$
# 
# Let $\fl{x_k} = x_k(1 + \tau_k)$, then $|\tau_k|\le \textrm{u}$. From the $n$-th iteration step 
# 
# $$
# P_n = P_{n-1}~\boxed{\ast}~ \fl{x_n} = \fl{\fl{P_{n-1}}\fl{x_n}} = \fl{P_{n-1}}\fl{x_n}(1 + \delta_n)
# $$ 
# 
# such that $|\delta_n|\le \textrm{u}$. Since $P_{n-1}\in \mathbb{F}$, $\fl{P_{n-1}} = P_{n-1}$, then 
# 
# $$
# \begin{aligned}
#     P_n &= P_{n-1} \fl{x_n}(1 + \delta_n) = P_{n-2} \fl{x_{n-1}}(1 + \delta_{n-1}) (1 + \delta_n) = \cdots \\
#     &=  \fl{x_1} \fl{x_2}\cdots \fl{x_n} \prod_{j=2}^n (1 + \delta_j) \\
#     &= \prod_{j=1}^n x_j (1 + \tau_j)  \prod_{j=2}^n (1 + \delta_j) \\
#     &\le \left(\prod_{j=1}^n x_j \right) (1 + \eta_n),
# \end{aligned}
# $$
# where $|\eta_n|\le \frac{(2n-1)\textrm{u}}{1 - (2n-1)\textrm{u}}$. 
# 
# (SSEC-ER-AC-AD)=
# ### Error Accumulation: Addition
# For the naive floating point summation $S_n$ of $n$ real numbers $\{x_j \}_{j=1}^n$ by the iteration 
# 
# $$
# \begin{cases}
#     S_k = 0 & k=0,\\
#     S_k = S_{k-1}~\boxed{+}~ \fl{x_k} & k\ge 1,
# \end{cases}
# $$
# 
# we can carry out a similar analysis. Let $S^{\ast}_j = \sum_{k=1}^j x_k$ and  $\fl{x_k} = x_k ( 1 + \tau_k)$ for $|\tau_k|\le \textrm{u}$, denote $\Delta S_j =     S_j^{\ast} - S_j$, then
# 
# $$
# \begin{aligned}
#     \Delta S_j &= S_j^{\ast} - S_j \\&= S_{j}^{\ast} - ( S_{j-1}~\boxed{+}~ \fl{x_j}) \\
#     &=  S_{j}^{\ast}  - ( S_{j}^{\ast} - \Delta S_{j-1} + x_j\tau_j)(1 + \delta_j) \\
#     &= \Delta S_{j-1}(1 + \delta_j) - \delta_j S_{j}^{\ast} - x_j(\tau_j)(1+\delta_j),
# \end{aligned}
# $$
# 
# where $|\delta_j|\le \textrm{u}$. Therefore
# 
# $$
# \begin{aligned}
#     |\Delta S_j|&\le |\Delta S_{j-1}|(1 + \textrm{u}) + (\sum_{k=1}^j |x_k|)\textrm{u} + |x_j| \textrm{u}(1 + \textrm{u}) \\
#     &\le |\Delta S_{j-2}|(1 + \textrm{u})^2 + (\sum_{k=1}^j |x_k|)\textrm{u}(1+\textrm{u}) + |x_j| \textrm{u}(1 + \textrm{u})^2 \\
#     &\quad\quad +  (\sum_{k=1}^{j-1} |x_k|)\textrm{u} + |x_{j-1}| \textrm{u}(1 + \textrm{u})\\
#     &\le |\Delta S_{j-3}|(1 + \textrm{u})^3 + (\sum_{k=1}^j |x_k|)\textrm{u}(1+\textrm{u})^2 + |x_j| \textrm{u}(1 + \textrm{u})^3 \\
#     &\quad\quad +  (\sum_{k=1}^{j-1} |x_k|)\textrm{u}(1 + \textrm{u}) + |x_{j-1}| \textrm{u}(1 + \textrm{u})^2 \\
#     &\quad\quad +  (\sum_{k=1}^{j-2} |x_k|)\textrm{u}+ |x_{j-2}| \textrm{u}(1 + \textrm{u})\\
#     &\le  \cdots \\
#     &\le  \sum_{l=1}^{j} \left( (\sum_{k=1}^l |x_k|) \textrm{u}(1 + \textrm{u})^{l-1} + |x_{l}|\text{u}(1 + \textrm{u})^{l-1} \right) \\
#     & \le \sum_{l=1}^{j} \left( (\sum_{k=1}^j |x_k|) \textrm{u}(1 + \textrm{u})^{l-1} \right) + \left(\sum_{l=1}^j |x_{l}|\right)\text{u}(1 + \textrm{u})^{j-1}\\
#     &= \left(\sum_{l=1}^j |x_{l}|\right) ((1 + \textrm{u})^j - 1)
# \end{aligned}
# $$
# 
# Using the previous lemma will provide an estimate of $ ((1 + \textrm{u})^j - 1)$ as long as $j\textrm{u} < 1$.
# 
# (SEC-CH1-EXER)=
# ## Exercises
# 
# Assume $n\in \mathbb{N}$ and $n\textrm{u} < 1$, let $\{x_j\}_{j=1}^n$ be a sequence of real numbers, $\textrm{fl}:\overline{\mathbb{F}}\mapsto \mathbb{F}$ is the rounding operation, $\textrm{u}$ is the machine epsilon. In the following, we briefly discuss the rounding error for floating point summation (sum reduction). 
# 
# ````{prf:definition}
# A reduction $\Pi$ of the floating point summation 
#     
# $$
# \fl{x_1} ~\boxed{+}~ \fl{x_2} ~\boxed{+}~\cdots ~\boxed{+}~ \fl{x_n}
# $$ 
#     
# is an evaluation order for the $~\boxed{+}~$ operations. For example, $n = 4$, then the reduction $\Pi = (3,1,2)$ is corresponding to the following calculation
# 
# $$
# \fl{x_1} ~\boxed{+}~ \fl{x_2} + \underbrace{(\fl{x_3} ~\boxed{+}~ \fl{x_4})}_{=y_1}\to \underbrace{(\fl{x_1} ~\boxed{+}~ \fl{x_2})}_{=y_2} ~\boxed{+}~ y_1 \to \underbrace{(y_2 ~\boxed{+}~ y_1)}_{=y_3},
# $$
# 
# where $y_i$ denotes the result of $i$-th $\boxed{+}$ operation in the reduction. Obviously the final summation will be $y_{n-1}$. We denote $T_n(\Pi)$ be the result of floating point summation with reduction order $\Pi$.
# ````
# 
# (SSEC-CH1-EXER-THEORY)=
# ### Theoretical Part
# 
# ```{admonition} Problem 1
# Prove that the native summation has following error bound
# 
# $$
# |T_n(\Pi) - S_n|\le \left( \frac{\textrm{u} n }{1 - \textrm{u} n } \right)\sum_{j=1}^n |x_j|,
# $$
# 
# where $S_n = \sum_{j=1}^n x_j$ and $\Pi = (1,2,\dots, (n-1))$.
# ```
# 
# ```{admonition} Problem 2
# Prove that
# 
# $$
# \min_{\Pi}|T_n(\Pi) - S_n| \le \left(\frac{ \textup{u} H }{1 - \textup{u} H  }\right)\sum_{j=1}^n |x_j|,
# $$
# 
# where $H = \lceil\log_2 n \rceil$ and $T_n(\Pi)$ is the floating point summation with reduction order $\Pi$.
# ```
# 
# ```{admonition} Problem 3 (Horner's scheme)
# The evaluation of polynomial 
# 
# $$
# p(x) = a_0 + a_1 x + \dots + a_n x^n
# $$
# 
# is mostly using Horner's scheme, which writes the polynomial in 'nested' form:
# 
# $$
# p(x) = a_0 + x (a_1 + x(a_2 + \dots x(a_{n-1} + x(a_n)))).
# $$
#     
# Please find an upper bound of the rounding error for this scheme.         
# ```
# 
# (SSEC-CH1-EXER-COMP)=
# ###  Computational Part
# 
# ```{admonition} Problem 4 (Pairwise summation)
# Based on theoretical part, implement an algorithm for the summation $\sum_{j=1}^n x_j$ which has $\mathcal{O}(\textrm{u}\log_2 n)$ rounding error.
# ```
# 
# ````{admonition} Problem 5 (Kahan compensated summation)
# Suppose $a, b\in\mathbb{R}$, the rounding error for the sum $s = \fl{\fl{a}+\fl{b}}$ ($a\ge b$) can be computed using 
# 
# 
# $$
# \fl{\fl{s - \fl{a}}-\fl{b}}
# $$
# 
# Based on this property, one can keep tracking of the rounding error. 
# 
# ```{prf:algorithm} Kahan Compensated Summation
# :label: AL-KA-CO-SU
# 
# **Inputs** $\{x_j\}_{j=1}^n\subset \mathbb{R}$
# 
# **Outputs** $s_n = \sum_{j=1}^n x_j$ 
# 
# 1. $j\gets 1$, $e_j \gets 0$, $s_j \gets x_j$ //    initialization; 
# 2. While $j < n$:
# 
#     1. $j\gets j + 1$
#     2. $y_j = x_j - e_{j-1}$ //remove compensated error;
#     3. $s_j = s_{j-1} + y_j$ //perform summation;
#     4. $e_j = (s_j - s_{j-1}) - y_j$ //restore the rounding error;
# 
# ```
# Implement the {prf:ref}``AL-KA-CO-SU`` described above and compare the accuracy with naive summation and pairwise summation with test cases.
# ````
# 
# ```{admonition} Problem 6
# Suppose the inputs $\{x_j\}_{j=1}^n\subset \mathbb{R}$ are randomly distributed (say $x_j\sim U(0,1)$ i.i.d), what is growth of the expected rounding error with respect to total number of inputs $n$ for the naive summation and pairwise summation? Please provide an explanation of your result. You can use Kahan sum as the accurate result approximately. 
# ```

# ## Extended Reading
# 
# See {cite}`higham2019new,higham1993accuracy,muller2006elementary`.
# 
# 
# ```{bibliography}
# ```
