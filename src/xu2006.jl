module Xu

#=
This is an implementation of the algorithm presented in

"Alternative mathematical technique to determine LS spectral terms"

by Xu Renjun and Dai Zhenwen, published in JPhysB, 2006.
doi:10.1088/0953-4075/39/16/007

In the following, \(S'=2S\in\mathbb{Z}\) and
\(M_S'=2M_S\in\mathbb{Z}\), as in the original article.

\[\begin{aligned}
X(N,\ell,S',L) = A(N,\ell,\ell,S',L) - A(N,\ell,\ell,S',L+1)\\
+A(N,\ell,\ell,S'+2,L+1) - A(N,\ell,\ell,S'+2,L)
\end{aligned}
\quad(1)\]

with \(A(N,\ell,\ell_b,M_S',M_L)\) obeying four different cases:


Case 1, \(M_S'=1\), \(|M_L|\leq\ell\), and \(N=1\):

\[A(1,\ell,\ell_b,1,M_L) = 1\qquad(2)\]


Case 2, \(\{M_S'\}={2-N,4-N,...,N-2}\), \(|M_L| \leq
f\left(\frac{N-M_S'}{2}-1\right)+f\left(\frac{N+M_S'}{2}-1\right)\),
and \(1<N\leq 2\ell+1\):

\[\begin{aligned}
A(N,\ell,\ell,M_S',M_L) =
\sum_{\displaystyle M_{L-}=\max\left\{-f\left(\frac{N-M_S'}{2}-1\right),M_L-f\left(\frac{N+M_S'}{2}-1\right)\right\}}
^{\displaystyle\min\left\{f\left(\frac{N-M_S'}{2}-1\right),M_L+f\left(\frac{N+M_S'}{2}-1\right)\right\}}\\
\Bigg\{A\left(\frac{N-M_S'}{2},\ell,\ell,\frac{N-M_S'}{2},M_{L-}\right)
\times
A\left(\frac{N+M_S'}{2},\ell,\ell,\frac{N+M_S'}{2},M_L-M_{L-}\right)\Bigg\}
\end{aligned}\quad (3)\]


Case 3, \(M_S'=N\), \(|M_L|\leq f(N-1)\), and \(1<N\leq 2\ell+1\):

\[A(N,\ell,\ell_b,N,M_L) =
\sum_{\displaystyle M_{L_I} = \left\lfloor{\frac{M_L-1}{N}+\frac{N+1}{2}}\right\rfloor}
^{\displaystyle\min\{\ell_b,M_L+f(N-2)\}}
A(N-1,\ell,M_{L_I}-1,N-1,M_L-M_{L_I})
\]


Case 4, else:
\[A(N,\ell,\ell_b,M_S',M_L) = 0 \qquad (5)\]

\[f(n)=\begin{cases}
\displaystyle\sum_{m=0}^n(\ell-m), & n\geq0\\
0, & n<0
\end{cases}\]

=#

f(n::I,ℓ::I) where {I<:Integer} = n >= 0 ? sum(ℓ-m for m in 0:n) : 0

function A(N::I, ℓ::I, ℓ_b::I, M_S′::I, M_L::I) where {I<:Integer}
    if M_S′ == 1 && abs(M_L) <= ℓ && N == 1
        1 # Case 1
    elseif 1 < N && N <= 2ℓ + 1
        # Cases 2 and 3
        # N ± M_S' is always an even number
        a = (N-M_S′) >> 1
        b = (N+M_S′) >> 1
        fa = f(a-1,ℓ)
        fb = f(b-1,ℓ)

        if M_S′ in 2-N:2:N-2 && abs(M_L) <= fa + fb
            mapreduce(+, max(-fa, M_L - fb):min(fa, M_L + fb)) do M_Lm
                A(a,ℓ,ℓ,a,M_Lm)*A(b,ℓ,ℓ,b,M_L-M_Lm)
            end
        elseif M_S′ == N && abs(M_L) <= f(N-1,ℓ)
            mapreduce(+, floor(Int, (M_L-1)//N + (N+1)//2):min(ℓ_b,M_L + f(N-2,ℓ))) do M_LI
                A(N - 1, ℓ, M_LI - 1, N - 1, M_L-M_LI)
            end
        else
            0 # Case 4
        end
    else
        0 # Case 4
    end
end

X(N::I, ℓ::I, S′::I, L::I) where {I<:Integer} = A(N,ℓ,ℓ,S′,L) - A(N,ℓ,ℓ,S′,L+1) + A(N,ℓ,ℓ,S′+2,L+1) - A(N,ℓ,ℓ,S′+2,L)

end
