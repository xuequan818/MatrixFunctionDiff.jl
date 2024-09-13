var documenterSearchIndex = {"docs":
[{"location":"api/#Matrix-Functions-Differentiation-API","page":"API","title":"Matrix Functions Differentiation API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = MatrixFunctionDiff","category":"page"},{"location":"api/#Fréchet-derivatives-of-f(x::Real)::Number","page":"API","title":"Fréchet derivatives of f(x::Real)::Number","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"MatrixFunctionDiff.mat_fun_frechet","category":"page"},{"location":"api/#MatrixFunctionDiff.mat_fun_frechet","page":"API","title":"MatrixFunctionDiff.mat_fun_frechet","text":"mat_fun_frechet(f, eigs, Ψ::AbstractMatrix, h::Vector{AbstractMatrix})\nmat_fun_frechet(f, H::AbstractMatrix, h::Vector{AbstractMatrix})\n\nReturn the n-th order Fréchet derivative d^nf(H)h[1]…h[n], assuming f is called as f(x).\n\n\n\n\n\n","category":"function"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"The mat_fun_frechet function computes the higher order Fréchet derivatives rm d^nf(H)h1cdots hn. It requires that the matrix H can be diagonalizable, the function f be a scalar function, and h be a vector consisting of diagonalizable matrices. The order of the derivative is equal to the length of h.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"julia> using MatrixFunctionDiff\n\njulia> f(x) = 1.0 / (1.0 + exp(100 * (x - 0.1))); # Fermi Dirac function\n\njulia> order = 3; # 3rd order Fréchet derivative\n\njulia> N = 10; # dimension of matrix\n\njulia> X = rand(N, N);\n\njulia> H = 0.5 * (X + X');\n\njulia> h = [rand(N, N) for i = 1:order];\n\njulia> map!(x->0.5 * (x + x'), h, h);","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"MatrixFunctionDiff can also compute the Fréchet derivatives by the eigenpairs of H:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using LinearAlgebra\n\njulia> eigs, Ψ = eigen(H);\n\njulia> mat_fun_frechet(f, eigs, Ψ, h);","category":"page"},{"location":"background/#Fréchet-derivative-of-Matrix-functions","page":"Background","title":"Fréchet derivative of Matrix functions","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"Let mathcalH=mathbbC^Ntimes N_rm diag be the vector space of Ntimes N diagonalizable matrices. For HinmathcalH and h_1h_ninmathcalH, and f is an n times continuously differentiable function on a subset of mathbbC containing the spectrum of H+t_1h_1+cdots + t_nh_n, the n-th order Fréchet derivative of f(H) is","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    rm d^nf(H)h_1cdots h_n = frac12pi ioint_mathcalC f(z) sum_pinmathcalP_n(z-H)^-1h_p(1)(z-H)^-1cdots(z-H)^-1h_p(n)(z-H)^-1 dz\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"where mathcalC is a contour in the complex plane enclosing all the eigenvalues of H, and pinmathcalP_n is an arbitrary permutation of 1cdotsn. This can be proved by induction.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"For n = 1, we have ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    rm d f(H)h=lim_tto 0 fracf(H+th)-f(H)t=frac12pi ioint_mathcalC f(z)lim_tto 0frac(z-H-th)^-1-(z-H)^-1tdz\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Note that","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    (z-H-th)^-1  = (I-t(z-H)^-1h)^-1(z-H)^-1 = (I+t(z-H)^-1h+O(t^2))(z-H)^-1\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"we have ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    rm d f(H)h=frac12pi ioint_mathcalC f(z)(z-H)^-1h(z-H)^-1dz\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"which satisfies the formula.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Assume the n-1-th order derivative satisfies the formula. Then we have the n-th order derivative","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    rm d^nf(H)h_1cdots h_n =lim_tto 0 fracrm d^n-1f(H+th_n)h_1cdots h_n-1-rm d^n-1f(H)h_1cdots h_n-1t\n    =frac12pi ioint_mathcalC f(z)lim_tto 0sum_pinmathcalP_n-1frac1tBig((z-H-th_n)^-1h_p(1)cdots h_p(n-1)(z-H-th_n)^-1\n    qquad -(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1Big)dz\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Similarly, we have","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n   (z-H-th_n)^-1h_p(1)cdots h_p(n-1)(z-H-th_n)^-1\n    =(I+t(z-H)^-1h_n)(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1(I+th_n(z-H)^-1) + O(t^2)\n    =(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1 \n    quad+ tBig((z-H)^-1h_n(z-H)^-1h_p(1)cdots h_p(n-1)(z-H)^-1\n    qquadquad+(z-H)^-1h_p(1)(z-H)^-1h_ncdots h_p(n-1)(z-H)^-1\n    qquadquad+(z-H)^-1h_p(1)(z-H)^-1h_p(2)cdots h_n(z-H)^-1Big) + O(t^2)\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Therefore, we can obtain the formula.","category":"page"},{"location":"background/#Divided-difference-form","page":"Background","title":"Divided difference form","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"Let H=Phi LambdaPhi^-1=sum_i^Nlambda_iphi_iphi_i^-1, where phi_i is the i-th column of Phi and phi_i^-1 is the i-th row of Phi^-1, and assume there is no degeneration. Then for pinmathcalP_n we have","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    (z-H)^-1h_p(1)(z-H)^-1cdots(z-H)^-1h_p(n)(z-H)^-1\n    =sum_i_0cdotsi_n=1^Nphi_i_0(h_p(1))_i_0i_1cdots (h_p(n))_i_n-1i_nphi_i_n^-1(z-lambda_i_0)^-1cdots (z-lambda_i_n)^-1\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"where (h_p(k))_ij=phi_i^-1h_p(k)phi_j. We let","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"    (z-lambda_i_0)^-1cdots (z-lambda_i_n)^-1 = sum_k=0^nC_k (z-lambda_i_k)^-1","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"then","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"sum_k=0^nC_k prod_ellneq k(z-lambda_i_ell)=1","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Let z=lambda_i_k, we can obtain","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"C_k = frac1prod_ellneq k(lambda_i_k-lambda_i_ell)","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Therefore, we have","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginalign*\n    frac12pi ioint_mathcalC f(z)(z-lambda_i_0)^-1cdots (z-lambda_i_n)^-1dz = sum_k=0^nfracf(lambda_i_k)prod_ellneq k(lambda_i_k-lambda_i_ell)=flambda_i_0cdotslambda_i_n\nendalign*","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Finally, we obtain","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginequation\n    rm d^nf(H)h_1cdots h_n =sum_i_0cdotsi_n=1^Nphi_i_0Bigg(sum_pinmathcalP_n(h_p(1))_i_0i_1cdots (h_p(n))_i_n-1i_nBigg)flambda_i_0cdotslambda_i_nphi_i_n^-1\nendequation","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"which can be also written as ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"beginequation\n(Phi^-1rm d^nf(H)h_1cdots h_nPhi)_kell=sum_i_1cdotsi_n-1=1^NBigg(sum_pinmathcalP_n(h_p(1))_ki_1cdots (h_p(n))_i_n-1ellBigg)flambda_klambda_i_1cdotslambda_i_n-1lambda_ell\nendequation","category":"page"},{"location":"background/#Array-operations","page":"Background","title":"Array operations","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"Use array operations to efficiently compute the Fréchet derivative. For simplicity, just consider the no permutation case and define","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"(F_n)_kℓ=_i_1i_n-1=1^N(h_1)_ki_1 (h_n)_i_n-1ℓΛ^01n-1n_ki_1i_n-1ℓ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"where Λ^0n_i_0i_n = fλ_i_0λ_i_n. It is immediately to obtain that ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"F_1 =  h_1  Λ^01","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"and ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"F_2 = _i=1^N (mathfrakh^12  Λ^012)_i","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"with mathfrakh^12_i = (h_1)_i(h_2)_i. For n  3, first compute ","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"mathfrakF^02n_j_3j_n = _i=1^N (mathfrakh^12  Λ^01n_j_3j_n)_i","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"and permute such that the 0-dimension is at the end mathfrakF^2n0. Then for m  3, there is the recursion","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"mathfrakF^mn0_j_mj_n = _i=1^N (h_m  mathfrakF^m-1n0_j_mj_n)_i","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"and F_n = (mathfrakF^n0)^T.","category":"page"},{"location":"#MatrixFunctionDiff.jl","page":"Home","title":"MatrixFunctionDiff.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A julia package for computing Fréchet derivatives of scalar functions with matricies as variables. The higher order Fréchet derivatives are first written in a formula similar to the Daleskii-Krein theorem, and then computed by the divided difference.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MatrixFunctionDiff.jl is a registered package, so it can be installed by running","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> Pkg.add(\"MatrixFunctionDiff\")","category":"page"},{"location":"#Related-packages","page":"Home","title":"Related packages","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"DividedDifferences.jl: Divided difference for Julia","category":"page"}]
}
