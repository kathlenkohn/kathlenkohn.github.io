--This is Macaulay2 code

loadPackage "SegreClasses"

--we work with 3x3 matrices
n=3

--we create the ring whose unkowns are the entries u_(i,j) of covariance matrices 
inds = flatten apply(n, i -> apply(i+1, j -> (i,j)))
R = QQ[apply(inds, I -> u_I)]

--a general covariance matrix:
Sigma = matrix apply(n, i -> apply(n, j -> if j > i then u_(j,i) else u_(i,j)))

--the following function computes the adjugate matrix of a matrix M
adj = (M) -> matrix apply(n, i -> apply(n, j -> (-1)^(i+j)*(determinant submatrix'(M,{i},{j}))))

--the following function compute the ideal of the reciprocal variety of an LSSM spanned by the matrices in a given list L
Linverse = (L) -> (
	Rtemp = QQ[apply(#L, i -> z_i)];
	Ltemp = apply(L, B -> sub(B,Rtemp));
	Mtemp = adj sum apply(#L, i -> z_i*(Ltemp#i));
	Ftemp = map(Rtemp, R, apply(inds, I -> Mtemp_I));
	kernel Ftemp
)

--the following function compute the ideal of the polar linear space of an LSSM spanned by the matrices in a given list L
Lpolar = (L) -> ideal apply(L, B -> trace(Sigma*B))

--the following code computes the ML degree of an LSSM spanned by the matrices in a given list L
mlDegree = (L) -> (
	s = apply(#inds, i -> random(QQ));
	subList = apply(#inds, i -> u_(inds#i) => s#i);
	S = sub(Sigma, subList);
	Crit = Linverse(L) + ideal apply(L, B -> trace((Sigma-S)*B));
	degree Crit
)


----------------------------------------------------
--Example 4.4

--the following is a basis for L
B1 = matrix{{1,0,0},{0,-1,0},{0,0,0}}
B2 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,1},{0,1,0}}
L = {B1,B2,B3,B4}

Linv = Linverse L
Lpol = Lpolar L

--we first compute the ring of P^5 x P^5
Rproduct = makeProductRing ({5,5})

--the first factor has the unknowns a,b,c,d,e,f
--the second factor has the unknowns g,h,i,j,k,l
F1 = map(Rproduct,R,{a,b,c,d,e,f})
F2 = map(Rproduct,R,{g,h,i,j,k,l})

--we embed the reciprocal variety in the first factor, and the polar space in the second factor
Linv = F1 Linv
Lpol = F2 Lpol

--we compute the ideal of the product (reciprocal variety) x (polar space)
V = Linv + Lpol

--we compute the ideal of the intersection as the diagonal of the above product
Diag = minors(2, matrix{{a,b,c,d,e,f},{g,h,i,j,k,l}})
Z = Linv + Lpol + Diag

--finally we obtain the Segre class in the Chow ring of P^5 x P^5
segre(Z,V)


----------------------------------------------------
--Example 4.5

--the following is a basis for L
B1 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B2 = matrix{{0,0,0},{0,0,1},{0,1,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,0},{0,0,1}}
L = {B1,B2,B3,B4}

--as in Example 4.4, we compute the Segre class in the Chow ring of P^5 x P^5
Linv = Linverse L
Lpol = Lpolar L
Rproduct = makeProductRing ({5,5})
F1 = map(Rproduct,R,{a,b,c,d,e,f})
F2 = map(Rproduct,R,{g,h,i,j,k,l})
Linv = F1 Linv
Lpol = F2 Lpol
V = Linv + Lpol
Diag = minors(2, matrix{{a,b,c,d,e,f},{g,h,i,j,k,l}})
Z = Linv + Lpol + Diag
segre(Z,V)
--the output is -7*H_1^5*H_2^5 + 2*H_1^5*H_2^4 + 2*H_1^4*H_2^5
--this corresponds to the Segre class -7 \zeta^5 + 2 \zeta^4 in the Chow ring of the diagonal, as described in Example 4.5 in the article


----------------------------------------------------
--Section 5.4
--we apply the function mlDegree defined above to one representative per congruence class
L = {B1,B2,B3,B4}
mlDegree L

--polar of [1 1 1]
B1 = matrix{{1,0,0},{0,-1,0},{0,0,1}}
B2 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,1},{0,1,0}}

--polar of [2 1]
B1 = matrix{{1,0,0},{0,-1,0},{0,0,0}}
B2 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,0},{0,0,1}}

--polar of [(1 1) 1]
B1 = matrix{{1,0,0},{0,-1,0},{0,0,0}}
B2 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,1},{0,1,0}}

--polar of [3]
B1 = matrix{{1,0,0},{0,0,0},{0,0,0}}
B2 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B3 = matrix{{0,0,1},{0,-2,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,0},{0,0,1}}

--polar of [(2 1)]
B1 = matrix{{1,0,0},{0,0,0},{0,0,0}}
B2 = matrix{{0,1,0},{1,0,0},{0,0,-2}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,1},{0,1,0}}

--polar of [; 1 ;]
B1 = matrix{{1,0,0},{0,0,0},{0,0,0}}
B2 = matrix{{0,0,0},{0,1,0},{0,0,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,0},{0,0,1}}

--polar of [1 1 ; ; 1]
B1 = matrix{{0,1,0},{1,0,0},{0,0,0}}
B2 = matrix{{0,0,0},{0,0,1},{0,1,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,0},{0,0,1}}

--polar of [2 ; ; 1]
B1 = matrix{{1,0,0},{0,0,0},{0,0,0}}
B2 = matrix{{0,0,0},{0,0,1},{0,1,0}}
B3 = matrix{{0,0,1},{0,0,0},{1,0,0}}
B4 = matrix{{0,0,0},{0,0,0},{0,0,1}}
