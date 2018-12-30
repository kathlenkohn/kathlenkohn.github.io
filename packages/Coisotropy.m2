newPackage(
    	"Coisotropy",
	Version => "0.1",
	Date => "June 2016",
	Authors => {
	    	{Name => "Kathlen Kohn", Email => "", HomePage=>""}
		},
	Headline => "Computes Chow forms, Hurwitz forms, etc., and recovers underlying varieties",
	Configuration => {
		"path" => "",
		"fig2devpath" => "",
		"keepfiles" => false
	},
        PackageExports => {},
	DebuggingMode => true
)

export {
  "primalToDual",
  "dualToPrimal",
  "isCoisotropic",
  "recoverVar",
  "dualVariety",
  "coisotropicForm",
  "polarDegrees",
  "Smooth",
  "PolarDegrees"
}


------------------------------------------------------------------------------
-- CODE
------------------------------------------------------------------------------
primalToDual = method(TypicalValue => Thing)
--computes polynomial in dual Plücker coordinates, given a polynomial in primal Plücker coordinates

primalToDual (Thing,ZZ,ZZ) := (Q,k,n) -> (
-- first, bring given polynomial into the right format
	I := ideal(Q);
	J:=parsePolynomial(I,n-k-1,n);
--check for errors returned by parse function
	if (instance(J,String)) then return J;
	var := flatten entries vars ring J;
--compute new ring;
	p := local p;
	R := QQ[apply(sort subsets(n+1,k+1), s -> p_s)];
	eq := (mingens J)_(0,0);
--compute new variables with plus minus
	sets := sort subsets(n+1,n-k);
	complements := apply(sets, s -> setComplement(s,n));
	pl := apply(#sets, i -> p_(complements#i)*sgn(sets#i,complements#i));
--substitute: old variables => new variables
	subst := apply(#var, i -> var#i => pl#i); 
	eq = sub(eq, subst);
	return sub(eq, R);
)

dualToPrimal = method(TypicalValue => Thing)
--computes polynomial in primal Plücker coordinates, given a polynomial in dual Plücker coordinates

dualToPrimal (Thing,ZZ,ZZ) := (Q,k,n) -> (
-- first, bring given polynomial into the right format
	I := ideal(Q);
	J:=parsePolynomial(I,k,n);
--check for errors returned by parse function
	if (instance(J,String)) then return J;
	var := flatten entries vars ring J;
--compute new ring;
	p := local p;
	R := QQ[apply(sort subsets(n+1,n-k), s -> p_s)];
	eq := (mingens J)_(0,0);
--compute new variables with plus minus
	sets := sort subsets(n+1,k+1);
	complements := apply(sets, s -> setComplement(s,n));
	pl := apply(#sets, i -> p_(complements#i)*sgn(complements#i,sets#i));
--substitute: old variables => new variables
	subst := apply(#var, i -> var#i => pl#i); 
	eq = sub(eq, subst);
	return sub(eq, R);
)

setComplement = (s,n) ->(
--computes {0, ..., n}\s
	k := #s;
	t := apply(n+1, i -> i);
	scan(s, e -> t = delete(e,t));
	scan(subsets(n+1, n-k+1), u -> if (u==t) then t=u);
	return t
)

sgn = (s,t) -> (
--computes the sign of the permutation (s,t)
	perm := s|t;
--compute first how many elements to the left of the i-th element are greater than the i-th element
	inversionVector := apply(#perm, i -> number(take(perm,i),m -> m > perm#i));
--then sum these up and check parity
	length := sum inversionVector;
	if (even length) then (return 1) else (return -1);
)

--TODO: better method?
iteratedColons = (I, saturateBy, elimVar) -> (
	J := I;
	saturateIdeal := ideal saturateBy;
	local saturateByEq;
	elimIdeal := eliminate(elimVar, J);
	while (numgens(elimIdeal)==0 and #saturateBy>0) do {
		saturateByEq = first saturateBy;
		saturateBy = drop(saturateBy,1);
		print ("compute colon by "|toString(saturateByEq));
		J = J:ideal(saturateByEq);
		elimIdeal = eliminate(elimVar, J);
	};
	if (numgens(elimIdeal) == 0) then (
		print "have to compute saturation!";
		J = saturate(I, saturateIdeal);
		elimIdeal = eliminate(elimVar, J);
	);
	return elimIdeal;
)


chowForm = method(TypicalValue => Ideal)
--computes the Chow form of a given ideal

chowForm (Ideal) := (I)->(
--compute new ring
	R := ring I;
	n := (numgens R)-1;
	k := dim I-1;
	x := new MutableList;
	scan(n+1,i -> x#i = (vars R)_(0,i));
	a := local a;
	R = R[a_(0,0)..a_(n,k)];
	(S,M) := flattenRing R;
--compute linear forms that cut out linear subspace
	l := apply(k+1,j -> sum(0..n, i -> a_(i,j)*x#i));
--add these linear forms to given ideal
	J := substitute(I + ideal l,S);

--saturate by irrevalent ideal (s.t. all-zero-point does not count), and eliminate old variables
--for efficiency, compute colon ideals elementwise und check if ideal still becomes empty after elimination	
	scan(n+1,i -> x#i = substitute(x#i,S));
	x = toList x;
	dual := iteratedColons(J, x, x);
--check if resulting ideal is principal and of correct degree
	ch := flatten entries mingens dual;
	if (#ch > 1) or ((degree(ch#0))#0 != degree(I)*#l)then (
		print "have to compute saturation!";
		J = saturate(J,ideal x);
		elimIdeal := eliminate(x, J);
		ch = flatten entries mingens elimIdeal;
	);
--compute new ring for ch
	R = QQ[a_(0,0).. a_(n,k), MonomialOrder=>{Lex}];
	ch=sub(ch#0,R);
--this gives chow form in primal Stiefel coordinates; compute now primal Pluecker coordinates
	return toPluecker(ch,n,k)
)    

isCoisotropic = method(TypicalValue => Boolean)
--checks if a hypersurface, given by a principal radical ideal, is coisotropic

isCoisotropic (Thing, ZZ, ZZ) := (CH, k, n)->(
	I := ideal CH;
--check if CH is squarefree (this is equivalent to that the ideal is radical)
	if (not isSquarefree(I)) then return "given polynomial not squarefree!";
--bring I into right format
	I=parsePolynomial(I,n-k-1,n);
--check for errors coming from parse function
	if (instance(I,String)) then return I;
--define new ring in affine coordinates
	a := local a;
	R := QQ[a_(0,n-k) .. a_(n-k-1,n)];
	P := (mingens I)_(0,0);
--compute minors of standard affine chart
	A := id_(R^(n-k)) | transpose genericMatrix(R,k+1,n-k);
	pl := apply(sort subsets(n+1,n-k), s -> determinant submatrix(A,s));
--evaluate given polynomial at minors
	var := flatten entries vars ring P;
	P = sub(P, apply(#pl, i -> var#i => pl#i));
--compute Jacobian matrix
	J := matrix apply(n-k, i -> apply(k+1, j -> diff(a_(i,j+n-k),P)));
--check if Jacobian matrix has rank at most one everywhere on the given hypersurface
	return isSubset(minors(2,J),ideal(P))
)

recoverVar = method(TypicalValue => Ideal)
--recovers the underlying projective variety of a coisotropic hypersurface

recoverVar (Thing, ZZ, ZZ) := (CH,k,n) -> (
	if (# gens ring CH != #subsets(n+1,k+1)) then return "number of variables in ring of input polynomial does not equal (n+1) choose (k+1)";
	I := ideal CH;
--for efficiency: check if k < n-k-1
	local useDualCoords;
	local l;
	if (k < n-k-1) then (useDualCoords = true; l = k) else (useDualCoords = false; l=n-k-1);
--bring I into right format
	local J;	
	if (useDualCoords) then (
		J = ideal primalToDual(CH,k,n);
	) else (
		J = parsePolynomial(I,n-k-1,n);
	);
--check for errors coming from parse function
	if (instance(J,String)) then return J;
	a := local a;
--define new ring in Stiefel coordinates
	R := QQ[a_(0,0) .. a_(l,n)];
	P := (mingens J)_(0,0);
--compute minors of Stiefel matrix
	A := transpose genericMatrix(R,n+1,l+1);
	pl := apply(sort subsets(n+1,l+1), s -> determinant submatrix(A,s));
--evaluate given polynomial at minors
	var := flatten entries vars ring P;
	P = sub(P, apply(#pl, i -> var#i => pl#i));
--since given hypersurface must be coisotropic, its dual (in Stiefel coordinates) factors as X x P^l
	X := dualVarietyFactored(ideal P, l, n, false, true);

--check if i-th coisotropic form of X is CH
	print "check result";
	i := dim(X)-1-l;
	ch := coisotropicForm(X,i);
	I = parsePolynomial(ideal ch,l,n);
	f := map(ring I, ring J, vars ring I);
	if not (I == f J) then (
		print "have to compute saturation";		
		X = dualVarietyFactored(ideal P, l, n, false, false);
	);
--if dual coordinates were used, then the current result is the dual of the required result
	if (useDualCoords) then (X = dualVariety X);
	return X;
)

polarDegrees = method(TypicalValue => List, Options => {
	Smooth => false
	})
--computes polar degrees as multidegree of conormal variety

polarDegrees (Ideal) := opt -> (I) -> (
	R := ring I;
	oldVar := gens R;
	n := #oldVar-1;
	cod := codim I;
--compute large ring
	v := local v;
	R = R[v_0 .. v_n];
	(S,D) := flattenRing R;
	var := gens R;
--compute hyperplanes that are tangent: in fact, compute {(x,H) | H is tangent to x}
	A := (jacobian I) | (matrix apply(#var, i -> {var#i}));
	Con := sub(I + minors(cod+1,A),S);
--use graded ring to compute multidegree of Con
	degreeList := (apply(#var, i -> {1,0})|apply(#var, i -> {0,1}));
	R = QQ[oldVar|var, Degrees => degreeList];
	Con = sub(Con,R);
--check if variety defined by I is mooth
	isSmooth := opt.Smooth;
	if (not isSmooth) then isSmooth = (0 == dim singularLocus I);
	mult := local mult;
	if (isSmooth) then (
		mult = multidegree Con;
--we just need to remove monomials of the form t^(n+1) from mult for the extraneous component {0}xP^n
		t := gens ring mult;
		T := QQ[t];
		mult = sub(mult,T);
		mult = mult % sub((t#0)^(n+1),T);
		mult = mult % sub((t#1)^(n+1),T);
	) else (
		Sing := ideal singularLocus I;
		Sing = sub(Sing, R);
		Con = saturate(Con,Sing);
		mult = multidegree Con
	);
	(M,C) := coefficients mult;
	mu := numgens target C;                        --this computes dim(X)-codim(X^v)+2
	return apply(mu, i -> sub(C_(mu-i-1,0),ZZ))
)

dualVariety = method(TypicalValue => Ideal, Options => {
	Smooth => false
	})
--computes the projectively dual variety of a given variety

dualVariety (Ideal) := opt -> (I) -> (
--this is just the special case where the dual variety factors as X x P^0
	isSmooth := opt.Smooth;
	if (not isSmooth) then isSmooth = (0 == dim singularLocus I);
	D := dualVarietyFactored(I, 0, (numgens ring I)-1, isSmooth, true);
--check biduality theorem
	print "check result";
	isSmoothD := (0 == dim singularLocus D);
	J := dualVarietyFactored(D, 0, (numgens ring D)-1, isSmoothD, true);
	f := map(ring I, ring J, vars ring I);
	if not (I == f J) then (
		print "have to compute saturation!";
		D = dualVarietyFactored(I, 0, (numgens ring I)-1, isSmooth, false);
	);
	return D;
)

dualVarietyFactored = (I,k,n,isSmooth,useIteratedColons) -> (
--computes the projectively dual variety of a given variety, when one knows that the dual factors as X x P^k
--idea: introduce new variables for X, but choose fixed prime numbers for P^k instead of more variables 
--first compute basic parameters
	R := ring I;
	oldVar := gens R;
	size := (k+1)*(n+1);
	cod := codim I;
--compute large ring
	v := local v;
	R = R[v_0 .. v_n];
	(S,M) := flattenRing R;
--find k prime numbers
	primes := {};
	integ := 11;
	while (#primes < k) do{
		if (isPrime(integ)) then primes = primes|{integ};
		integ = integ+1;
	};
	primes = {1}|primes;
--multiply variables for X with found primes
	var := gens R;
	var = flatten apply(k+1, i -> (primes#i)*var);
--compute hyperplanes that are tangent: in fact, compute {(x,H) | H is tangent to x}
	A := (jacobian I) | (matrix apply(size, i -> {var#i}));
	Con := sub(I + minors(cod+1,A),S);
	oldVar = apply(oldVar, v -> sub(v,S));
	saturateByIdeal := local saturateByIdeal;
	elimIdeal := local elimIdeal;
	if isSmooth then saturateByIdeal = ideal oldVar else saturateByIdeal = sub(minors(cod, jacobian I), S);
	if useIteratedColons then (
		elimIdeal = iteratedColons(Con, flatten entries mingens saturateByIdeal, oldVar)
	) else (
		Con = saturate(Con, saturateByIdeal);
		elimIdeal = eliminate(oldVar, Con)
	);
--compute dual variety in new ring
	R = QQ[gens R];
	return sub(elimIdeal,R)
)

coisotropicForm = method(TypicalValue => Thing, Options => {
	Smooth => false,
	PolarDegrees => {}
	})
--computes the i-th coisotropic form of a given variety

coisotropicForm (Ideal, ZZ) := opt -> (I, i)->(
	PD := opt.PolarDegrees;
	if not instance(PD,List) then return "optional input Polar has to be a list!";
	if (i<0) then return "i has to be a positive integer!";
	if (i==0) then return chowForm I;
--check if i is in the range s.t. the i-coisotropic variety is a hypersurface	
	if (#PD == 0) then  PD = polarDegrees(I,Smooth => opt.Smooth);
	if (not i<#PD) then (
	numberString := i|"-th";
	if (i==1) then numberString = "first"
	else if (i==2) then numberString = "second"
	else if (i==3) then numberString = "third";
	return "the "|numberString|" coisotropic variety is not a hypersurface!";
	);
--compute new ring
	R := ring I;
	n := (numgens R)-1;
	k := dim I-1;
	x := new MutableList;
	scan(n+1,i -> x#i = (vars R)_(0,i));
	a := local a;
	R = R[a_(0,0)..a_(n,k-i)];
	(S,M) := flattenRing R;
--compute linear forms that cut out linear subspace
	l := apply(k+1-i,j -> sum(0..n, i -> a_(i,j)*x#i));
--add these linear forms to given ideal
	J := I + ideal l;
--add the condition that dim(TxX \cap H) \geq i
	M = (jacobian I)|(transpose genericMatrix(R,k+1-i,n+1));
	J = J+ minors(n+1-i,(jacobian I)|(transpose genericMatrix(R,k+1-i,n+1)));
	J = sub(J,S);

--saturate by singular locus, and eliminate old variables
--for efficiency, compute colon ideals elementwise und check if ideal still becomes empty after elimination	
	scan(n+1,i -> x#i = substitute(x#i,S));
	x = toList x;
	isSmooth := opt.Smooth;
	if (not isSmooth) then isSmooth = (0 == dim singularLocus I);
	local saturateBy;
	if isSmooth then saturateBy = x else saturateBy = flatten entries mingens sub(minors(codim I, jacobian I), S);
	elimIdeal := iteratedColons(J, saturateBy, x);

--check if resulting ideal is principal and of correct degree
	ch := flatten entries mingens elimIdeal;
	if (#ch > 1) or ((degree(ch#0))#0 != (PD#i)*#l)then (
		print "have to compute saturation!";
		J = saturate(J,ideal saturateBy);
		elimIdeal = eliminate(x, J);
		ch = flatten entries mingens elimIdeal;
	);
--compute new ring for ch
	R = QQ[a_(0,0).. a_(n,k-i), MonomialOrder=>{Lex}];
	ch=sub(ch#0,R);
--this gives chow form in primal Stiefel coordinates; compute now primal Pluecker coordinates
	return toPluecker(ch,n,k-i);
)

parsePolynomial = (I,k,n) -> (
	--helping function: takes given polyomials and writes them in primal Plücker coordinates of G(k+1,n+1), where the Plücker variables are of the form p_{0, 1, ...}
	R := ring I;
--sort variables of R first!
	var := flatten entries vars R;
	varString := apply(var, v -> toString v);
	varSorted := sort varString;
	oldVar := new MutableList;
	scan(#var, i -> oldVar#i = var#(position(varString, v -> v==(varSorted#i))));
	p := local p;
--compute new variables
	newVar := sort apply(subsets(n+1,k+1), s -> p_s);
	if (#var != #newVar) then return "number of variables in ring does not equal (n+1) choose (k+1)";	
--compute new ring
	R = R[newVar];
	(S,M) := flattenRing R;
--recompute newVar, s.t. it is known in S
	newVar = flatten entries vars R;
--add to I all equations (newVar = var) and eliminate old variables
	J:= ideal apply(#var, i -> oldVar#i-newVar#i);
	J = J + I;
	scan(#var, i -> oldVar#i = sub(oldVar#i,S));
	J = sub(J,S);
	J = eliminate(toList oldVar,J);
	R = QQ[newVar];
	return sub(J,R)
)

isSquarefree = (I) ->(
--checks if the generating polynomial of a given principal ideal is squarefree
	return isSubset(radical(I),I)
)

toPluecker = (chA,n,k)->(
--Bernds algorithm 3.2.8: compute polynomial in Plücker coordinates from a given polynomial in Stiefel coordinates
--compute new ring	
	p := local p;
	allSubsets := subsets(n+1,k+1);
	R := QQ[apply(allSubsets, s -> p_s)];
--compute Stiefel matrix
	M := genericMatrix(ring chA, k+1, n+1);
--1) if chA==0, then summand is 0 and algorithm is done
	chP := 0_R;
	local d;
	local T;
	local coef;
	local s;
	local plueckerList;
	while(chA!=0) do {
--2) m = init(chA)
		m := leadMonomial(chA);		
--3) compute diagonal form of m
--determine size of tableaux T: d x (k+1)
		d = (degree m)#0%(k+1);
		if (d!=0) then return "for some monomial: degree not divisible by (dim X + 1)";
		d = (degree m)#0//(k+1);
		T = diagonalForm(m,k);
--4) summand is coefficient(m)*T, and replace chA
--compute list of Pluecker coordinates
		plueckerList={};
		scan(d, i -> (
			s = apply(k+1, j -> T#j#i);
			scan(allSubsets, t -> if (t==s) then plueckerList = plueckerList|{t});
		));
--add new summand
		coef=leadCoefficient(chA);
		chP=chP + coef*product(apply(plueckerList, t -> p_t));
--replace chA
		chA = chA - coef*product(apply(plueckerList, t -> determinant submatrix(M,t)));
	};
	return chP
)

diagonalForm = (m,k) ->(
--initialize tableaux
	T := new MutableList;
	scan(k+1, j -> T#j = new List);
--if the support of m is empty, algorithm is done
	supp := support m;
	local var;
	local i;
	local j;
	local ind;
	while((#supp)>0) do {
--take a variable occuring in m and determine its indices (i,j): index var = i*(k+1)+j
		var = supp#0;
		ind = index var;
		j = ind % (k+1);
		i = (ind-j)/(k+1);
--add i into the j-th list in T
		T#j = T#j|{i};
--divide m by the considered variable and update support
		m = m/var;
		supp = support m;
	};
	scan(k+1, j -> T#j = sort T#j);
	return T;
)




------------------------------------------------------------------------------
-- DOCUMENTATION
------------------------------------------------------------------------------
beginDocumentation()
doc ///
    Key
    	Coisotropy
    Headline
    	a M2 package for coisotropic hypersurfaces in Grassmannians
    Description
    	Text
	    This package determines if a given hypersurface in the Grassmannian is coisotropic and can recover its underlying variety. Given a variety, it can compute all its coisotropic hypersurfaces.
///

doc ///
    Key
	polarDegrees
	(polarDegrees, Ideal)
	[polarDegrees, Smooth]
    Headline
    	computes the degrees of all coisotropic hypersurfaces of a given projective variety
    Usage
    	polarDegrees(I)
	polarDegrees(I, Smooth => false)
    Inputs
    	I:Ideal 
	Smooth=>Boolean
    Outputs
    	L:List
    Description
	Text
	    This function computes a list whose i-th entry is the degree of the i-th coisotropic hypersurface of the projective variety given by the ideal I. These degrees are exacty the well-studied polar degrees. If you know that the variety given by I is smooth, you can use the optional input Smooth=>true to speed up the computation.
	Example
	    R = QQ[x_0..x_3];
	    I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	    polarDegrees(I)
	    polarDegrees(I, Smooth => true)
///

doc ///
    Key
	isCoisotropic
	(isCoisotropic, Thing, ZZ, ZZ)
    Headline
    	checks if a hypersurface in a Grassmannian, given by its defining polynomial, is coisotropic
    Usage
    	isCoisotropic(CH,k,n)
    Inputs
    	CH:Thing
    	k:ZZ 
    	n:ZZ 
    Outputs
    	b:Boolean
    Description
	Text
	    This function checks if a hypersurface in G(k+1,n+1) = {(k+1)-dimensional linear subspaces of C^(n+1)} is coisotropic. The hypersurface has to be given by its defining polynomial CH in primal Plücker coordinates p_{i_0, ..., i_(n-k-1)}. If you have your polynomial in dual Plücker coordinates q_{j_0, ..., j_k}, you can use the method dualToPrimal(CH,k,n) before calling isCoisotropic to parse your polynomial. The parameters (k,n) determine the used Grassmannian G(k+1,n+1). In particular, the number of variables in the ring of CH has to be equal to (n+1) choose (k+1).
      	Example
	    R = QQ[p01,p02,p03,p12,p13,p23];
	    isCoisotropic(p12^2+p03^2-2*p02*p13-2*p01*p23, 1, 3)
///

doc ///
    Key
	dualToPrimal
	(dualToPrimal, Thing, ZZ, ZZ)
    Headline
    	computes a polynomial in primal Plücker coordinates, given a polynomial in dual Plücker coordinates
    Usage
    	dualToPrimal(Q,k,n)
    Inputs
    	Q:Thing
    	k:ZZ 
    	n:ZZ 
    Outputs
    	P:Thing
    Description
	Text
	    This function transforms the given polynomial Q in dual Plücker coordinates q_{j_0, ..., j_k} of G(k+1,n+1) = {(k+1)-dimensional linear subspaces of C^(n+1)} to a polynomial P in primal Plücker coordinates p_{i_0, ..., i_(n-k-1)} of G(k+1,n+1). The parameters (k,n) determine the used Grassmannian G(k+1,n+1). In particular, the number of variables in the ring of Q has to be equal to (n+1) choose (k+1).
      	Example
	    R = QQ[p01,p02,p03,p12,p13,p23];
	    dualToPrimal(p01-p02+p03+p12-p13+p23, 1, 3)
///

doc ///
    Key
	primalToDual
	(primalToDual, Thing, ZZ, ZZ)
    Headline
    	computes a polynomial in dual Plücker coordinates, given a polynomial in primal Plücker coordinates
    Usage
    	primalToDual(P,k,n)
    Inputs
    	P:Thing 
    	k:ZZ 
    	n:ZZ 
    Outputs
    	Q:Thing
    Description
	Text
	    This function transforms the given polynomial P in primal Plücker coordinates p_{i_0, ..., i_(n-k-1)} of G(k+1,n+1) = {(k+1)-dimensional linear subspaces of C^(n+1)} to a polynomial Q in dual Plücker coordinates q_{j_0, ..., j_k} of G(k+1,n+1). The parameters (k,n) determine the used Grassmannian G(k+1,n+1). In particular, the number of variables in the ring of P has to be equal to (n+1) choose (k+1).
      	Example
	    R = QQ[p01,p02,p03,p12,p13,p23];
	    primalToDual(p01+p02+p03+p12+p13+p23, 1, 3)
/// 

doc ///
    Key
	recoverVar
	(recoverVar, Thing, ZZ, ZZ)
    Headline
    	recovers the underlying projective variety of a given coisotropic hypersurface
    Usage
    	recoverVar(CH,k,n)
    Inputs
    	CH:Thing 
    	k:ZZ 
    	n:ZZ 
    Outputs
    	I:Ideal
    Description
	Text
	    This function computes the ideal I of the underlying projective variety of a given coisotropic hypersurface in G(k+1,n+1) = {(k+1)-dimensional linear subspaces of C^(n+1)}. The hypersurface has to be given by its defining polynomial CH in primal Plücker coordinates p_{i_0, ..., i_(n-k-1)}. If you have your polynomial in dual Plücker coordinates q_{j_0, ..., j_k}, you can use the method dualToPrimal(CH,k,n) before calling recoverVar to parse your polynomial. Moreover, you have to know that the given hypersurface is coisotropic, since this function does not check it (due to efficiency). You can manually check if the given hypersurface is coisotropic by calling isCoisotropic(CH,k,n) before calling recoverVar. The parameters (k,n) determine the used Grassmannian G(k+1,n+1). In particular, the number of variables in the ring of CH has to be equal to (n+1) choose (k+1).
      	Example
	    R = QQ[p01,p02,p03,p12,p13,p23];
	    recoverVar(p12^2+p03^2-2*p02*p13-2*p01*p23, 1, 3)
///

doc ///
    Key
	coisotropicForm
	(coisotropicForm, Ideal, ZZ)
	[coisotropicForm, Smooth]
	[coisotropicForm, PolarDegrees]
    Headline
    	computes the i-th coisotropic form of a given projective variety
    Usage
    	coisotropicForm(I,i)
    	coisotropicForm(I,i, Smooth => false)
    	coisotropicForm(I,i, PolarDegrees => {})
    Inputs
    	I:Ideal 
    	i:ZZ 
	Smooth=>Boolean
	PolarDegrees=>List
    Outputs
    	CH:Thing
    Description
	Text
	    This function computes the i-th coisotropic form CH of the projective variety given by the ideal I. If the outcome would not be a hypersurface, then a message is returned. If the outcome is a hypersurface, the polynomial CH is written in primal Plücker coordinates p_{i_0, ..., i_(n-k-1)} of G(k+1,n+1) = {(k+1)-dimensional linear subspaces of C^(n+1)}. If you know that the variety given by I is smooth, you can use the optional input Smooth=>true to speed up the computation. Moreover, if you already have the list L of polar degrees of the variety given by I, you can use the optional input PolarDegrees=>L such that the method coisotropicForm does not need to compute the polar degrees again. 
      	Example
	    R = QQ[a..d];
	    I = ideal(a*d-b*c);
	    coisotropicForm(I, 1)	    
	    coisotropicForm(I, 1, Smooth => true)	    
	    coisotropicForm(I, 1, PolarDegrees => {2,2,2})	    
	    coisotropicForm(I, 1, Smooth => true, PolarDegrees => {2,2,2})	    
/// 

doc ///
    Key
	dualVariety
	(dualVariety, Ideal)
	[dualVariety, Smooth]
    Headline
    	computes the projectively dual variety of a given projective variety
    Usage
    	dualVariety(I)
    	dualVariety(I, Smooth => false)
    Inputs
    	I:Ideal 
	Smooth=>Boolean
    Outputs
    	J:Ideal
    Description
	Text
	    This function computes the ideal J of the projectively dual variety of the projective variety given by the ideal I. If you know that the variety given by I is smooth, you can use the optional input Smooth=>true to speed up the computation.
      	Example
	    R = QQ[a..d];
	    dualVariety ideal(a*d-b*c)	    
/// 

doc///
    Key
        Smooth
    Headline
                option to ensure the input variety is smooth
    Usage
        polarDegrees(I,Smooth=>false)
        dualVariety(I,Smooth=>false)
        coisotropicForm(I,i,Smooth=>false)
    Description
                Text
			This might speed up the computation since the methods do not have to check for smoothness of the input. The best speed up is expected for the method polarDegrees.
                Example
	    		R = QQ[x_0..x_3];
	    		I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	    		polarDegrees(I, Smooth => true)
///


doc///
    Key
        PolarDegrees
    Headline
                option to include the polar degrees of the input variety
    Usage
        coisotropicForm(I,i,PolarDegrees=>{})
    Description
                Text
			This might speed up the computation since the method coisotropicForm does not have to re-compute the polar degrees when they are already known.
                Example
	    		R = QQ[a..d];
	    		I = ideal(a*d-b*c);
	    		coisotropicForm(I, 1, PolarDegrees => {2,2,2})	    
///




--------------------- TESTS ----------------------------------


--primalToDual and dualToPrimal
TEST ///
	R = QQ[x_0..x_3];
	I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	ch = coisotropicForm(I,0, Smooth => true, PolarDegrees => {3,4});
	chDual = primalToDual(ch,1,3);
	ch2 = dualToPrimal(chDual,1,3);
	G = gens ring ch;
	f = map(ring ch, ring ch2, matrix{{G#0, G#1, G#3, G#2, G#4, G#5}});
	assert(ideal ch == ideal f ch2)
///
 


--isCoisotropic
TEST ///
	G = Grassmannian (1,3, CoefficientRing => QQ);
	P = gens ring G;
	assert(isCoisotropic (P#3,1,3) == true)
	assert(isCoisotropic (P#3-P#2,1,3) == false)

	R = QQ[x_0..x_3];
	I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	ch = coisotropicForm(I,0, Smooth => true, PolarDegrees => {3,4});
	assert(isCoisotropic(ch,1,3))
///

 
--dualVariety
TEST ///
	R = QQ[x_0..x_3];
	I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	D= dualVariety(I, Smooth => true);
	I2 = dualVariety(D);
	f = map(R, ring I2, vars R);
	assert(I == f I2)

	assert(polarDegrees(D) == reverse polarDegrees(I))

	ch = coisotropicForm(I,0, Smooth => true);
	hur = coisotropicForm(D,1);
	chDual = primalToDual(ch,1,3);
	G = gens ring chDual;
	f = map(ring chDual, ring hur, matrix{{G#0, G#1, G#3, G#2, G#4, G#5}});
	assert(ideal chDual == ideal f hur)
///

--coisotropicForm and recoverVar
TEST ///
	R = QQ[x_0..x_3];
	I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	ch1 = coisotropicForm(I,0, Smooth => true, PolarDegrees => {3,4});
	ch2 = coisotropicForm(I,0);
	f = map(ring ch1, ring ch2, vars ring ch1);
	assert(ch1 == f ch2)

	I2 = recoverVar(ch1,1,3);
	f = map(R, ring I2, vars R);
	assert(I == f I2)
///

--polarDegrees
TEST ///
	R = QQ[x_0..x_3];
	I = ideal(x_1^2-x_0*x_2, x_2^2-x_1*x_3, x_1*x_2-x_0*x_3);
	assert(polarDegrees(I, Smooth => true) == {3,4})
///
