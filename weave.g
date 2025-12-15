
#############################################################################
##
##  Compute the action of the mxmxm cube group on the m^3 cube-positions
##  as permutations.
## 
##  The main functions we need are to generate G := R:T:
##  CubeRotationMatrixGroup() yields the cube group as a matrix group 
##  CubeRotationGroup(m)   computes  R
##  CubeTranslationGroup(m)  computes T
##  CubeFullSymmetryGroup(m) computes G = R:T
##  FindCycleLengths( m )  computes cycle lengths of sigma T for
##  sigma the representatives of the conjugacy classes of R
##

# encode (x,y,z) with 0 <= x,y,z < m into [1..m^3]
CubeEncode := function(x,y,z,m)
    return 1 + x + m*y + m^2*z;
end;

# decode index n (1..m^3) into x,y,z in 0..m-1
CubeDecode := function(n,m)
    local k, z, y, x;
    k := n - 1;
    z := QuoInt(k, m^2);
    k := RemInt(k, m^2);
    y := QuoInt(k, m);
    x := RemInt(k, m);
    return [x, y, z];
end;

ApplyMat := function(A, v)
    local i, j, res;
    res := [ 0, 0, 0 ];
    for i in [1..3] do
        res[i] := 0;
        for j in [1..3] do
            res[i] := res[i] + A[i][j] * v[j];
        od;
    od;
    return res;
end;

ApplyMatTrans := function(A,tau, v, m)
        local res, centre, poslist, v0, nx, ny, nz, rotated;

        centre := (m-1)/2;
        poslist := List([0..m-1], i -> i - centre);
        v0 := [ v[1] - centre, v[2] - centre, v[3] - centre ];

        rotated := ApplyMat(A,v0);

        # find positions (should match exactly one element in poslist)
        nx := Position(poslist, rotated[1]);
        ny := Position(poslist, rotated[2]);
        nz := Position(poslist, rotated[3]);

        if nx = fail or ny = fail or nz = fail then
            Error("Rotation produced coordinate not on grid; m = ", m,
                  " rotated = ", rotated, " poslist = ", poslist);
        fi;

        return [nx-1, ny-1, nz-1] + tau;
end;

CubeRotationMatrixGroup := function()

        local mats;

    mats := [];

    # 90 degrees about z: (x,y,z) -> (-y, x, z)
    Add(mats, [ [ 0, -1, 0 ],
                 [ 1,  0, 0 ],
                 [ 0,  0, 1 ] ]);

    # 90 degrees about x: (x,y,z) -> (x, -z, y)
    Add(mats, [ [ 1, 0,  0 ],
                 [ 0, 0, -1 ],
                 [ 0, 1,  0 ] ]);

    # 90 degrees about y: (x,y,z) -> (z, y, -x)
    Add(mats, [ [ 0, 0, 1 ],
                 [ 0, 1, 0 ],
                 [ -1,0,0 ] ]);

    # 120 degrees rotation about the (1,1,1) diagonal: (x,y,z) -> (y,z,x) (order 3)
    Add(mats, [ [ 0, 1, 0 ],
                 [ 0, 0, 1 ],
                 [ 1, 0, 0 ] ]);

     return mats;

end;

PermRotationMat := function(A, m )

        local i, perm, v0, nx, ny, nz, idx, centre, rotated, poslist;

        centre := (m-1)/2;
        poslist := List([0..m-1], i -> i - centre);

        perm := [];
        for i in [1..m^3] do
            v0 := CubeDecode(i, m);   # returns x,y,z in 0..m-1
            # centred coordinates
            v0 := [ v0[1] - centre, v0[2] - centre, v0[3] - centre ];
            rotated := ApplyMat(A, v0); 

            nx := Position(poslist, rotated[1]);
            ny := Position(poslist, rotated[2]);
            nz := Position(poslist, rotated[3]);

            if nx = fail or ny = fail or nz = fail then
                Error("Rotation produced coordinate not on grid; m = ", m,
                      " rotated = ", rotated, " poslist = ", poslist);
            fi;

            # convert back to 0..m-1 indices 
            idx := CubeEncode(nx-1, ny-1, nz-1, m);
            Add(perm, idx);
        od;
        return PermList(perm);

end;

RotationGeneratorsCentered := function(m)
    local gens, mats, A, perm, i, v0, centre, poslist, nx, ny, nz, rotated, idx;

    gens := [];

    # centre as rational: (m-1)/2
    centre := (m-1) / 2;

    poslist := List([0..m-1], i -> i - centre);  

    # Choose rotation matrices 
    mats := [];

    # 90 degrees about z: (x,y,z) -> (-y, x, z)
    Add(mats, [ [ 0, -1, 0 ],
                 [ 1,  0, 0 ],
                 [ 0,  0, 1 ] ]);

    # 90 degrees about x: (x,y,z) -> (x, -z, y)
    Add(mats, [ [ 1, 0,  0 ],
                 [ 0, 0, -1 ],
                 [ 0, 1,  0 ] ]);

    # 90 degrees about y: (x,y,z) -> (z, y, -x)
    Add(mats, [ [ 0, 0, 1 ],
                 [ 0, 1, 0 ],
                 [ -1,0,0 ] ]);

    # 120 degrees rotation about the (1,1,1) diagonal: (x,y,z) -> (y,z,x) (order 3)
    Add(mats, [ [ 0, 1, 0 ],
                 [ 0, 0, 1 ],
                 [ 1, 0, 0 ] ]);

    # For each matrix, produce permutation of 1..m^3
    for A in mats do
        perm := [];
        for i in [1..m^3] do
            v0 := CubeDecode(i, m);   # returns x,y,z in 0..m-1
            # centred coordinates
            v0 := [ v0[1] - centre, v0[2] - centre, v0[3] - centre ];
            rotated := ApplyMat(A, v0);  

            nx := Position(poslist, rotated[1]);
            ny := Position(poslist, rotated[2]);
            nz := Position(poslist, rotated[3]);

            if nx = fail or ny = fail or nz = fail then
                Error("Rotation produced coordinate not on grid; m = ", m,
                      " rotated = ", rotated, " poslist = ", poslist);
            fi;

            # convert back to 0..m-1 indices 
            idx := CubeEncode(nx-1, ny-1, nz-1, m);
            Add(perm, idx);
        od;
        Add(gens, PermList(perm));
    od;

    return Set(gens);   
end;


PermTranslationVector := function (tau, m)
    local   i, perm, vec, newvec;

    perm := [];
    for i in [1..m^3] do
        vec := CubeDecode(i, m);    # 0..m-1
        newvec := [ (vec[1]+tau[1]) mod m,
                    (vec[2]+tau[2]) mod m,
                    (vec[3]+tau[3]) mod m ];
        Add(perm, CubeEncode(newvec[1], newvec[2], newvec[3], m));
    od;

    return PermList(perm);

end;



TranslationGenerators := function(m)
    local gens, move, i, perm, vec, newvec;
    gens := [];
    for move in [ [1,0,0], [0,1,0], [0,0,1] ] do
        perm := [];
        for i in [1..m^3] do
            vec := CubeDecode(i, m);    # 0..m-1
            newvec := [ (vec[1]+move[1]) mod m,
                        (vec[2]+move[2]) mod m,
                        (vec[3]+move[3]) mod m ];
            Add(perm, CubeEncode(newvec[1], newvec[2], newvec[3], m));
        od;
        Add(gens, PermList(perm));
    od;
    return gens;
end;


CubeRotationGroup := function(m)
    local Rgens;
    Rgens := RotationGeneratorsCentered(m);
    return Group(Rgens);
end;

CubeTranslationGroup := function(m)
    local Tgens;
    Tgens := TranslationGenerators(m);
    return Group(Tgens);
end;

CubeFullSymmetryGroup := function(m)
    local Tgens, Rgens;
    Tgens := TranslationGenerators(m);
    Rgens := RotationGeneratorsCentered(m);
    return Group(Concatenation(Tgens, Rgens));
end;


FixedPoints := function ( g, m )
    return Difference( Set( [ 1 .. m ^ 3 ] ), MovedPoints( g ) );
end;


CheckoutGroup := function(n)

    local rgrp, rele, tgrp, cc, tele, centlist, sigma, sigmacoset, fixsigmacoset, i, j, t;


    rgrp :=CubeRotationGroup(n);
    rele := Elements(rgrp);;
    tgrp := CubeTranslationGroup(n);
    tele := Elements(CubeTranslationGroup(n));;
    cc := ConjugacyClasses(rgrp );;
    Print(" Testing Group for n=", n,  ", \n");
    Print(" Rot has conjugacy classes of orders " );
    Print( List(cc, c -> Order( Representative(c))), "\n");
    centlist := List( cc, c -> Centralizer( tgrp, Representative(c) ) );
    Print("Their centraalizers in T have orders " );
    Print( List( centlist, Size ), "\n");
    for i in [ 1 .. Length(cc) ] do
        sigma := Representative( cc[i] );
        if sigma = sigma^0 then continue; fi;
        Print("|sigma| = ", Order(sigma), " number of fixed points = ");
        Print(Length( FixedPoints( Representative( cc[i] ), n)), "\n" );
        sigmacoset := List(tgrp, t-> sigma * t );
        fixsigmacoset := List(sigmacoset, x -> FixedPoints(x,n) );;
        Print("fixed points of coset : ", Collected( List( fixsigmacoset, Length)), "\n");
    od;

end;


FindCycleLengths := function ( n )

	local tele, i, cos, rgrp, cc, tgrp, sigma, res;

        rgrp := CubeRotationGroup(n);
        tgrp := CubeTranslationGroup(n);
	tele := Elements(tgrp);
        cc := ConjugacyClasses( rgrp );
        cc := [List( cc, c -> Representative(c) )[4]];
        res := [];

        Print( List(cc, Order ), "\n");
        for sigma in cc do
		cos := sigma * tele;
                Add( res,  Collected( List( cos, x -> 
                Collected(CycleLengths(x,[1..n^3])))) );
        od;

	return res;

end;



TestTwoa := function( n )


         local cgrp, tgrp, cc, sigma, tele, sigT;

         cgrp := CubeRotationGroup(n);
         tgrp := CubeTranslationGroup(n);
         cc := ConjugacyClasses(cgrp);
         cc := List( cc, c->Representative(c));;
         sigma := cc[3];
         tele := Elements(tgrp);;
         sigT := List(tele, i-> sigma * i );;
         return
         Collected(List( sigT, x-> Collected(CycleLengths(x,[1..n^3]))));

end;


TestTwob := function( n )
         local cgrp, tgrp, cc, sigma, tele, sigT;

         cgrp := CubeRotationGroup(n);
         tgrp := CubeTranslationGroup(n);
         cc := ConjugacyClasses(cgrp);
         cc := List( cc, c->Representative(c));;
         sigma:=cc[5];
         tele := Elements(tgrp);;
         sigT := List(tele, i-> sigma * i );;
         return
         Collected(List( sigT, x-> Collected(CycleLengths(x,[1..n^3]))));

end;



TestProp := function( n )

         local d, m, res, cya, cyb;
         res := [];

         cya := function(m)
                 local s, a;

                 return [[2*m, n^3/(2*m)]];
                 s := "(2*";
                 s := Concatenation(s,String(m));
                 s := Concatenation(s,")^");
                 s := Concatenation(s,String(n^3/(2 * m)));

             return s;
         end;

         cyb := function(m)
                 local s, a;

                 if IsEvenInt(n) then a := 2 * n;
                 else a := n;
                 fi;
                 return [[m,a/m],[2*m,(n^3-a)/(2*m)]];
                 s := String(m);
                 s := Concatenation(s,"^");
                 s := Concatenation(s,String(a/m));
                 s := Concatenation(s,"(2*");
                 s := Concatenation(s,String(m));
                 s := Concatenation(s,")^");
                 s := Concatenation(s,String((n^3-a)/(2 * m)));

             return s;
         end;


         d := DivisorsInt(n);
         for m in d do
                 if IsEvenInt(n) then
                         if IsEvenInt(m) and n/2 mod m = 0 then
                                 Add(res, [cya(m), n^2/2 * Phi(2 * m)]); ### Here changed
                         else
                                 Add(res, [cya(m), n^2/2 * Phi(2*m)]);
                         fi;
                 fi;
                 if not IsEvenInt(m) then
                         if IsEvenInt(n) then
                                 Add(res, [cyb(m), n^2/2 * Phi(m)]);
                         else
                                 Add(res, [cyb(m), n^2 * Phi(m)]);
                         fi;
                 fi;
         od;


         Sort(res);
         return res;

end;



TestPropb := function( n )

         local d, m, res, cya, cyb;
         res := [];

         cya := function(m)
                 local s, a;

                 return [[2*m, n^3/(2*m)]];##n^3/(2*m);

                 s := "(2*";
                 s := Concatenation(s,String(m));
                 s := Concatenation(s,")^");
                 s := Concatenation(s,String(n^3/(2 * m)));

             return s;
         end;

         cyb := function(m)
                 local s, a;

                 if IsEvenInt(n) then
                     return [[m,4*n/m],[2*m,(n^3-4*n)/(2*m)]];
                 else 
	               return [[m,n/m],[2*m,(n^3-n)/(2*m)]];
                 fi;

                 s := String(m);
                 s := Concatenation(s,"^");
                 s := Concatenation(s,String(a/m));
                 s := Concatenation(s,"(2*");
                 s := Concatenation(s,String(m));
                 s := Concatenation(s,")^");
                 s := Concatenation(s,String((n^3-a)/(2 * m)));

             return s;
         end;


         d := DivisorsInt(n);
         for m in d do
                 if IsEvenInt(n) then
	
                         if IsEvenInt(m) and n/2 mod m = 0 then
                                 Add(res, [cya(m), n^2 * Phi(2*m)]); 
                         elif n/2 mod m = 0 then
                                 Add(res, [cya(m), 7*n^2/4*Phi(m)]);
                         fi;
                 fi;
                 if not IsEvenInt(m) then
                         if IsEvenInt(n) then
                                 Add(res, [cyb(m), n^2/4 * Phi(m)]);  ##/4 statt durch /2
                         else
                                 Add(res, [cyb(m), n^2 * Phi(m)]);
                         fi;
                 fi;
         od;


         Sort(res);
         return res;

end;




TestTell:=function(n,l)
    local tgrp,i,tl,t,temp,res,cgrp,cc,sigT,sigma,cyclen;
   


    cgrp := CubeRotationGroup(n);
    cc := ConjugacyClasses(cgrp);
    cc := List( cc, c->Representative(c));;
    sigma := cc[3];

    tgrp := CubeTranslationGroup(n);
    tgrp := Elements(tgrp);;
    res :=[];
    for t in tgrp do
        sigT:=sigma*t;
        cyclen:=Set(CycleLengths(sigT,[1..n^3]));
        temp:=Filtered(cyclen,i -> l mod i =0);
        if Length(temp)<> 0 then
	      Add(res,t);  
	  fi;
    od;

    if IsEvenInt(n) then
	tl:=n^2/2*Gcd(l,n);
    else
	tl:=n^2*Gcd(l,n);
    fi;
##return res;
    return [tl,Length(res)];
end;



Test2:=function(n,l)
    local tgrp,i,tl,t,temp,res,cgrp,cc,sigT,sigma,cyclen;
   


    cgrp := CubeRotationGroup(n);
    cc := ConjugacyClasses(cgrp);
    cc := List( cc, c->Representative(c));;
    sigma := cc[3];

    tgrp := CubeTranslationGroup(n);
    tgrp := Elements(tgrp);;
    res :=[];
    for t in tgrp do
        sigT:=sigma*t;
	  Add(res,Order(sigT));  
    od;
    return Set(res);
end;







