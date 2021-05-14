default(parisizemax, 100m);
g = Mod(6, 682492462409094395392022581537473179285250139967739310024802121913471471);
A = 245036439927702828116237663546936021015004354074422410966568949608523157;

n = 682492462409094395392022581537473179285250139967739310024802121913471471;



\\ On a A = g^a et on cherche 'a' par le log discret.
\\ On suppose a = x0 + x1 B.

\\ Calcul des g^x0.

init_Baby(g, racine) = {
  baby = vector(racine);
  baby[1] = lift(g);

  for(i = 2, racine,
    baby[i] = lift( baby[i-1] * g )
  );

  return(baby);
};



\\ Calcul des 'A x G^x1'.
\\ On s'arrête dès qu'on a une collision avec la liste des g^x0.

BabyStep_GiantStep(A, g, baby) = {

  my(B, ordre, G, x1);

  B = #baby;

  \\ baby steps

  ordre = vecsort(baby, , 1);             \\ Ordonne la liste baby.
  baby = vector(B, i, baby[ordre[i]] );

  G = g^(-B);

  \\ giant steps

  for(j = 0, B - 1,
    i = vecsearch( baby, lift(A) );
    if(i,
      return( ordre[i] + B * j) );      \\ On a donc a = x0 + x1 * B.

    A = A * G;
  );

  return(0);
};



\\ POHLIG-HELLMAN
\\ puissance du facteur premier p est 1

Puiss_unique(p, A, g, n) = {
  my(d, baby);

  d = n/p;
  A = A^d;

  if(A == 1,
    return(0)
  );

  g = g^d;
  B = sqrtint(p) + 1;
  baby = init_Baby(g, B);

  a = Mod( BabyStep_GiantStep(A, g, baby), p );

  return(lift(Mod(d, p)^-1) * d * lift(a));
};


\\ POHLIG-HELLMAN
\\ puissance du facteur premier p est e != 1

Puiss_multiple(p, e, A, g, n) = {

  my(d, B, g_tmp, baby, ai, Ai, res);

  d = n/(p^e);
  A = A^d;
  g = g^d;
  B = sqrtint(p) +1;

  g_tmp = g^( p^(e-1) ) ;

  baby = init_Baby(g_tmp, B) ;

  ai = vector(e);
  ai[1] = BabyStep_GiantStep( A^( p^(e-1) ), g_tmp, baby ) ;

  Ai = g^ai[1] ;

  for(i = 0, e-2,
    ai[i+2] = BabyStep_GiantStep( (A * (Ai^(-1)))^( p^(e-i-2) ), g_tmp, baby );
    Ai = Ai * g^( ai[i+2] * p^(i+1) );
  );

  res = Mod(0, p^e);
  for(i = 1, e,
    res += ai[i] * p^(i-1);
  );

  return(lift(Mod(d, p^e)^-1) * d * lift(res));
};






Log_Discret() = {

  my(n, g, A, facteurs, res);

  n = 682492462409094395392022581537473179285250139967739310024802121913471471;

  g = Mod(6, n);
  A = Mod(245036439927702828116237663546936021015004354074422410966568949608523157, n);

  facteurs = factor(n-1);
  res = Mod(0, n-1);

  for(i = 1, #facteurs[, 1],
    if(facteurs[i, 2] == 1,
      res += Puiss_unique(facteurs[i, 1], A, g, n-1),

    \\ else
      res += Puiss_multiple(facteurs[i, 1], facteurs[i, 2], A, g, n-1);
    );
  );

  return(lift(res));

};

print(Log_Discret());
