<TeXmacs|1.0.7.18>

<style|generic>

<\body>
  <math|S<around*|(|k|)>> is defined as:

  <\equation*>
    S<around*|(|k|)>\<equiv\><big|sum><rsub|r>e<rsup|-\<beta\> k
    \<varepsilon\><rsub|r>>,
  </equation*>

  where the sum is over all single-particle states <math|r>.

  For a <math|d>-dimensional isotropic harmonic oscillator with trap
  frequency <math|\<omega\>>:

  <\equation*>
    \<varepsilon\>=<around*|(|n<rsub|x<rsub|1>>+\<cdots\>+n<rsub|x<rsub|d>>+d/2|)>\<hbar\>\<omega\>,<space|2em>n<rsub|x<rsub|k>>=0,1,2,\<ldots\>
    .
  </equation*>

  We will measure the temperature <math|T=1/\<beta\>> in units of
  <math|\<hbar\>\<omega\>>, so <math|\<beta\>*\<varepsilon\><rsub|r>> is
  dimensionless.

  For this system:

  <\eqnarray*>
    <tformat|<table|<row|<cell|S<around*|(|k|)>>|<cell|=>|<cell|<around*|[|<big|sum><rsub|n=0><rsup|\<infty\>>e<rsup|-\<beta\>
    k <around*|(|n+1/2|)>>|]><rsup|d>>>|<row|<cell|>|<cell|=>|<cell|<around*|[|<frac|e<rsup|-\<beta\>*k/2>|1-e<rsup|-\<beta\>*k>>|]><rsup|d>>>|<row|<cell|>|<cell|=>|<cell|<block|<tformat|<table|<row|<cell|<around*|[|2*sinh<around*|(|\<beta\>*k/2|)>|]><rsup|-d>>>>>>.>>>>
  </eqnarray*>

  Therefore:

  <\eqnarray*>
    <tformat|<table|<row|<cell|S<rprime|'><around*|(|k|)>\<equiv\><frac|\<partial\>S<around*|(|k|)>|\<partial\>\<beta\>>>|<cell|=>|<cell|<frac|-k*d*|2<rsup|d+1>>*<around*|[|sinh<around*|(|\<beta\>*k/2|)>|]><rsup|-<around*|(|d+1|)>>*cosh<around*|(|\<beta\>k/2|)>>>|<row|<cell|>|<cell|=>|<cell|<block|<tformat|<table|<row|<cell|<frac|-k*d|2>*S<around*|(|k|)>*coth<around*|(|\<beta\>k/2|)>>>>>>,>>|<row|<cell|S<rprime|''><around*|(|k|)>\<equiv\><frac|\<partial\><rsup|2>S<around*|(|k|)>|\<partial\>\<beta\><rsup|2>>>|<cell|=>|<cell|<frac|-k*d|2>*<around*|[|S<rprime|'><around*|(|k|)>*coth<around*|(|\<beta\>k/2|)>-S<around*|(|k|)>*<frac|k|2>
    csch<rsup|2><around*|(|\<beta\>k/2|)>|]>>>|<row|<cell|>|<cell|=>|<cell|d<around*|(|<frac|k|2>|)><rsup|2>S<around*|(|k|)><around*|[|d**coth<rsup|2><around*|(|\<beta\>k/2|)>+csch<rsup|2><around*|(|\<beta\>k/2|)>|]>>>|<row|<cell|>|<cell|=>|<cell|<block|<tformat|<table|<row|<cell|d<around*|(|<frac|k|2>|)><rsup|2>S<around*|(|k|)><frac|d*cosh<rsup|2><around*|(|\<beta\>k/2|)>+1|sinh<rsup|2><around*|(|\<beta\>k/2|)>>>>>>>.>>>>
  </eqnarray*>

  All of these formulas are numerically stable, assuming accurate
  calculations of the hyperbolic functions.

  However, the calculation of the <math|N>-particle partition function
  <math|Z<rsub|N>> is unstable at low temperature and requires high-precision
  arithmetic. (An alternative method at low temperature would be to use
  particle-number projection.)

  \;

  \;

  \;
</body>