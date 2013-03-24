<TeXmacs|1.0.7.18>

<style|generic>

<\body>
  <section|Calculation of <math|S<around*|(|k|)>>>

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

  <section|Calculation of observables>

  <\equation*>
    Z<rsub|><around*|(|N|)>=<frac|1|N><big|sum><rsub|k=1><rsup|N><around*|(|\<um\>|)><rsup|k+1>S<around*|(|k|)>Z<around*|(|N-k|)><space|0.2spc>,<space|2em>Z<around*|(|0|)>=1.
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|E<rsub|><around*|(|N|)>>|<cell|=>|<cell|<frac|-1|Z<rsub|><around*|(|N|)>>*<frac|\<partial\>Z|\<partial\>\<beta\>>=:<frac|-1|Z<around*|(|N|)>>
    Z<rprime|'><rsub|><around*|(|N|)>.>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|Z<rprime|'><rsub|><around*|(|N|)>>|<cell|=>|<cell|<frac|1|N><big|sum><rsub|k=1><rsup|N><around*|(|\<um\>|)><rsup|k+1><around*|[|S<rprime|'><around*|(|k|)>Z<around*|(|N-k|)>+S<around*|(|k|)>Z<rprime|'><around*|(|N-k|)>|]>.>>>>
  </eqnarray*>

  <section|Particle-number projection>

  <\equation*>
    Z<rsub|N>=<frac|1|2\<pi\>><big|int><rsub|0><rsup|2\<pi\>>d\<varphi\>*e<rsup|i\<varphi\>N>e<rsup|\<alpha\>N>*\<zeta\><around*|(|\<varphi\>,\<alpha\>|)>,
  </equation*>

  where

  <\equation*>
    \<zeta\><around*|(|\<varphi\>,\<alpha\>|)>=Tr<rsub|GC><around*|(|e<rsup|-\<beta\><wide|h|^>>e<rsup|-i\<varphi\><wide|N|^>>e<rsup|-\<alpha\><wide|N|^>>|)>=<big|prod><rsub|k=1><rsup|\<infty\>><around*|(|1+e<rsup|-\<beta\>*\<varepsilon\><rsub|k>>e<rsup|-i\<varphi\>>e<rsup|-\<alpha\>>|)><space|0.2spc>.
  </equation*>

  This formula is more appropriate at low temperatures, where only a few
  terms in the product are important (<math|\<beta\>> is large). At high
  temperatures, the product becomes more difficult to compute.

  \;
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Calculation
      of <with|mode|<quote|math>|S<around*|(|k|)>>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Calculation
      of observables> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>