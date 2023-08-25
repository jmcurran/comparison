## `comparison`: A package for computing likelihood ratios for univariate and multivariate evidence.

This package is for computing the *weight of the evidence*, i.e. the
*likelihood ratio (LR)* for trace evidence which has been quantified with
some instrument. For example a forensic scientist might be have determined
the refractive indices of fragments of glass taken from a crime scene and
fragments of glass recovered from the clothing of the suspected breaker. This package evaluates the probability (density) of the evidence, $E$, (the RI values from the two samples) under the hypothesis $H_p$ that they
originated from the same source, and alternatively under the hypothesis
$H_d$ that they originated from another source. The $LR$ is the ratio of
these two quantities, i.e. 
$$LR = \frac{p(E|H_p)}{p(E|H_d)}.$$
A $LR$ which is greater than one indicates that the evidence supports $H_p$, and a $LR$ which is less than one indicates that the evidence supports $H_d$.

The computation can use either univariate or multivariate observations of a
physical object. For example trace element measurements, and a similar set of uni/multivariate observations from another object, and calculates a
likelihood ratio for the propositions that the first item came from the same source as the second given some population data.
 
## Acknowledgements
In a package of functions such as these which have undergone a long
development over a number of years, it is inevitable that a number of people, besides those directly cited, have helped to correct and add to the code. These people are (in alphabetical order): Ivo Alberink (NFI), Anabel  Bolck(NFI), Sonja Menges (BKA), Geoff Morrison (Aston), Tereza Neocleous
(Glasgow), Anders Nordgaard (SKL), Brad Patterson (George Mason), Phil Rose
(ANU), Agnieszka Rzepecka (Jagiellonian), Marjan Sjerps (NFI) and Hanjing
Zhang (Edinburgh).

## References 
Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in
  the form of multivariate data. Applied Statistics: 53(1):109-122.
  
