# commander2.5
Generelization of Commander2 code for frequency-correlated noise

 1. Each processor handles a single main band, although importing every noise covariance 
    and computing residual/signal maps for all bands when a N^-1 d operation has to be 
    performed. Other quantities are expanded to every band, such as beams, foreground
    spectral responses and templates.
2. Ft N^-1 F computation is simplified, limiting the Ft N^-1 F solving to a vector-scalar
   product rather then a matrix-vector one.
3. Chi^2 computation formula has been changed, replacing the sqrt{N^-1} multiplica-
   tion with the standard N^-1 one. The reason is that the new routine treats 
   differently the rhs e and the lhs map terms, so that the multiplication of the two
   no longer equals lhs^2.


For more information check out the QUALITATIVE_LOG.md and QUANTITATIVE_LOG.md files
(given in order of how detailed the list of modifications is)
