This is the repository of ComplexSim package, which is an R package for simulating longitudial omics data with user defined structures. 

The simulation function includes the Gaussian and Negative Binomial distributions and takes into consideration the treatment effect, time effect, subject level random effect as well as the sample correlations. 

RNA sequencing data in discrete count unit can be described by a negative binomial distribution.  
 The following Delta method is introduced to illustrate the link between the mean and variance of the two distributions, which guarantees the equivalence of the two simulation functions.

Following the RNA sequencing procedure, we build the negative binomial model for the mapped read counts: 

$$ \mathrm{NB}(r,p):  \mathrm{Pr}(X=k) = \binom{k+r-1}{k}\cdot(1-p)^{r} p^{k} . $$ 

The mean and variance of this negative binomial distribution are 

$$ \mu= \frac{p r}{1-p} \quad \text{ and } \quad \sigma^2 = \frac{pr}{(1-p)^2}. $$

One important proposition of the gene expression count data is the over-despersion property. We can use the variance-mean ratio or taking the difference between the variance and the mean to quantify the dispersion. The ratio is given by

$$ \frac{\sigma^2}{\mu} = \frac{1}{1-p} >1, $$

and the difference is 
$$
 \sigma^2 - \mu = p\cdot \frac{pr}{(1-p)^2} = p \sigma^2 \geq 0\quad \text{or equivalently} \quad 
 \sigma^2 - \mu = \frac{1}{r} \cdot \frac{p^2 r^2}{(1-p)^2} = \frac{1}{r} \mu^2 \geq 0 . 
$$
From the ratio and the difference, it is obvious to see that the counts of mapped reads which follow the negative binomial distribution is over dispersed. 
The variance can be decomposed to

$$ \sigma^2 = \mu + \frac{1}{r} \mu^2 . $$

If $`r \to \infty`$, $`\sigma^2 = \mu`$ corresponds to Poisson case. 

Alternatively, one can normalize the count data by using the log-cpm described as follows. Let $`(r_{gi})_{G \times n}`$ be the matrix of counts for $`n`$ RNA samples and $`G`$ genes.
 For sample $`i`$, the total number of mapped reads is $`R_{i}=\sum_{g=1}^{G} r_{g i}`$, and we usually compute the log-counts per million (log-cpm): 
 
$$ y_{g i}=\log _{2}\left(\frac{r_{g i}+0.5}{R_{i}+1.0} \times 10^{-6}\right). $$ 

This re-scaling of the count data can be considered as the simplest way of normalization, and the normalized data is usually considered as Gaussian. We use both Gaussian and Negative Binomial distribution to simulate the normalized and raw gene expression data.

Consider a read count $`r`$ and library size $`R`$.  
The expected value is defined as  $`\lambda=\textbf{E}(r)`$, and we obtain the variance decomposition 

$$ \textbf{Var}(r)=\lambda+\phi \lambda^{2}, $$ 

which corresponds to $` \sigma^2 = \mu + \frac{1}{r} \mu^2 `$ for negative binomial distribution.
Therefore, the log-cpm value  can be approximated by

 $$ y \approx \log _{2}(r) -\log _{2}(R)+6 \log _{2}(10), $$
 
and has variance

 $$\textbf{Var}(y) \approx \operatorname{Var}\left(\log _{2}(r)\right) \approx \frac{\operatorname{Var}(r)}{\lambda^{2}}=\frac{1}{\lambda}+\phi$$  
 
 following the Delta-method.
The mean and variance of the normalized sequencing data is linked to the original mean and dispersion of the negative binomial data and can be modeled and estimated. 

In this R package, we allow simulation using both Gaussian and Negative Binomial distribution with different set ups of the parameters, that is, mean and variance of Gaussian, mean and dispersion parameter of Negative Binomial. One can switch between the two distributions to get equivalent simulations. 

