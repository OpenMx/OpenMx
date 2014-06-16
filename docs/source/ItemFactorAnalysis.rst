.. _Item Factor Analysis:

Item Factor Analysis
********************

What is Item Factor Analysis?
=============================

..
	library(rpf)
	library(reshape2)
	library(gridExtra)
	library(ggplot2)
	skill.n <- 500
	width <- 5
	skill <- sort(runif(skill.n,-width,width))
	item.p <- .4
	empirical.bins <- 20
	correct <- rep(TRUE, length(skill))
	skill.o <- skill + rnorm(length(skill), sd=1)
	correct[order(skill.o)[seq(1,(1-item.p) * length(skill))]] <- FALSE
	grid <- list()
	grid$correct <- correct
	grid$skill <- skill
	breaks <- seq(min(skill)-.001, max(skill)+.001, length.out=empirical.bins)
	bin <- cut(skill, breaks, labels=FALSE)
	bin.correct <- data.frame(at=breaks[-1] - diff(breaks)/2,
	 pct=vapply(1:max(bin), function(l) sum(correct[bin==l])/sum(bin==l), 0))
	bin.correct$pct <- sprintf("%.0f", 100 * bin.correct$pct)
	pl <- ggplot(as.data.frame(grid), aes(skill, correct)) +
	  geom_point(position=position_jitter(0,.05), size=1) +
	  geom_segment(data=data.frame(thresh=breaks),
	               aes(x=thresh, xend=thresh, y=TRUE, yend=FALSE), color="red") +
	  geom_text(data=bin.correct, aes(x=at,y=1.5,label=pct), color="blue",
	            angle=90) +
	  labs(y="% correct")
	pdf("cache/ifa-intervalize.pdf", height=2.5)
	print(pl)
	dev.off()
	png("cache/ifa-intervalize.png", height=180)
	print(pl)
	dev.off()

.. _figure-ifa-intervalize:
.. figure:: cache/ifa-intervalize.*

	    Percentage correct of responses by skill bin

	    Bins are demarcated by red vertical bars.
	    This procedure converts nominal or ordinal data (the black dots) into
	    interval scale data conditional on examinee skill.
	    The blue numbers in the middle of the plot are the interval data.

Educational assessments are often scored *correct* and *incorrect*.
Some items may also offer partial credit for partially correct answers.
One approach to analyze these ordinal data is Classical Test Theory (CTT).
CTT includes statistics like Cronbach's α for estimating reliability,
proportion correct as an estimate of item difficulty,
and item-total correlation as an estimate of item discrimination.
However, a disadvantage of CTT is the need for huge
samples to establish test norms.
Modern Test Theory, also known as Item Response Theory or Item Factor Analysis (IFA),
offers substantial advantages over CTT
such as the ability to obtain unbiased item parameter estimates
from an unrepresentative sample [Embretson1996]_.

..
	plot.icc <- function(item, param, width=3) {
	pm <- t(rpf.prob(item, param, seq(-width, width, .1)))
	icc <- as.data.frame(melt(pm, varnames=c("theta",'category')))
	icc$theta <- seq(-width, width, .1)
	icc$category <- as.factor(icc$category - 1)
	ggplot(icc, aes(theta, value)) +
	geom_line(aes(color=category, linetype=category)) +
	ylim(0,1) + xlim(-width,width) + labs(y="Probability", x="Theta")
	}
	i1 <- c(1, -1.5, logit(0), logit(1))
	i2 <- c(1, .9, logit(0), logit(1))
	pdf("cache/ifa-1pl-rpf.pdf", height=2.5)
	grid.arrange(plot.icc(rpf.drm(), i1, width=6) + labs(title="a."),
	      plot.icc(rpf.drm(), i2, width=6) + labs(title="b."), ncol=2)
	dev.off()
	png("cache/ifa-1pl-rpf.png", height=180)
	grid.arrange(plot.icc(rpf.drm(), i1, width=6) + labs(title="a."),
	      plot.icc(rpf.drm(), i2, width=6) + labs(title="b."), ncol=2)
	dev.off()

.. _figure-ifa-1pl-rpf:
.. figure:: cache/ifa-1pl-rpf.*

	Example plots of the logistic curve

	Example plots of the logistic curve (blue) and its complement (red).
	The interval data are fit against the blue curve.
	The *c* parameter is set to 1.5 for (a) and -0.9 for (b).
	Which curve would fit the data better? Answer: Curve (a) is a better fit
	because the 50% level occurs near skill=1.5 in :ref:`figure-ifa-intervalize`.

IFA is based on the idea of fitting models to data. The first step
is to convert nominal or ordinal data into interval scale data. To understand how this is
accomplished, assume that the true skill (traditionally :math:`\theta`) of participants is known. We can
plot item outcome by true skill, partition the skill distribution into bins, and count the
proportion of correct responses in each bin (:ref:`figure-ifa-intervalize`). We can fit a model to these data.
One popular model is the logistic model (also called the 1PL or Rasch model),

.. math::
  Pr(pick=0|c,\theta) = 1 - Pr(pick=1|c,\theta) \\
  Pr(pick=1|c,\theta) = \frac{1}{1+\exp(-(\theta-c))}

where :math:`\theta` is the participant's skill and *c* is the estimated parameter to describe the item's
difficulty. See :ref:`figure-ifa-1pl-rpf`. While the examinee's
true skill is assumed known in this simplified introduction, such an assumption is unnecessary.
[BockAitkin1981]_ can estimate both item parameters and person skills
simultaneously.

Simple Models
=============

This chapter assumes that you have already read the :ref:`BasicIntroduction`.

A Rasch model
-------------

Suppose you regularly administer the PANAS [WatsonEtal1988]_, but
instead of scoring participants by adding up the item scores, you want
to try IFA. Here is how you might do it. To simplify this example, we
will only consider the positive affect part of the scale.

.. code-block:: r

   library(OpenMx)
   library(rpf)

   PANASItem <- c("Very Slightly or Not at All",  "A Little",
		"Moderately", "Quite a Bit",	"Extremely")
   spec <- list()
   spec[1:10] <- rpf.grm(outcomes = length(PANASItem))

   # replace with your own data
   data <- rpf.sample(500, spec, sapply(spec, rpf.rparam))

   for (cx in 1:10) levels(data[[cx]]) <- PANASItem          # repair level labels
   colnames(data) <- c("interested", "excited", "strong", "enthusiastic", "proud",
                       "alert", "inspired", "determined", "attentive", "active")
   head(data)  # much easier to understand with labels
   origData <- data
   
   startingValues <- matrix(c(1, seq(1,-1,length.out=4)), ncol=length(spec), nrow=5)
   ip.mat <- mxMatrix(name="ItemParam", values=startingValues,
                   free=TRUE, dimnames=list(names(rpf.rparam(spec[[1]])), colnames(data)))
   rownames(ip.mat)[1] <- "posAff"
   ip.mat$labels[1,] <- 'slope'

   # assume a standard Normal latent distribution
   m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0,
		free=FALSE, dimnames=list('posAff', 'posAff'))
   cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1,
		free=FALSE, dimnames=list('posAff', 'posAff'))

   originalPanas1 <- mxModel(model="panas1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=spec, ItemParam="ItemParam", mean="mean", cov="cov"),
              mxFitFunctionML(),
	      mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()))
   panas1 <- mxRun(originalPanas1)

A PANAS item always has 5 possible outcomes, but best practice is to
list the outcomes labels and let R count them for you. If you work on a
new measure then you will appreciate the ease with which your script
can adapt to adding or removing outcomes. Even for PANAS, we could try
collapsing two outcomes and see how the model fit changes.

.. code-block:: r

   spec <- list()
   spec[1:10] <- rpf.grm(outcomes = length(PANASItem))

The ``rpf.grm`` function creates an ``rpf.base`` class object that
represents an item response function. An item response function
assigns probabilities to response outcomes. The ``grm`` in the
function name ``rpf.grm`` stands for *graded response model*. You can inspect the
mathematical formula for the graded response model by requesting the
manual page with ``?rpf.grm``.

.. code-block:: r

   data <- rpf.sample(500, spec, sapply(spec, rpf.rparam))

This line creates some fake data for 500 random participants based on
our list of item models and random item parameters.  Instead of this
line, you would typically read in your data using ``read.csv`` and
convert it to ordered factors using ``mxFactor``.

.. code-block:: r

   startingValues <- matrix(c(1, seq(1,-1,length.out=4)), ncol=length(spec), nrow=5)

We can input particular starting values. That is what we do here. The
graded response model is a little finicky; the threshold parameters
must be ordered. Alternately, a good way to obtain random starting
values is with,

.. code-block:: r

   startingValues <- mxSimplify2Array(lapply(spec, rpf.rparam))
   startingValues[1,] <- 1  # these parameters must be equal

This is convenient because it will work for a non-homogeneous list of
items. If you need to set some starting values to something specific then
you might start with random starting values and then override any rows
and columns as needed.

.. code-block:: r

   ip.mat <- mxMatrix(name="ItemParam", values=startingValues,
                   free=TRUE, dimnames=list(names(rpf.rparam(spec[[1]])), colnames(data)))
   rownames(ip.mat)[1] <- "posAff"

The ItemParam matrix contains response probability function parameters
in columns. This layout can be a little awkward when you estimate a
mixed format measure with different numbers of outcomes for each item.
Fortunately, all the PANAS items are the same. We must label all the
rows and columns. The first row must be labeled with the name of our
factor, ``posAff``. The remaining rows can take any label. Since all
of our items are the same, we can use the default item parameter
names. The column names must match the column names of the data.

.. code-block:: r

   ip.mat$labels[1,] <- 'slope'

Here we set the label of every parameter in the first row to
``slope``. This acts like an equality constraint. The effect is that
we assume all items work equally well at measuring the latent trait.
This constraint is what makes the difference between a Rasch model and
any other kind of IFA model. A Rasch model makes this assumption.

.. code-block:: r

   m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0,
		free=FALSE, dimnames=list('posAff', 'posAff'))
   cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1,
		free=FALSE, dimnames=list('posAff', 'posAff'))

We must specify the distribution of our latent factor. Following
tradition, we assume that ``posAff`` has a standard Normal
distribution. We will come back to discuss the latent distribution
when we consider multigroup models. In multigroup models, some of
these parameters can be estimated. For this model, our latent
parameters are fixed.

.. code-block:: r

   panas1 <- mxModel(model="panas1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=spec, ItemParam="ItemParam", mean="mean", cov="cov"),
              mxFitFunctionML(),
	      mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()))

Here we put everything together. There are a few things that are new.

.. code-block:: r

   mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson())

The ``mxComputeEM`` plan is a somewhat more sophisticated version of

.. code-block:: r

   mxComputeIterate(steps=list(
                    mxComputeOnce('expectation', 'scores'),
                    mxComputeNewtonRaphson(),
                    mxComputeOnce('expectation'),
                    mxComputeOnce('fitfunction', 'fit')))

In both compute plans, OpenMx will iterate until the change in the
fitfunction is less than some threshold. For every iteration, person
scores are predicted by the expectation, a Newton-Raphson optimization
takes place to improve the item parameter values, the predicted person
scores are discarded, and the fitfunction is evaluated. In comparison
to ``mxComputeIterate``, the ``mxComputeEM`` plan offers additional
options to speed up convergence and estimate standard errors.

After running the model, we can inspect the parameters estimates,

.. code-block:: r

   panas1 <- mxRun(originalPanas1)
   panas1$ItemParam$values
   # or
   summary(panas1)

At this point, you might notice something unsettling about the summary
output. No standard errors are reported. How do we know whether our
model converged? Excellent question. There are a few things that we
can check. We can look at the count of EM cycles and M-step
Newton-Raphson iterations.

.. code-block:: r

   panas1$compute$output

If the number of EM cycles is 2 or less then it is likely that the
parameters are still sitting at their starting values. Since our
starting values were mostly random, that's probably not the solution
we were looking for. Instead of digging into the ``MxComputeEM`` output, we
can also re-run the model with extra diagnostics enabled.

.. code-block:: r

   panas1 <- mxModel(model=originalPanas1,
                  mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson(), verbose=2L))
   panas1 <- mxRun(panas1)

Note the addition of ``verbose=2L`` to ``mxComputeEM``.
The ``mxComputeNewtonRaphson`` object takes a ``verbose`` parameter as
well if you want to examine the progress of the optimization in (much)
more detail. From this diagnostic output, we can discern that the
parameters are changing, but we still do not know whether the solution
is a candidate global optimum.

.. code-block:: r

  info1 <- mxModel(panas1,
                mxComputeSequence(steps=list(
                  mxComputeOnce('fitfunction', 'information', "meat"),
                  mxComputeStandardError(),
                  mxComputeHessianQuality())))
  info1 <- mxRun(info1)

As a starting point, we can estimate the information matrix by
taking the inverse of the covariance of the first derivatives.
This estimate is called ``meat`` because it forms the inside
of a sandwich covariance matrix [White1994]_.
The first thing to look at is the condition number of the
information matrix.

.. code-block:: r

  info1$output$conditionNumber
  # 525.8942     # yours may be different

A finite condition number implies that the information matrix is
positive definite.  Since the condition number is roughly closer to
zero than to positive infinity, there is a good chance that the parameters
are at a candidate global optimum. We can examine the standard errors.

.. code-block:: r

  summary(info1)

Further diagnostics are available from the `RPF package <http://cran.r-project.org/web/packages/rpf/index.html>`_.
Many of these
diagnostic functions are most convenient to use when all the relevant
information is packaged up into an IFA group object. An IFA group is
not an object in the usual R sense, but you can think of it like an
object. The problem with R objects is that they are a little
mysterious. IFA groups are deliberately designed as simple lists to
eliminate the mystery and encourage interoperability between IFA
software.

.. code-block:: r

   grp1 <- list(spec=panas1$expectation$ItemSpec,
            param=panas1$ItemParam$values,
            mean=panas1$mean$values,
            cov=panas1$cov$values,
            data=origData)

Note that we used ``origData`` instead of ``panas1$data$observed``.
That is because the observed data in the model has been sorted by ``mxRun``.

A fundamental assumption of IFA is that items are conditionally
independent. That is, the outcome on a given item only depends on its
item parameters and examinee skill, not on the outcome of other items.
At least some attempt should be made to check this assumption.

.. code-block:: r

   ChenThissen1997(grp1)

Item pairs that exhibit statistically significant local dependence
and positively correlated residuals should be investigated.
If ignored, local dependence exaggerates the accuracy
of measurement.
Standard errors will be smaller and
items will seem to fit the data better than they otherwise would [Yen1993]_.

.. code-block:: r

   score1 <- mxModel(panas1,
                  mxExpectationBA81(ItemSpec=spec, ItemParam="ItemParam",
		    mean="mean", cov="cov", scores="full"),
                  mxComputeOnce('expectation'))

   score1 <- mxRun(score1)
   head(score1$expectation$output$scores)

Since we have a single factor Rasch model, the residuals are easily
interpretable. We can examine Rasch fit statistics *infit* and
*outfit*.  Before we do that, however, we need to compute EAP
scores. To get EAP scores, we need to add ``scores="full"`` to
``mxExpectationBA81``. We could have done this from the beginning, but
sometimes the extra overhead of computing EAP scores is undesirable.
Note that the scores are in the original data order, not the
sorted data order.
 
.. code-block:: r

   grp1$scores <- score1$expectation$output$scores

   rpf.1dim.fit(group=grp1, margin=2)

..
	item.map <- function(grp, factor=1) {
	item.mask <- grp$param[factor,] > 0
	result <- NULL
	for (ix in rev(colnames(grp$param)[item.mask])) {
	  lev <- levels(grp$data[,ix])
	  for (ox in 1:length(lev)) {
	  mask <- grp$data[,ix]==lev[ox]
	  mask <- !is.na(mask) & mask
	  if (all(!mask)) next
	  result <- rbind(result, data.frame(item=ix,
                                         outcome=ox, outcome.name=lev[ox],
                                         score=mean(grp$scores[mask, factor], na.rm=TRUE)))
	  }
	}
	result
	}
	map1 <- item.map(grp1, 1)
	pl <- ggplot(map1, aes(x=score, y=item, label=outcome)) + geom_text(size=4, position=position_jitter(h=.25))
	pdf("cache/ifa-1pl-itemMap.pdf", height=2.5)
	print(pl)
	dev.off()
	png("cache/ifa-1pl-itemMap.png", height=180)
	print(pl)
	dev.off()

.. _figure-ifa-1pl-itemMap:
.. figure:: cache/ifa-1pl-itemMap.*

	Example item map

	Outcomes located at the mean of the ability of every examinee
	who picked that outcome.

This gives us item-wise statistics. For person-wise statistics, we can
replace ``margin=2`` with ``margin=1``. For some discussion on the
interpretation of these statistics, visit the `Infit and Outfit page
<http://www.rasch.org/rmt/rmt162f.htm>`_ at the Institute for
Objective Measurement.
Another way to look at the results is to create an item plot.
An item plot assigns to every outcome the mean of the ability of
every examinee who picked that outcome (:ref:`figure-ifa-1pl-itemMap`).

.. code-block:: r

   sumScoreEAP(grp1)

Finally, we can generate a sum-score EAP table. The ``posAff`` column
contains the interval-scale score corresponding to the row-wise
sum-score. You could use this table to score the PANAS instead of
merely using the sum-score. The EAP sum-score will likely provide
higher accuracy measurement of the latent trait *positive affect*.

A 2PL model
-------------

Suppose you are skeptical that all positive affect items work equally
well at measuring positive affect.
We can relax this assumption and let the optimizer estimate
how well each items is working.
Continuing our previous example,
all that is needed is to remove the label from the item parameter matrix.

.. code-block:: r

   panas2 <- mxModel(model=panas1,
                mxComputeSequence(list(
		mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()),
                mxComputeOnce('fitfunction', 'information', "meat"),
                mxComputeStandardError(),
                mxComputeHessianQuality())))
   panas2$ItemParam$labels[1,] <- NA
   panas2 <- mxRun(panas2)

The compute plan is the same as before except that we combined both
the ``mxComputeEM`` fit step and computation of standard errors
in a single ``mxComputeSequence``.
As usual, we start by inspecting the condition number of the Hessian
then look at the parameter estimates with standard errors.

.. code-block:: r

   panas2$output$conditionNumber
   summary(panas2)

Since these models are nested, we can conduct a likelihood ratio test
to determine whether ``panas2`` fits the data significantly better than
our Rasch constrained model ``panas1``.

.. code-block:: r

   mxCompare(panas2, panas1)

Rasch fit statistics are not appropriate for a 2PL model because
residuals from different items have different weights. However, we can
compare expected item outcome proportions against sum-scores.  First
we need to create IFA group object for ``panas2`` then we can run the
S test [OrlandoThissen2000]_.

.. code-block:: r

   grp2 <- list(spec=panas2$expectation$ItemSpec,
            param=panas2$ItemParam$values,
            mean=panas2$mean$values,
            cov=panas2$cov$values,
            data=origData,
	    free=panas2$ItemParam$free)
   SitemFit(grp2)
   
The additional ``free`` logical matrix is included in the group for
a degrees of freedom adjustment to the test.
However, this adjustment is somewhat controversial and
typically makes little difference in the outcome of the test.


Additional possible topics:

* 3PL with Bayesian priors (with likelihood-based CIs)
* multiple group models
* special features to cope with missing data
* automatic-ish item construction (especially for the nominal model)
* two-tier models
* simulation studies
* more than 1 primary factor
* more plots

.. [BockAitkin1981] Bock, R. D. & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters:
		    Application of an EM algorithm. Psychometrika, 46, 443–459.

.. [Embretson1996] Embretson, S. E. (1996). The new rules of measurement. Psychological Assessment, 8(4), 341-349.

.. [OrlandoThissen2000] Orlando, M. and Thissen, D. (2000). Likelihood-Based
			Item-Fit Indices for Dichotomous Item Response Theory Models.
			Applied Psychological Measurement, 24(1), 50-64.

.. [WatsonEtal1988] Watson, D., Clark, L. A., & Tellegen, A. (1988). Development and validation of brief
		    measures of positive and negative affect: The PANAS scales. Journal of Personality
		    and Social Psychology, 54 (6), 1063.

.. [White1994] Estimation, Inference and Specification
	       Analysis. Cambridge University Press, Cambridge.

.. [Yen1993] Yen, W. M. (1993). Scaling performance assessments:
             Strategies for managing local item dependence. Journal of
             Educational Measurement, 30, 187-213.
