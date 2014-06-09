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

   panas1 <- mxModel(model="panas1", ip.mat, m.mat, cov.mat,
              mxData(observed=data, type="raw"),
              mxExpectationBA81(ItemSpec=spec, ItemParam="ItemParam", mean="mean", cov="cov"),
              mxFitFunctionML(),
	      mxComputeEM('expectation', 'scores', mxComputeNewtonRaphson()))
   panas1 <- mxRun(panas1)

A PANAS item always has 5 possible outcomes, but best practice is to
list the outcomes labels and let R count them for you. If you work on a
new measure then you will appreciate the ease with which your script
can adapt to adding or removing outcomes. Even for PANAS, we could try
collapsing two outcomes and see how the model fit changes.

.. code-block:: r

   spec <- list()
   spec[1:10] <- rpf.grm(outcomes = length(PANASItem))

The ``rpf.grm`` function creates a ``rpf.base`` class object that
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
items. If you need to set some starting values to specific values then
you might start with random starting values and then override the rows
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

Here set the label of every parameter in the first row to
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

The custom compute plan is a somewhat more sophisticated version of

.. code-block:: r

   mxComputeIterate(list(
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

   panas1 <- mxRun(panas1)
   panas1$ItemParam$values
   # or
   summary(panas1)

.. [BockAitkin1981] Bock, R. D. & Aitkin, M. (1981). Marginal maximum likelihood estimation of item parameters:
		    Application of an EM algorithm. Psychometrika, 46, 443–459.

.. [Embretson1996] Embretson, S. E. (1996). The new rules of measurement. Psychological Assessment, 8(4), 341-349.

.. [WatsonEtal1988] Watson, D., Clark, L. A., & Tellegen, A. (1988). Development and validation of brief
		    measures of positive and negative affect: The PANAS scales. Journal of Personality
		    and Social Psychology, 54 (6), 1063.
