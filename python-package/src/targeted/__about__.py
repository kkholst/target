# coding: utf-8

__name__        = "targeted"
__version__     = "0.0.26"
__license__     = "Apache Software License"
__description__ = """
**targeted** provides a number of methods for semi-parametric
estimation.  The library also contains implementations of various
parametric models (including different discrete choice models) and
model diagnostics tools.

The implemention currently includes
- **Risk regression models* with binary exposure
  (Richardson et al. (2017) <doi:10.1080/01621459.2016.1192546>)
- **Augmented Inverse Probability Weighted** estimators for missing data and causal inference
  (Bang and Robins (2005) <doi:10.1111/j.1541-0420.2005.00377.x>)
- Model diagnostics based on **cumulative residuals** methods
- Efficient weighted **Pooled Adjacent Violator Algorithms**
- **Nested multinomial logit** models
"""
__author__      = u"Klaus Kähler Holst"
__email__       = "klaus.holst@maersk.com"
__url__         = "https://targetlib.org/python/"
__copyright__   = u"Copyright 2019-2020, Klaus Kähler Holst"
__summary__     = "Targeted inference"
__keywords__    = "semi-parametric inference, double robust estimator,\
 model diagnostics, cumulative residuals, generalized linear models"

__all__ = [
    "__name__", "__summary__", "__description__",
    "__url__", "__version__", "__author__",
    "__email__", "__license__", "__copyright__",
]
