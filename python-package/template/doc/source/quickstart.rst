.. _code_directive:

Quickstart
====================


Installation
--------------------

Install with ``pip``:

.. code-block:: console

    pip install targeted

.. note:: Presently binary wheels are available for Linux
	  systems and Mac OS X running Python>=3.6. Windows
	  installations must be built from source.


Risk regression
--------------------

As an illustration data is loaded from the package

.. ipython:: python

   import targeted as tg

   d = tg.getdata() # returns a Pandas DataFrame
   print(d.head())

Here ``y`` is the binary response variable, ``a`` binary the exposure,
and ``x``, ``z`` covariates.

Next we estimate the *relative risk* of the exposed vs the non-exposed

.. ipython:: python
   :okwarning:

   import numpy as np
   from patsy import dmatrices

   y, X = dmatrices('y ~ x+z', d)
   a = d['a']
   ones = np.ones((y.size,1))
   tg.riskreg(y=y, a=a, x1=ones, x2=X, type='rr')

Or using the formula syntax

.. ipython:: python

   from targeted.formula import riskreg

   riskreg(d, 'y~a', type='rr')

And adjusting for covariates

.. ipython:: python

   from targeted.formula import riskreg

   riskreg(d, 'y~a', nuisance='~x+z', type='rr')

.. ipython:: python

   from targeted.formula import riskreg

   riskreg(d, 'y~a', nuisance='~x+z', propensity='~1', type='rr')
