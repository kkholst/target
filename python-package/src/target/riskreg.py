# encoding: utf-8
#
# Copyright (c) 2019 Klaus K. Holst.  All rights reserved.

import target as tg
from target import datatype as inp
import numpy as np
import statsmodels.api as sm
from scipy import optimize

class riskreg:
    """Documentation for Riskreg

    """
    model = None
    propensity = None
    mle = None
    propensity_coef = None
    mle_coef = None
    modeltype = None
    estimate = None

    def stageone(self, pr):
        n = pr.size
        par = np.concatenate((self.mle_coef, self.propensity_coef))
        self.model.update(par)

        pp = self.model.pr()
        p0 = pp[:, 1].flatten()
        p1 = pp[:, 2].flatten()
        targ = pp[:, 3].flatten()
        if self.modeltype == 'rd':
            # E(A drho/(Pa(1-Pa))|V) = pr*drho/[p1(1-p1)]
            nom = pr*(1-targ**2)/(p1*(1-p1))
            # E(1/(Pa(1-Pa))|V) =  (1-pr)/[p0(1-p0)] + pr/[p1(1-p1)]
            denom = (1-pr)/(p0*(1-p0)) + pr/(p1*(1-p1))
            omega = nom/denom / (pr*p0*(1-p0))
        else:  # rr
            # E(A pa/(1-Pa) |V) = pr*drho/[p1(1-p1)]
            nom = pr*p1/(1-p1)
            # E(A pa/(1-Pa) |V) = pr*drho/[p1(1-p1)]
            denom = (1-pr)*p0/(1-p0) + pr*p1/(1-p1)
            omega = nom/denom / (pr*(1-p0))

        w = self.model.data(inp.w)
        w = np.multiply(w.flatten(), omega).reshape(n, 1)
        self.model.weights(w)
        self.model.update(par)

    def __init__(self, y, a, **kwargs):
        """Risk regression with binary exposure

        :param y: Response vector (0,1)
        :param a: Exposure vector (0,1)
        :param x1: Design matrix for linear interactions with exposure 'a'
        :param x2: Design matrix for nuisance parameter (odds-product)
        :param x3: Design matrix for propoensity modle
        :returns: Riskreg object
        :rtype: Riskreg

        """

        n = y.size
        one = np.matrix(np.repeat(1.0, n)).transpose()
        x1 = kwargs.get('x1', one)
        x2 = kwargs.get('x2', one)
        x3 = kwargs.get('x3', one)
        w = kwargs.get('weights', one)
        self.modeltype = kwargs.get('model', 'rr')
        self.model = tg.riskregmodel(y, a, x1, x2, x3, w, self.modeltype)
        self.propensity = sm.GLM(endog=a, exog=x3, weights=w, family=sm.families.Binomial())
        self.propensity_coef = self.propensity.fit().params
        self.mle = riskreg_mle(y, a, x2=x2, x1=x1, weights=w, model=self.modeltype)
        self.mle_coef = self.mle['x']

        pr = tg.expit(np.matmul(self.model.data(inp.x3), self.propensity_coef))
        self.stageone(pr)
        alpha0 = self.mle_coef[:x1.shape[1]].reshape(x1.shape[1], 1)

        def obj(alpha):
            return sum(sum(self.model.esteq(alpha, pr))**2)

        op = optimize.minimize(obj, alpha0, method='Nelder-Mead')
        self.estimate = op['x']

    def __repr__(self):
        return "Riskreg. Estimate: " + str(self.estimate)

    def __str__(self):
        return "Riskreg. Estimate: " + str(self.estimate)


def riskreg_mle(y, a, x2, *args, **kwargs):
    one = np.matrix(np.repeat(1.0, len(y))).transpose()
    x1 = kwargs.get('x1', one)
    w = kwargs.get('weights', one)
    model = kwargs.get('model', 'rr')
    m = tg.riskregmodel(y, a, x1, x2, one, w, model)

    def obj(theta):
        m.update(np.matrix(theta))
        return -m.loglik()

    def jac(theta):
        m.update(np.matrix(theta))
        return -m.score().flatten()

    p = x1.shape[1]+x2.shape[1]
    init = kwargs.get('init', np.repeat(0, p))
    op = optimize.minimize(obj, init, method='BFGS', jac=jac)
    if not op['success']:
        op = optimize.minimize(obj, init, method='CG', jac=jac)
    if not op['success']:
        op = optimize.minimize(obj, init, method='Nelder-Mead', jac=jac)
    if not op['success']:
        op = optimize.minimize(obj, init, method='TNC', jac=jac)

    return op
