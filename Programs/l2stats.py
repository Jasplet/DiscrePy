#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 15:55:20 2021

@author: ja17375
"""
from scipy import stats

def ftest(lam2min,ndf,k=2,alpha=0.05):
    """
    returns lambda2 value at 100(1-alpha)% confidence interval
    by default alpha = 0.05 = 95% confidence interval
    following Silver and Chan (1991) [modifications by walsh et al., 2013]
    As we are dealing with traces that have alreayd been passed through SHEBA,
    we do not need to check (or calculate) degrees of freedom as this has alreay
    been done.

    Needed for pair_stack to calc lam2alpha for SKS and SKKS
    """
    F = stats.f.ppf(1-alpha,k,ndf)
    lam2alpha = lam2min * ( 1 + (k/(ndf-k)) * F)
    return lam2alpha