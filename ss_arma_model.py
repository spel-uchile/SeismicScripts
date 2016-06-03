# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=3>

# ARMA example using sunpots data

# <codecell>

import numpy as np
from scipy import stats
import pandas
import matplotlib.pyplot as plt

import statsmodels.api as sm

# <codecell>

from statsmodels.graphics.api import qqplot

# <codecell>

print sm.datasets.sunspots.NOTE

# <codecell>

dta = sm.datasets.sunspots.load_pandas().data

# <codecell>

dta.index = pandas.Index(sm.tsa.datetools.dates_from_range('1700', '2008'))
del dta["YEAR"]

# <codecell>

dta.plot(figsize=(12,8));

# <codecell>

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
fig = sm.graphics.tsa.plot_acf(dta.values.squeeze(), lags=40, ax=ax1)
ax2 = fig.add_subplot(212)
fig = sm.graphics.tsa.plot_pacf(dta, lags=40, ax=ax2)

# <codecell>

arma_mod20 = sm.tsa.ARMA(dta, (2,0)).fit()
print arma_mod20.params

# <codecell>

arma_mod30 = sm.tsa.ARMA(dta, (3,0)).fit()

# <codecell>

print arma_mod20.aic, arma_mod20.bic, arma_mod20.hqic

# <codecell>

print arma_mod30.params

# <codecell>

print arma_mod30.aic, arma_mod30.bic, arma_mod30.hqic

# <markdowncell>

# * Does our model obey the theory?

# <codecell>

sm.stats.durbin_watson(arma_mod30.resid.values)

# <codecell>

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax = arma_mod30.resid.plot(ax=ax);

# <codecell>

resid = arma_mod30.resid

# <codecell>

stats.normaltest(resid)

# <codecell>

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
fig = qqplot(resid, line='q', ax=ax, fit=True)

# <codecell>

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
fig = sm.graphics.tsa.plot_acf(resid.values.squeeze(), lags=40, ax=ax1)
ax2 = fig.add_subplot(212)
fig = sm.graphics.tsa.plot_pacf(resid, lags=40, ax=ax2)

# <codecell>

r,q,p = sm.tsa.acf(resid.values.squeeze(), qstat=True)
data = np.c_[range(1,41), r[1:], q, p]
table = pandas.DataFrame(data, columns=['lag', "AC", "Q", "Prob(>Q)"])
print table.set_index('lag')

# <markdowncell>

# * This indicates a lack of fit.

# <markdowncell>

# * In-sample dynamic prediction. How good does our model do?

# <codecell>

predict_sunspots = arma_mod30.predict('1990', '2012', dynamic=True)
print predict_sunspots

# <codecell>

ax = dta.ix['1950':].plot(figsize=(12,8))
ax = predict_sunspots.plot(ax=ax, style='r--', label='Dynamic Prediction');
ax.legend();
ax.axis((-20.0, 38.0, -4.0, 200.0));

# <codecell>

def mean_forecast_err(y, yhat):
    return y.sub(yhat).mean()

# <codecell>

mean_forecast_err(dta.SUNACTIVITY, predict_sunspots)

# <headingcell level=3>

# Exercise: Can you obtain a better fit for the Sunspots model? (Hint: sm.tsa.AR has a method select_order)

# <headingcell level=3>

# Simulated ARMA(4,1): Model Identification is Difficult

# <codecell>

from statsmodels.tsa.arima_process import arma_generate_sample, ArmaProcess

# <codecell>

np.random.seed(1234)
# include zero-th lag
arparams = np.array([1, .75, -.65, -.55, .9])
maparams = np.array([1, .65])

# <markdowncell>

# * Let's make sure this models is estimable.

# <codecell>

arma_t = ArmaProcess(arparams, maparams)

# <codecell>

arma_t.isinvertible()

# <codecell>

arma_t.isstationary()

# <rawcell>

# * What does this mean?

# <codecell>

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.plot(arma_t.generate_sample(size=50));

# <codecell>

arparams = np.array([1, .35, -.15, .55, .1])
maparams = np.array([1, .65])
arma_t = ArmaProcess(arparams, maparams)
arma_t.isstationary()

# <codecell>

arma_rvs = arma_t.generate_sample(size=500, burnin=250, scale=2.5)

# <codecell>

fig = plt.figure(figsize=(12,8))
ax1 = fig.add_subplot(211)
fig = sm.graphics.tsa.plot_acf(arma_rvs, lags=40, ax=ax1)
ax2 = fig.add_subplot(212)
fig = sm.graphics.tsa.plot_pacf(arma_rvs, lags=40, ax=ax2)

# <rawcell>

# * For mixed ARMA processes the Autocorrelation function is a mixture of exponentials and damped sine waves after (q-p) lags.
# * The partial autocorrelation function is a mixture of exponentials and dampened sine waves after (p-q) lags.

# <codecell>

arma11 = sm.tsa.ARMA(arma_rvs, (1,1)).fit()
resid = arma11.resid
r,q,p = sm.tsa.acf(resid, qstat=True)
data = np.c_[range(1,41), r[1:], q, p]
table = pandas.DataFrame(data, columns=['lag', "AC", "Q", "Prob(>Q)"])
print table.set_index('lag')

# <codecell>

arma41 = sm.tsa.ARMA(arma_rvs, (4,1)).fit()
resid = arma41.resid
r,q,p = sm.tsa.acf(resid, qstat=True)
data = np.c_[range(1,41), r[1:], q, p]
table = pandas.DataFrame(data, columns=['lag', "AC", "Q", "Prob(>Q)"])
print table.set_index('lag')

# <headingcell level=3>

# Exercise: How good of in-sample prediction can you do for another series, say, CPI

# <codecell>

macrodta = sm.datasets.macrodata.load_pandas().data
macrodta.index = pandas.Index(sm.tsa.datetools.dates_from_range('1959Q1', '2009Q3'))
cpi = macrodta["cpi"]

# <headingcell level=4>

# Hint:

# <codecell>

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax = cpi.plot(ax=ax);
ax.legend();

# <rawcell>

# P-value of the unit-root test, resoundly rejects the null of no unit-root.

# <codecell>

print sm.tsa.adfuller(cpi)[1]



