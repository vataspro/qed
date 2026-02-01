import numpy as np
import matplotlib.pyplot as plt
from glob import glob

files = glob("../data/*.txt")

def ms_autocorr_time(x, c=5.0):
    """
    Compute the Madras–Sokal integrated autocorrelation time and its analytic error.
    """
    x = np.asarray(x, np.float64)
    n = x.size
    if n < 8:
        return np.nan, np.nan

    x = x - np.mean(x)
    var = np.var(x)
    if not np.isfinite(var) or var <= 1e-15:
        return np.nan, np.nan

    # FFT-based autocovariance
    nfft = 1 << (2 * n - 1).bit_length()
    fx = np.fft.rfft(x, nfft)
    acov = np.fft.irfft(fx * np.conjugate(fx), nfft)[:n]
    rho = np.real(acov / acov[0])
    rho[np.isnan(rho)] = 0.0
    rho = np.clip(rho, -1.0, 1.0)

    # Cut negative tail
    neg = np.where(rho < 0)[0]
    if len(neg) > 0:
        rho[neg[0] :] = 0.0

    tau = 0.5
    for _ in range(1000):
        W = int(max(1, np.floor(c * tau)))
        if W >= n:
            break
        new_tau = 0.5 + np.sum(rho[1 : W + 1])
        if abs(new_tau - tau) < 1e-5:
            tau = new_tau
            break
        tau = new_tau

    if not np.isfinite(tau) or tau <= 0.0:
        return np.nan, np.nan

    tau_err = tau * np.sqrt((4 * (2 * W + 1)) / n)
    return float(tau), float(tau_err)

# Initialise plot lists
x = []
y1 = []
errors1 = []
y2 = []
errors2 = []

# Loop over files and measure observables
for file_ in files:

    beta = file_.split("beta")[1].split(".txt")[0]
    x.append(float(beta))

    plaq, mon = np.loadtxt(file_).T
    y1.append(np.mean(plaq))
    y2.append(np.mean(mon))

    
    for err_list, obs in zip([errors1, errors2], [plaq, mon]):

        tau, _ = ms_autocorr_time(plaq)
        meas = []

        for bs_num in range(200):
            idx = np.random.randint(len(obs), size=int(len(obs)//tau))
            meas.append(np.mean(obs[idx]))
        err_list.append(np.std(meas))


fig, ax = plt.subplots(2, 1)
ax[0].errorbar(x, y1, errors1, ls="none", marker='o')
ax[1].errorbar(x, y2, errors2, ls="none", marker='o')

ax[0].set_ylabel(r"$\langle P \rangle$")
ax[1].set_ylabel(r"$\langle M \rangle$")
ax[1].set_xlabel(r"$\beta$")
plt.show()

