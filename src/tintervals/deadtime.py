# from intervals import intervals2weights
from tintervals.intervals import intervals2weights
import numpy as np


def unc_fft(
    vals1,
    vals2,
    wpm=0.0,
    fpm=0.0,
    wfm=0.0,
    ffm=0.0,
    rwfm=0.0,
    fwfm=0.0,
    step=1.0,
    scale=1.0,
    ext_factor=10,
    return_dict=False,
    return_fft=False,
):
    """Calculate dead time uncertainty using the Fourier transform strategy [1].

    Parameters
    ----------
    vals1 : 2d array
            starting intervals in the form start,stop
    vals2 : 2d array
            ending intervals in the form start, stop
    wpm : float, optional
            white phase modulation noise, Adev normalized, by default 0.
    fpm : float, optional
            flicker phase modulation noise, Adev normalized, by default 0.
    wfm : float, optional
            white frequency modulation noise, Adev normalized, by default 0.
    ffm : float, optional
            flicker frequency modulation noise, Adev normalized, by default 0.
    rwfm : float
            random walk frequency modulation noise, Adev normalized, optional, by default 0.
    fwfm : float, optional
            flicker walk frequency modulation noise, Hadamar dev normalized, by default 0.
    step : float, optional
            minimum time interval in seconds, by default 1
    scale : float or 'day', optional
            scale of the input intervals, by default 1.
            If 'day', it is converted to 86400, for example for intervals given in MJD.
    ext_factor : float, optional
            extension factor (padding) of the FFT, by default 10
    return_dict : bool, optional
            If True, also return a tally of uncertainty by noise type, by default False
    return_fft : bool, optional
            If True, also return FFT arrays, by default False

    Returns
    -------
    unc  : float
            Dead time uncertainty between vals1 and vals2
    par  : dict
            Dictionary with each component of the uncertainty for 'wpm', 'fpm', 'wfm', 'ffm', 'rwfm', 'fwfm'
            Only provided if `return_dict` is True.
    fft_freq, fft_sens, fft_psd : arrays
            Arrays of the frequency, sensitivity function and noise psd used in the calculations.
            Only provided if `return_fft` is True.

    Notes
    -----
    Noise type 'wpm', 'fpm', 'wfm', 'ffm', 'rwfm' corresponds to the usual power-law noise and should be given
    as the corresponding value of the Allan deviation of the noise at 1 s. See for example [2].
    Noise type 'fwfm' corresponds to flicker walk frequency noise and should be given
    as the corresponding value of the Hadamar deviation of the noise at 1 s. See for example [3][4].

    Example
    -------

    >>> maser_noise = {'wfm': 4e-14, 'ffm': 3e-16}
    >>> initial = np.array([[0,1], [2,3]])
    >>> final = np.array([[0,5]])
    >>> ti.deadtime.unc_fft(initial, final, **maser_noise, scale='day')
    2.0026233736533257e-16

    References
    ----------
    .. [1] Grebing et al., Optica, 3, 563-569 (2016)
    .. [2] Dawkins et al., IEEE Trans. Ultrason., Ferroelect., Freq. Cont., 54, 918-925 (2007)
    .. [3] Pizzocaro et al., Nature Physics, 17, 223-227 (2021)
    .. [4] Nemitz et al., Metrologia,  58, 025006 (2021)



    """

    vals1 = np.atleast_2d(vals1)
    vals2 = np.atleast_2d(vals2)

    both = np.concatenate((vals1, vals2))
    start = np.amin(both)
    stop = np.amax(both)

    if scale == "day":
        scale = 86400.0

    t1, w1 = intervals2weights(vals1, step=step / scale, min=start, max=stop, norm=True)
    t2, w2 = intervals2weights(vals2, step=step / scale, min=start, max=stop, norm=True)

    # note timetags should be aligned
    # t = t1
    wt = w2 - w1
    wt = wt - np.mean(wt)  # remove mean to avoid 0 freq components (and leaking in neighbourg bins)

    ft_length = len(wt) * ext_factor

    # note np.fft pad the input array to ft_length
    # ranges for freq > 0
    ft_wt = np.fft.fft(wt, n=ft_length)[1 : ft_length // 2]
    ft_freq = np.fft.fftfreq(n=ft_length, d=step)[1 : ft_length // 2]
    ft_sens = np.abs(ft_wt) ** 2
    df = ft_freq[0]

    tau0 = step
    fH = 0.5 / (tau0)

    # S_y(f) = h_a * f^a
    # then sigma_y(tau)**2 = prefactor * h_a * tau**b
    # inspired by noise_kasdin_01.py from Nils Nemitz
    noises = {"wpm": wpm, "fpm": fpm, "wfm": wfm, "ffm": ffm, "rwfm": rwfm, "fwfm": fwfm}
    prefactor = {
        "wpm": (3.0 / 4.0) * fH / (np.pi**2),
        "fpm": (3.0 * np.euler_gamma + np.log(4.0) + 3.0 * np.log(fH * np.pi * tau0)) / (4.0 * np.pi**2),
        "wfm": 0.5,
        "ffm": 2.0 * np.log(2.0),
        "rwfm": (2.0 / 3.0) * np.pi**2,
        "fwfm": 16.0
        * (np.pi**2)
        * np.log((3.0 / 4.0) * 3.0 ** (11.0 / 16.0))
        / 6.0,  # prefactor for Hadamar var instead of allan
    }
    slope = {"wpm": 2, "fpm": 1, "wfm": 0, "ffm": -1, "rwfm": -2, "fwfm": -3}

    res = {}
    rvar = 0.0
    tot_psd = 0.0

    for noise in noises:
        psd = noises[noise] ** 2 / prefactor[noise] * ft_freq ** slope[noise]

        tot_psd += psd

        var = np.sum(psd * ft_sens * df)
        rvar += var
        res[noise] = np.sqrt(var)

    ret = (np.sqrt(rvar),)
    if return_dict:
        ret += (res,)
    if return_fft:
        ret += (ft_freq, ft_sens, tot_psd)

    if len(ret) == 1:
        return ret[0]
    else:
        return ret


def masercov_fft(
    vals, wpm=0.0, fpm=0.0, wfm=0.0, ffm=0.0, rwfm=0.0, fwfm=0.0, step=1.0, scale=1.0, ext_factor=10, return_dict=False
):
    """Calculate covariance between intervals sampling the same maser noise.
    Note that this is not yet a dead-time covariance, that can be calculated by the law of propagation of uncertainties.

    Parameters
    ----------
    vals : list of 2d array
            list of sampling intervals in the form start,stop
    wpm : float, optional
            white phase modulation noise, Adev normalized, by default 0.
    fpm : float, optional
            flicker phase modulation noise, Adev normalized, by default 0.
    wfm : float, optional
            white frequency modulation noise, Adev normalized, by default 0.
    ffm : float, optional
            flicker frequency modulation noise, Adev normalized, by default 0.
    rwfm : float
            random walk frequency modulation noise, Adev normalized, optional, by default 0.
    fwfm : float, optional
            flicker walk frequency modulation noise, Hadamar dev normalized, by default 0.
    step : float, optional
            minimum time interval in seconds, by default 1
    scale : float or 'day', optional
            scale of the input intervals, by default 1.
            If 'day', it is converted to 86400, for example for intervals given in MJD.
    ext_factor : float, optional
            extension factor (padding) of the FFT, by default 10
    return_dict : bool, optional
            If True, also return a tally of covariance by noise type, by default False

    Returns
    -------
    cov  : 2d array
            Covariance from the maser noise sampled by different intervals.
    par  : dict
            Dictionary each covariance for 'wpm', 'fpm', 'wfm', 'ffm', 'rwfm', 'fwfm'.
            Only provided if `return_dict` is True.

    Notes
    -----
    Noise type 'wpm', 'fpm', 'wfm', 'ffm', 'rwfm' corresponds to the usual power-law noise and should be given
    as the corresponding value of the Allan deviation of the noise at 1 s. See for example [2].
    Noise type 'fwfm' corresponds to flicker walk frequency noise and should be given
    as the corresponding value of the Hadamar deviation of the noise at 1 s. See for example [3][4].


    References
    ----------
    .. [1] Grebing et al., Optica, 3, 563-569 (2016)
    .. [2] Dawkins et al., IEEE Trans. Ultrason., Ferroelect., Freq. Cont., 54, 918-925 (2007)
    .. [3] Pizzocaro et al., Nature Physics, 17, 223-227 (2021)
    .. [4] Nemitz et al., Metrologia,  58, 025006 (2021)

    """

    vals = [np.atleast_2d(v) for v in vals]
    N = len(vals)

    all = np.concatenate(vals)
    start = np.amin(all)
    stop = np.amax(all)

    if scale == "day":
        scale = 86400.0

    # 2d array of weights for each val
    tws = np.array([intervals2weights(v, step=step / scale, min=start, max=stop, norm=True)[1] for v in vals])

    # all t, w should have the same length
    ft_length = tws.shape[1] * ext_factor

    # note np.fft pad the input array to ft_length
    # ranges for freq > 0
    # do not remove mean here - to avoid 0 freq components (and leaking in neighbourg bins)
    # it may mess with the padding
    # TODO: investigate proper windowing
    ft_wts = np.fft.fft(tws, axis=1, n=ft_length)[:, 1 : ft_length // 2]
    ft_freq = np.fft.fftfreq(n=ft_length, d=step)[1 : ft_length // 2]
    # ft_sens = np.abs(ft_wt)**2
    df = ft_freq[0]

    tau0 = step
    fH = 0.5 / (tau0)

    # S_y(f) = h_a * f^a
    # then sigma_y(tau)**2 = prefactor * h_a * tau**b
    # inspired by noise_kasdin_01.py from Nils Nemitz
    noises = {"wpm": wpm, "fpm": fpm, "wfm": wfm, "ffm": ffm, "rwfm": rwfm, "fwfm": fwfm}
    prefactor = {
        "wpm": (3.0 / 4.0) * fH / (np.pi**2),
        "fpm": (3.0 * np.euler_gamma + np.log(4.0) + 3.0 * np.log(fH * np.pi * tau0)) / (4.0 * np.pi**2),
        "wfm": 0.5,
        "ffm": 2.0 * np.log(2.0),
        "rwfm": (2.0 / 3.0) * np.pi**2,
        "fwfm": 16.0
        * (np.pi**2)
        * np.log((3.0 / 4.0) * 3.0 ** (11.0 / 16.0))
        / 6.0,  # prefactor for Hadamar var instead of allan
    }
    slope = {"wpm": 2, "fpm": 1, "wfm": 0, "ffm": -1, "rwfm": -2, "fwfm": -3}

    res = {}
    rcov = np.zeros((N, N))
    tot_psd = 0.0

    for noise in noises:
        psd = noises[noise] ** 2 / prefactor[noise] * ft_freq ** slope[noise]
        tot_psd += psd

        # numpy-esque math
        # ft_wts is a N x ft_length array, with each fft for each row

        # ft_sens is a N x N x ft_length array
        # with maser noise sensitivity
        ft_sens = np.real(ft_wts * np.conj(ft_wts)[:, np.newaxis])

        # integrate psd and ft_sens over last axis
        # to get the NxN covariance matrix
        cov = np.sum(psd * ft_sens * df, axis=-1)

        rcov += cov
        res[noise] = cov

    ret = (rcov,)
    if return_dict:
        ret += (res,)

    if len(ret) == 1:
        return ret[0]
    else:
        return ret


def cov_fft(
    vals1,
    vals2,
    wpm=0.0,
    fpm=0.0,
    wfm=0.0,
    ffm=0.0,
    rwfm=0.0,
    fwfm=0.0,
    step=1.0,
    scale=1.0,
    ext_factor=10,
    return_dict=False,
):
    """Calculate covariance for pairs of intervals sampling the same maser noise.
    From an input of two list of intervals of length N, returns a NxN covariance matrix.

    Parameters
    ----------
    vals1 : list of 2d array
            list of sampling intervals in the form start,stop
    vals2 : list of 2d array (same length as vals1)
            list of sampling intervals in the form start,stop
    wpm : float, optional
            white phase modulation noise, Adev normalized, by default 0.
    fpm : float, optional
            flicker phase modulation noise, Adev normalized, by default 0.
    wfm : float, optional
            white frequency modulation noise, Adev normalized, by default 0.
    ffm : float, optional
            flicker frequency modulation noise, Adev normalized, by default 0.
    rwfm : float
            random walk frequency modulation noise, Adev normalized, optional, by default 0.
    fwfm : float, optional
            flicker walk frequency modulation noise, Hadamar dev normalized, by default 0.
    step : float, optional
            minimum time interval in seconds, by default 1
    scale : float or 'day', optional
            scale of the input intervals, by default 1.
            If 'day', it is converted to 86400, for example for intervals given in MJD.
    ext_factor : float, optional
            extension factor (padding) of the FFT, by default 10
    return_dict : bool, optional
            If True, also return a tally of covariance by noise type, by default False

    Returns
    -------
    cov  : 2d array
            Covariance of the dead time extrapolation.
    par  : dict
            Dictionary each covariance for 'wpm', 'fpm', 'wfm', 'ffm', 'rwfm', 'fwfm'.
            Only provided if `return_dict` is True.

    Notes
    -----
    Noise type 'wpm', 'fpm', 'wfm', 'ffm', 'rwfm' corresponds to the usual power-law noise and should be given
    as the corresponding value of the Allan deviation of the noise at 1 s. See for example [2].
    Noise type 'fwfm' corresponds to flicker walk frequency noise and should be given
    as the corresponding value of the Hadamar deviation of the noise at 1 s. See for example [3][4].


    References
    ----------
    .. [1] Grebing et al., Optica, 3, 563-569 (2016)
    .. [2] Dawkins et al., IEEE Trans. Ultrason., Ferroelect., Freq. Cont., 54, 918-925 (2007)
    .. [3] Pizzocaro et al., Nature Physics, 17, 223-227 (2021)
    .. [4] Nemitz et al., Metrologia,  58, 025006 (2021)

    """

    if len(vals1) != len(vals2):
        raise ValueError("vals1 and vals2 need to have the same length.")

    cv, noise_cv = masercov_fft(
        vals1 + vals2,
        wpm=wpm,
        fpm=fpm,
        wfm=wfm,
        ffm=ffm,
        rwfm=rwfm,
        fwfm=fwfm,
        step=step,
        scale=scale,
        ext_factor=ext_factor,
        return_dict=True,
    )

    N = len(vals1)
    A = np.matrix(np.hstack((np.eye(N), -np.eye(N))))

    rcov = 0.0
    res = {}

    for ntype in noise_cv.keys():
        C = np.matrix(noise_cv[ntype])
        cav = np.array(A * C * A.T)
        rcov += cav
        res[ntype] = cav

    ret = (rcov,)
    if return_dict:
        ret += (res,)

    if len(ret) == 1:
        return ret[0]
    else:
        return ret


if __name__ == "__main__":
    from matplotlib.pyplot import *

    close("all")
    ion()

    import statsmodels.stats.moment_helpers as help

    vals1 = np.array([[3, 4]])  # np.array([[0.2,0.6],[0.8,0.9]])
    vals2 = np.array([[0, 1]])
    noise_dict = {"wpm": 1.5e-13, "wfm": 3.5e-14, "ffm": 5.5e-16, "fwfm": 0 * 2.5e-22}

    unc, par, ft_freq, ft_sens, psd = unc_fft(
        vals1, vals2, **noise_dict, step=100, scale="day", return_dict=True, return_fft=True, ext_factor=10
    )

    print("Uncertainty = ", unc)
    print("Tally = ", par)

    cov, noise_cov = masercov_fft(
        [vals1, vals2, np.array([[3, 4]])], **noise_dict, step=100, scale="day", return_dict=True, ext_factor=10
    )

    covunc = np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1])

    print("Maser Covariance = ", cov)
    print("Maser Correlation = ", help.cov2corr(cov))
    print("Dead time uncertainty from covariance = ", covunc)

    # target intervals = "GPS"
    supj = [np.array([[0, 20]]), np.array([[0, 20]]), np.array([[0, 20]]), np.array([[30, 40]])]

    # measured intervals = "clock uptime"
    supi = [np.array([[5, 15]]), np.array([[5, 16]]), np.array([[0, 5], [15, 20]]), np.array([[31, 36]])]

    cov2, noise_cov2 = cov_fft(supj, supi, **noise_dict, step=100, scale="day", return_dict=True, ext_factor=10)

    print("Dead uncs 2 = ", np.sqrt(np.diag(cov2)))
    print("Covariance2 = ", cov2)
    print("Correlation2 = ", help.cov2corr(cov2))

    # fig, ax1 = subplots()
    # ax2  =ax1.twinx()

    # ax1.loglog(ft_freq, ft_sens, 'b-', label='Sensitivity')
    # ax2.loglog(ft_freq, psd, 'r.', label='Noise')

    # ax1.set_xlabel('Freq /Hz')
    # ax1.set_ylabel('Sesnitivity')
    # ax2.set_ylabel('PSD /(1/Hz)')

    # lines, labels = ax1.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # ax2.legend(lines + lines2, labels + labels2, loc=0)

    # fig.tight_layout()
