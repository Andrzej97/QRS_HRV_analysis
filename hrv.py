import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
from scipy.integrate import trapz

MAX_RR_DURATION = 2     # seconds


def time_domain(rrs, signal_freq):
    output = {}

    rr_s = list(map(lambda x: x[0] / signal_freq, rrs))

    nn_int = np.diff(rr_s)  # intervals in [s] between beats
    nn_int = list(filter(lambda x: x <= MAX_RR_DURATION, nn_int))

    output['avg_nn_int'] = np.mean(nn_int)
    output['max_nn_int'] = np.max(nn_int)
    output['min_nn_int'] = np.min(nn_int)
    output['range_nn_int'] = np.max(nn_int) - np.min(nn_int)
    output['SDNN'] = np.std(nn_int)

    successive_differences = np.diff(nn_int)
    output['RMSSD'] = np.sqrt(np.mean(np.array(successive_differences) ** 2))  # root mean square of successive differences
    output['SDSD'] = np.std(successive_differences)  # standard deviation of successive differences
    output['NN50'] = np.sum(successive_differences > 0.050)   # the number of pairs of successive NNs that differ by more than 50 ms
    output['NN20'] = np.sum(successive_differences > 0.020)   # the number of pairs of successive NNs that differ by more than 20 ms
    output['pNN50'] = output['NN50'] / len(nn_int)
    output['pNN20'] = output['NN20'] / len(nn_int)

    output['avg_hr'] = 60 / output['avg_nn_int']
    output['max_hr'] = 60 / output['min_nn_int']
    output['min_hr'] = 60 / output['max_nn_int']

    print(output)
    return output


def geometric(rrs, signal_freq):
    rr_s = list(map(lambda x: x[0] / signal_freq, rrs))
    nn_int = np.diff(rr_s)  # intervals in [s] between beats
    nn_int = list(filter(lambda x: x <= MAX_RR_DURATION, nn_int))

    plt.title('Sample density distribution of NN interval durations')
    plt.xlabel('NN intervals duration [s]')
    plt.ylabel('Interval duration count')
    plt.hist(nn_int, bins=30)
    plt.show()

    successive_differences = np.diff(nn_int)
    plt.title('Sample density distribution of differences between adjacent NN intervals')
    plt.xlabel('Adjacent NN intervals duration [s]')
    plt.ylabel('Adjacent interval duration count')
    plt.hist(successive_differences, bins=30)
    plt.show()

    plt.title('Poincare plot')
    plt.xlabel('NN_n [s]')
    plt.ylabel('NN_n+1 [s]')
    plt.scatter(nn_int[:-1], nn_int[1:])
    plt.show()


def frequency_domain(rrs, signal_freq):    
    '''
    http://web.mit.edu/~gari/www/papers/GDCliffordThesis.pdf
    ULF: 0.0001-0.003Hz
    VLF: 0.003-0.04Hz 
    LF: 0.04-0.15Hz 
    HF: 0.15-0.4Hz
    '''

    rr_s = list(map(lambda x: x[0] / signal_freq, rrs))
    nn_int = np.diff(rr_s)  # intervals in [s] between beats
    nn_int = list(filter(lambda x: x <= MAX_RR_DURATION, nn_int))

    x = np.cumsum(nn_int)
    f = interp1d(x, nn_int, kind='cubic')

    fs = 4.0
    steps = 1 / fs

    xx = np.arange(1, np.max(x), steps)
    rr_interpolated = f(xx)


    # based on https://www.kaggle.com/stetelepta/exploring-heart-rate-variability-using-python
    # Estimate the spectral density using Welch's method
    fxx, pxx = signal.welch(x=rr_interpolated, fs=fs)
    cond_vlf = (fxx >= 0.003) & (fxx < 0.04)
    cond_lf = (fxx >= 0.04) & (fxx < 0.15)
    cond_hf = (fxx >= 0.15) & (fxx < 0.4)
    
    # calculate power in each band by integrating the spectral density 
    vlf = trapz(pxx[cond_vlf], fxx[cond_vlf])
    lf = trapz(pxx[cond_lf], fxx[cond_lf])
    hf = trapz(pxx[cond_hf], fxx[cond_hf])
    
    # sum these up to get total power
    total_power = vlf + lf + hf

    # find which frequency has the most power in each band
    peak_vlf = fxx[cond_vlf][np.argmax(pxx[cond_vlf])]
    peak_lf = fxx[cond_lf][np.argmax(pxx[cond_lf])]
    peak_hf = fxx[cond_hf][np.argmax(pxx[cond_hf])]

    lf_nu = 100 * lf / (lf + hf)
    hf_nu = 100 * hf / (lf + hf)
    
    output = {}
    output['Power VLF (ms2)'] = vlf
    output['Power LF (ms2)'] = lf
    output['Power HF (ms2)'] = hf   
    output['Power Total (ms2)'] = total_power

    output['LF/HF'] = (lf/hf)
    output['Peak VLF (Hz)'] = peak_vlf
    output['Peak LF (Hz)'] = peak_lf
    output['Peak HF (Hz)'] = peak_hf

    output['Fraction LF (nu)'] = lf_nu
    output['Fraction HF (nu)'] = hf_nu


    plt.plot(fxx, pxx)
    plt.title("Welch's periodogram")

    psd_f = interp1d(fxx, pxx)
    x_vlf = np.linspace(0.003, 0.04, 100)
    x_lf = np.linspace(0.04, 0.15, 100)
    x_hf = np.linspace(0.15, 0.4, 100)

    plt.gca().fill_between(x_vlf, psd_f(x_vlf), alpha=0.2, label="VLF")
    plt.gca().fill_between(x_lf, psd_f(x_lf), alpha=0.2, label="LF")
    plt.gca().fill_between(x_hf, psd_f(x_hf), alpha=0.2, label="HF")

    plt.gca().set_xlim(0, 0.5)
    plt.gca().set_ylim(0)
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Density")
    plt.legend()
    plt.show()

    return output


def non_linear(rrs, signal_freq):
    pass