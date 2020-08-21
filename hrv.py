import numpy as np
import matplotlib.pyplot as plt

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
    pass

def non_linear(rrs, signal_freq):
    pass