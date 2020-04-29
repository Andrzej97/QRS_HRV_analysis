from multiprocessing.dummy import freeze_support
import wfdb
import os
import matplotlib.pyplot as plt
import numpy as np

DESTINATION_PATH = os.getcwd() + '/db'
DB = 'mitdb'
REF_SAMPLES = 40
SEARCH_SAMPLES = 72
THRESHOLD = 0.48

def download_all_files():
    wfdb.dl_database(DB, DESTINATION_PATH, records=['100', '107', '108', '200', '203', '207', '222', '233'])

def find_R_peaks(ecg):
    ref_samples = list(ecg[0:REF_SAMPLES])
    ref_samples_sum = sum(ref_samples)
    search_samples_left = 0
    max_diff_sample = (0, 0)
    r_peaks = []
    for sample in enumerate(ecg):
        ref_signal = ref_samples_sum / REF_SAMPLES
        signal = abs(sample[1] - ref_signal)
        if signal > THRESHOLD:
            if signal > max_diff_sample[1]:
                max_diff_sample = sample
            if search_samples_left <= 0:
                search_samples_left = SEARCH_SAMPLES
        if search_samples_left == 1:
            r_peaks.append(max_diff_sample)
            max_diff_sample = (0, 0)
        search_samples_left -= 1

        ref_samples_sum -= ref_samples.pop(0)
        ref_samples_sum += sample[1]
        ref_samples.append(sample[1])
    return r_peaks


if __name__ == '__main__':
    # freeze_support()
    # download_all_files()
    record = wfdb.rdrecord('db/100', sampto=5000)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    ecg = np.array(signal_ch0)
    peaksR = find_R_peaks(ecg)
    peaksY = list(map(lambda x: x[1], peaksR))
    peaksX = list(map(lambda x: x[0], peaksR))
    plt.plot(peaksX, peaksY, 'ro')
    plt.plot(signal_ch0)
    plt.show()

