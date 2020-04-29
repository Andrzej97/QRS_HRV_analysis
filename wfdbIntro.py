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
DETECTION_RANGE = 53

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


def calculate_stats(annotatedX, detectedX):
    fPos = []; fNeg = []; tPos = []
    print('annotated:', len(annotatedX), '/ detected:', len(detectedX))

    for anno in annotatedX:
        inRange = list(filter(lambda x: x >= anno - DETECTION_RANGE and x <= anno + DETECTION_RANGE, detectedX))
        if len(inRange) > 0:
            tPos.append(inRange[0])
        else:
            fNeg.append(anno)

    for det in detectedX:
        inRange = list(filter(lambda x: x >= det - DETECTION_RANGE and x <= det + DETECTION_RANGE, annotatedX))
        if len(inRange) == 0:
            fPos.append(det)
    
    print('tPos:', len(tPos), 'fPos:', len(fPos), 'fNeg: ', len(fNeg))
    return (tPos, fPos, fNeg)


if __name__ == '__main__':
    # freeze_support()
    # download_all_files()
    filepath = 'db/100'
    record = wfdb.rdrecord(filepath, sampto=5000)
    anno = wfdb.rdann(filepath, 'atr', sampto=5000)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    ecg = np.array(signal_ch0)
    peaksR = find_R_peaks(ecg)
    peaksY = list(map(lambda x: x[1], peaksR))
    peaksX = list(map(lambda x: x[0], peaksR))

    annoPeaksX = anno.sample
    tPos, fPos, fNeg = calculate_stats(annoPeaksX, peaksX)
    
    plt.plot(signal_ch0)
    # plt.plot(peaksX, peaksY, 'ro')
    plt.plot(annoPeaksX, ecg[annoPeaksX], 'bo')
    plt.plot(tPos, ecg[tPos], 'go')
    plt.plot(fPos, ecg[fPos], 'yo')
    plt.plot(fNeg, ecg[fNeg], 'rx')

    plt.show()

