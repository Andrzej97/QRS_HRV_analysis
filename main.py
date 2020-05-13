from multiprocessing.dummy import freeze_support
import wfdb
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

FILEPATH = 'db/107'
DESTINATION_PATH = os.getcwd() + '/db'
DB = 'mitdb'
REF_SAMPLES = 20
SEARCH_SAMPLES = 72
THRESHOLD = 0.48
DETECTION_X_RANGE = 53
DETECTION_Y_RANGE = 0.1
R_SYMBOLS = ['N', 'R', 'J', 'A', 'E', 'j', '/', 'Q']


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


def calculate_stats(signal_ch0, annotated_x, detected_x):
    f_pos = []
    f_neg = []
    t_pos = []
    print('annotated:', len(annotated_x), '/ detected:', len(detected_x))

    for anno in annotated_x:
        in_range = list(filter(lambda x: anno - DETECTION_X_RANGE <= x <= anno + DETECTION_X_RANGE, detected_x))
        in_y_range_found= False
        if len(in_range) > 0:
            for x in in_range:
                if abs(signal_ch0[anno] - signal_ch0[x]) <= DETECTION_Y_RANGE:
                    t_pos.append(x)
                    in_y_range_found = True
                    break
            if not in_y_range_found:
                f_neg.append(anno)
        else:
            f_neg.append(anno)

    for det in detected_x:
        in_range = list(filter(lambda x: det - DETECTION_X_RANGE <= x <= det + DETECTION_X_RANGE, annotated_x))
        if len(in_range) == 0:
            f_pos.append(det)
        else:
            in_y_range_found = False
            for x in in_range:
                if abs(signal_ch0[det] - signal_ch0[x]) <= DETECTION_Y_RANGE:
                    in_y_range_found = True
                    break
            if not in_y_range_found:
                f_pos.append(det)
    print('t_pos:', len(t_pos), 'f_pos:', len(f_pos), 'f_neg: ', len(f_neg))
    return t_pos, f_pos, f_neg

def get_r_samples(ann):
    return list(filter(lambda x: x[1] in R_SYMBOLS, zip(ann.sample, ann.symbol)))

def get_plot_data():
    record = wfdb.rdrecord(FILEPATH, sampto=5000)
    ann = wfdb.rdann(FILEPATH, 'atr', sampto=5000)
    annotations = get_r_samples(ann)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    ecg = np.array(signal_ch0)
    peaks_r = find_R_peaks(ecg)
    peaks_y = list(map(lambda x: x[1], peaks_r))
    peaks_x = list(map(lambda x: x[0], peaks_r))
    anno_peaks_x = list(map(lambda x: x[0], annotations))
    return signal_ch0, peaks_r, peaks_y, peaks_x, anno_peaks_x, ecg


def plot_data(subplot):
    signal_ch0, peaks_r, peaks_y, peaks_x, anno_peaks_x, ecg = get_plot_data()
    t_pos, f_pos, f_neg = calculate_stats(signal_ch0, anno_peaks_x, peaks_x)
    subplot.cla()
    subplot.plot(signal_ch0)
    subplot.plot(peaks_x, peaks_y, 'ro')
    subplot.plot(anno_peaks_x, ecg[anno_peaks_x], 'bo')
    subplot.plot(t_pos, ecg[t_pos], 'go')
    subplot.plot(f_pos, ecg[f_pos], 'yo')
    subplot.plot(f_neg, ecg[f_neg], 'rx')


def sliders_on_changed(val):
    global REF_SAMPLES, SEARCH_SAMPLES, THRESHOLD, DETECTION_X_RANGE
    REF_SAMPLES = int(ref_slider.val)
    SEARCH_SAMPLES = int(search_slider.val)
    THRESHOLD = threshold_slider.val
    DETECTION_X_RANGE = int(detection_slider.val)
    plot_data(ax)


if __name__ == '__main__':
    # freeze_support()
    # download_all_files()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom=0.3)

    plot_data(ax)

    ref_samples_ax = plt.axes([0.25, 0.15, 0.65, 0.03])
    ref_slider = Slider(ref_samples_ax, 'ref samples', 1, 100, valinit=REF_SAMPLES, valfmt='%0.0f')

    # Draw another slider
    search_saples_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
    search_slider = Slider(search_saples_ax, 'search samples', 1, 100, valinit=SEARCH_SAMPLES, valfmt='%0.0f')

    threshold_ax = plt.axes([0.25, 0.20, 0.65, 0.03])
    threshold_slider = Slider(threshold_ax, 'threshold', 0.01, 1.0, valinit=THRESHOLD)

    # Draw another slider
    detections_samples_ax = plt.axes([0.25, 0.05, 0.65, 0.03])
    detection_slider = Slider(detections_samples_ax, 'detection_samples', 1.0, 100.0, valinit=DETECTION_X_RANGE,
                              valfmt='%0.0f')

    ref_slider.on_changed(sliders_on_changed)
    search_slider.on_changed(sliders_on_changed)
    threshold_slider.on_changed(sliders_on_changed)
    detection_slider.on_changed(sliders_on_changed)
    plt.show()
