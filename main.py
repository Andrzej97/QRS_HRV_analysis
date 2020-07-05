from multiprocessing.dummy import freeze_support
import wfdb
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider
import math

FILEPATH = 'db/217'
DESTINATION_PATH = os.getcwd() + '/db'
DB = 'mitdb'
REF_SAMPLES = 40
SEARCH_SAMPLES = 72
THRESHOLD = 0.48
DETECTION_X_RANGE = 53
DETECTION_Y_RANGE = 0.5
R_SYMBOLS = ['N']




def download_all_files():
    wfdb.dl_database(DB, DESTINATION_PATH)


def find_R_peaks(ecg):
    ref_samples = list(ecg[0:REF_SAMPLES])
    ref_samples_sum = sum(ref_samples)
    search_samples_left = 0
    max_signal = 0
    max_diff_sample = (0, 0)
    r_peaks = []
    for sample in enumerate(ecg):
        if sample[0] < REF_SAMPLES:
            continue
        ref_signal = ref_samples_sum / REF_SAMPLES
        signal = sample[1] - ref_signal
        if signal > THRESHOLD:
            if signal > max_signal:
                max_signal = signal
                max_diff_sample = sample
            if search_samples_left <= 0:
                search_samples_left = SEARCH_SAMPLES
        if search_samples_left == 1:
            r_peaks.append(max_diff_sample)
            max_signal = 0
            max_diff_sample = (0, 0)
        search_samples_left -= 1

        ref_samples_sum -= ref_samples.pop(0)
        ref_samples_sum += sample[1]
        ref_samples.append(sample[1])
    return r_peaks


def find_R_peaks10(ecg):
    ref_samples = list(ecg[0:REF_SAMPLES])
    ref_samples_sum = sum(ref_samples)
    search_samples_left = 0
    max_signal = 0
    max_diff_sample = (0, 0)
    r_peaks = []
    for sample in enumerate(ecg):
        if sample[0] < REF_SAMPLES:
            continue
        ref_signal = ref_samples_sum / REF_SAMPLES
        signal = sample[1] - ref_signal
        if signal > THRESHOLD:
            if signal > max_signal:
                max_signal = signal
                max_diff_sample = sample
            if search_samples_left <= 0:
                search_samples_left = SEARCH_SAMPLES
        if search_samples_left == 1:
            r_peaks.append(max_diff_sample)
            if 10 == len(r_peaks):
                return r_peaks
            max_signal = 0
            max_diff_sample = (0, 0)
        search_samples_left -= 1

        ref_samples_sum -= ref_samples.pop(0)
        ref_samples_sum += sample[1]
        ref_samples.append(sample[1])
    return r_peaks


weights = [0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05]
PEAK_THRESHOLD = 0.3
def find_R_peaks_weights(ecg):
    old_r_peaks = find_R_peaks(ecg)
    last_ten_peaks = old_r_peaks[:10]
    last_ten_peaks.reverse()
    r_peaks = []
    max_signal = 0
    max_diff_sample = (0,0)
    search_samples_left = 0
    expected_peak = get_expected_peak(last_ten_peaks)
    for sample in enumerate(ecg):
        if abs(sample[1] - expected_peak) < PEAK_THRESHOLD:
            if sample[1] > max_signal:
                max_signal = sample[1]
                max_diff_sample = sample
            if search_samples_left <= 0:
                search_samples_left = SEARCH_SAMPLES
        if search_samples_left == 1:
            r_peaks.append(max_diff_sample)
            last_ten_peaks.pop()
            last_ten_peaks.insert(0,max_diff_sample)
            # print(expected_peak, max_diff_sample)
            max_signal = 0
            max_diff_sample = (0, 0)
            expected_peak = get_expected_peak(last_ten_peaks)
        search_samples_left -= 1

    return r_peaks


def find_R_peaks_weights2(ecg, weight, threshold):
    old_r_peaks = find_R_peaks10(ecg)
    last_ten_peaks = old_r_peaks[:10]
    last_ten_peaks.reverse()
    r_peaks = []
    max_signal = 0
    max_diff_sample = (0,0)
    search_samples_left = 0
    last_peak_sample = -72
    expected_peak = get_expected_peak(last_ten_peaks)
    for sample in enumerate(ecg):
        if abs(sample[1] - expected_peak) < (threshold * expected_peak):
            if sample[1] > max_signal:
                max_signal = sample[1]
                max_diff_sample = sample
            if search_samples_left <= 0:
                search_samples_left = SEARCH_SAMPLES
        if search_samples_left == 1:
            if (sample[0] - last_peak_sample) >= 72:
                r_peaks.append(max_diff_sample)
                last_peak_sample = sample[0]
            else:
                print('???')
            last_ten_peaks.pop()
            last_ten_peaks.insert(0,max_diff_sample)
            # print(expected_peak, max_diff_sample)
            max_signal = 0
            expected_peak = (1 - weight) * expected_peak + weight * max_diff_sample[1]
            max_diff_sample = (0, 0)
        search_samples_left -= 1

    return r_peaks


def get_expected_peak(last_ten):
    expected = 0.0
    for i in range(len(last_ten)):
        expected += last_ten[i][1] * weights[i]

    return expected/sum(weights)


def correct_r_peaks(peaks, ecg):
    corrected_peaks = []
    for p in peaks:
        fr = max(0, p - 10)
        to = min(len(ecg) - 1, p + 10)
        sub = ecg[fr : to]
        m = fr + max(range(len(sub)), key=sub.__getitem__)
        # m = max(sub, key=lambda x: ecg[x])
        corrected_peaks.append(m)


    return corrected_peaks


def get_shift(current_window):
    max_elem_idx = max(range(len(current_window)), key=current_window.__getitem__)
    return len(current_window) - max_elem_idx - 1


def ff(ecg):
    # Running median elements
    N = 3 
    Nd = 2

    # Preprocessing variables
    padding = max(N, Nd, 10)
    i = 0
    current_window = []
    derivated_window = []

    # Threshold stage constants
    qrs_interval = 21 # 60ms
    rr_min = 72 # 200ms
    p_threshold = 0.7 * 360 / 128 + 4.7

    # Threshold stage variables
    counter = 0
    state = 1
    max_peak = (-1, -1, -1) # (index, value, shift)
    r_peaks = []
    r_peaks_pos = []
    th = 0

    for value in ecg: 

        # collecting initial data
        if i < padding:
            current_window.append(value)
            derivated_window.append(value)
            i += 1
            continue

        # Preprocessing
        current_window.append(value)
        derivated_window.append(value - current_window[-Nd])

        integrated_value = np.sum(derivated_window[-1:-1-N:-1]) * 1 / (N - 1)
        processed_value = integrated_value ** 2

        # Analysis
        if state == 1:
            counter += 1
            if processed_value > max_peak[1]:
                max_peak = (i, processed_value, get_shift(current_window))

            if counter > rr_min + qrs_interval:
                counter = i - max_peak[0]
                state = 2

                r_peaks.append(max_peak[1])
                r_peaks_pos.append(max_peak[0] - max_peak[2])

        elif state == 2:
            counter += 1
            if counter > rr_min:
                th = np.mean(r_peaks)
                state = 3

        elif state == 3:
            if processed_value > th:
                counter = 0
                state = 1
                max_peak = (-1, -1)
            else:
                th = th * math.exp(-p_threshold / 360)

        current_window.pop(0)
        derivated_window.pop(0)
        i += 1

    # plt.plot(ecg)
    # plt.plot(r_peaks_pos, ecg[r_peaks_pos], 'rx')

    return list(map(lambda x: (x, ecg[x]), r_peaks_pos))



def ff2(ecg):
    N = 8
    Nd = N - 1

    y0 = []
    for i in enumerate(ecg):
        if i[0] < Nd:
            y0.append(i[1])
            continue

        y0.append(i[1] - ecg[i[0] - Nd])
        # y0.append(1 - math.pow(i[1], -Nd))


    y1 = []
    for i in enumerate(y0):
        if i[0] - N + 1 < 0:
            y1.append(i[1])
            continue

        s = 0
        for q in range(N):
            s += y0[i[0] - q]

        x = 1 / (N - 1) * s
        y1.append(x)

    y = []
    for i in y1:
        y.append(i * i)

    qrs_interval = 21 # 60ms
    rr_min = 72 # 200ms
    counter = 0
    state = 1
    max_peak = (-1, -1) # (index, value)
    r_peaks = []
    r_peaks_pos = []

    th = 0
    pth = 0.7 * 360 / 128 + 4.7
    for idx, val in enumerate(ecg):
        if state == 1:
            counter += 1
            if val > max_peak[1]:
                max_peak = (idx, val)

            if counter > rr_min + qrs_interval:
                counter = idx - max_peak[0]
                state = 2
                r_peaks.append(max_peak[1])
                r_peaks_pos.append(max_peak[0])

        elif state == 2:
            counter += 1
            if counter > rr_min:
                th = np.mean(r_peaks)
                state = 3

        elif state == 3:
            # plt.plot(idx, th, 'bx')
            if val > th:
                counter = 0
                state = 1
                max_peak = (-1, -1)

            else:
                th = th * math.exp(-pth / 360)
                # th = 0


    # plt.plot(y0)
    # plt.plot(y1)

    # plt.plot(y)
    # r_peaks_pos = list(map(lambda x: x - N))
    y = np.array(y)
    # y0 = np.array(y0)
    y1 = np.array(y1)
    plt.plot(ecg, 'g-')
    # plt.plot(y0)

    r_peaks_pos = correct_r_peaks(r_peaks_pos, ecg)

    plt.plot(r_peaks_pos, ecg[r_peaks_pos], 'rx')
    # plt.plot(r_peaks_pos, y[r_peaks_pos], 'cx')

    # print(len(y0))
    # print(len(y1))
    # print(len(y))
    # print(len(ecg))
    # print(r_peaks)
    return list(map(lambda x: (x, y[x]), r_peaks_pos))


def calculate_stats(signal_ch0, annotated_x, detected_x):
    f_pos = []
    f_neg = []
    t_pos = []
    print(annotated_x)
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

def calculate_stats_for_tests(signal_ch0, annotated_x, detected_x):
    f_neg = 0
    t_pos = 0
    print('annotated:', len(annotated_x), '/ detected:', len(detected_x))

    for anno in annotated_x:
        in_range = list(filter(lambda x: anno - DETECTION_X_RANGE <= x <= anno + DETECTION_X_RANGE, detected_x))
        in_y_range_found= False
        if len(in_range) > 0:
            for x in in_range:
                if abs(signal_ch0[anno] - signal_ch0[x]) <= DETECTION_Y_RANGE:
                    t_pos += 1
                    in_y_range_found = True
                    break
            if not in_y_range_found:
                f_neg += 1
        else:
            f_neg += 1

    f_pos = len(detected_x) - t_pos
    print('t_pos:', t_pos, 'f_pos:', f_pos, 'f_neg: ', f_neg)
    return t_pos, f_pos, f_neg

def calculate_stats_for_tests_bitmap(annotated_x, detected_x):
    t_pos = 0
    print('annotated:', len(annotated_x), '/ detected:', len(detected_x))

    if len(annotated_x) > 0:
        bitmap_len = annotated_x[len(annotated_x) - 1] + REF_SAMPLES
        anno_bitmap = np.zeros(bitmap_len)
        det_bitmap = np.zeros(bitmap_len)
        for x in annotated_x:
            for i in range(0, REF_SAMPLES):
                if x - i > 0:
                    anno_bitmap[x - i] = 1
                anno_bitmap[x + i] = 1
        for x in detected_x:
            if x < bitmap_len:
                det_bitmap[x] = 1

        go_to_next_peak = False
        for i in range(0, bitmap_len):
            if anno_bitmap[i]:
                if go_to_next_peak:
                    continue
                if det_bitmap[i]:
                    t_pos += 1
                    go_to_next_peak = True
            else:
                go_to_next_peak = False

        f_neg = len(annotated_x) - t_pos
        f_pos = len(detected_x) - t_pos
    else:
        t_pos = 0
        f_neg = 0
        f_pos = len(detected_x)
    print('t_pos:', t_pos, 'f_pos:', f_pos, 'f_neg: ', f_neg)
    return t_pos, f_pos, f_neg

def get_r_samples(ann):
    return list(filter(lambda x: x[1] in R_SYMBOLS, zip(ann.sample, ann.symbol)))

def get_plot_data():
    record = wfdb.rdrecord(FILEPATH, sampto=5000)
    ann = wfdb.rdann(FILEPATH, 'atr', sampto=5000)
    annotations = get_r_samples(ann)
    anno_peaks_x = list(map(lambda x: x[0], annotations))
    signal_ch0 = list(map(lambda x: x[1], record.p_signal))
    ecg = np.array(signal_ch0)

    # if FILEPATH in ['db/108']:
    #     ecg = np.array(-ecg)

    peaks_r = ff(ecg)
    peaks_y = list(map(lambda x: x[1], peaks_r))
    peaks_x = list(map(lambda x: x[0], peaks_r))



    print(record.__dict__)
    return signal_ch0, peaks_r, peaks_y, peaks_x, anno_peaks_x, ecg


def plot_data():
    signal_ch0, peaks_r, peaks_y, peaks_x, anno_peaks_x, ecg = get_plot_data()
    t_pos, f_pos, f_neg = calculate_stats_for_tests_bitmap(anno_peaks_x, peaks_x)
    # subplot.cla()
    # plt.plot(signal_ch0)
    # plt.plot(peaks_x, peaks_y, 'ro')
    plt.plot(anno_peaks_x, ecg[anno_peaks_x], 'bo')
    # plt.plot(t_pos, ecg[t_pos], 'go')
    # plt.plot(f_pos, ecg[f_pos], 'yo')
    # plt.plot(f_neg, ecg[f_neg], 'rx')
    plt.show()


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

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # fig.subplots_adjust(bottom=0.3)

    plot_data()

    # ref_samples_ax = plt.axes([0.25, 0.15, 0.65, 0.03])
    # ref_slider = Slider(ref_samples_ax, 'ref samples', 1, 100, valinit=REF_SAMPLES, valfmt='%0.0f')

    # # Draw another slider
    # search_saples_ax = plt.axes([0.25, 0.1, 0.65, 0.03])
    # search_slider = Slider(search_saples_ax, 'search samples', 1, 100, valinit=SEARCH_SAMPLES, valfmt='%0.0f')

    # threshold_ax = plt.axes([0.25, 0.20, 0.65, 0.03])
    # threshold_slider = Slider(threshold_ax, 'threshold', 0.01, 1.0, valinit=THRESHOLD)

    # # Draw another slider
    # detections_samples_ax = plt.axes([0.25, 0.05, 0.65, 0.03])
    # detection_slider = Slider(detections_samples_ax, 'detection_samples', 1.0, 100.0, valinit=DETECTION_X_RANGE,
    #                           valfmt='%0.0f')

    # ref_slider.on_changed(sliders_on_changed)
    # search_slider.on_changed(sliders_on_changed)
    # threshold_slider.on_changed(sliders_on_changed)
    # detection_slider.on_changed(sliders_on_changed)
    # plt.show()

