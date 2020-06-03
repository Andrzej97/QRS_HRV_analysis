from multiprocessing.dummy import freeze_support
import wfdb
import numpy as np
import csv
import main
from timeit import default_timer as timer

SAMPLES_PER_SECOND = 360
FILES = ['100', '107', '108', '200', '203', '207', '222', '233']
FILES_FRAGMENTS = {'100': [('00:00:00.000', '00:30:05.600')],
                   '207': [('00:00:00.000', '00:29:08.500')]}


def get_real_r_peaks(filename):
    anno = wfdb.rdann(filename, 'atr')
    return main.get_r_samples(anno)


def get_found_r_peaks(filename, weight, threshold):
    record = wfdb.rdrecord(filename)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    ecg = np.array(signal_ch0)
    # return main.find_R_peaks(ecg)
    return main.find_R_peaks_weights2(ecg, weight, threshold)


def test_detection(file_number, weight, threshold):
    filename = 'db/' + file_number
    start_time = timer()
    found_r_peaks = get_found_r_peaks(filename, weight, threshold)
    end_time = timer()
    # print('get_found_r_peaks took: ', end_time - start_time)
    found_r_peaks_x = list(map(lambda x: x[0], found_r_peaks))
    start_time = timer()
    real_r_peaks = get_real_r_peaks(filename)
    end_time = timer()
    # print('get_real_r_peaks took: ', end_time - start_time)
    real_r_peaks = filter_real_r_peaks(real_r_peaks, file_number)
    real_r_peaks_x = list(map(lambda x: x[0], real_r_peaks))
    # print('file: ' + filename)
    record = wfdb.rdrecord(filename)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    start_time = timer()
    t_pos, f_pos, f_neg = main.calculate_stats_for_tests_bitmap(signal_ch0, real_r_peaks_x, found_r_peaks_x)
    end_time = timer()
    # print('calculate_stats took: ', end_time - start_time)
    return len(real_r_peaks), t_pos, f_pos, f_neg


def write_stats_row(writer, file, weight, threshold, anno, t_pos, f_pos, f_neg):
    writer.writerow([file, weight, threshold, anno, t_pos, "{:.4f}".format(t_pos / anno), f_pos, "{:.4f}".format(f_pos / anno),
                     f_neg, "{:.4f}".format(f_neg / anno)])


def filter_real_r_peaks(real_r_peaks, filenumber):
    if filenumber not in list(FILES_FRAGMENTS.keys()):
        return real_r_peaks
    filtered_r_peaks = []
    for peak in real_r_peaks:
        if is_peak(peak[0], filenumber):
            filtered_r_peaks.append(peak)
    return filtered_r_peaks


def is_peak(x, filenumber):
    return is_in_some_fragment(x, FILES_FRAGMENTS[filenumber])


def is_in_some_fragment(x, ranges):
    for range in ranges:
        start_sample = convert_time_to_sample(range[0])
        end_sample = convert_time_to_sample(range[1])
        if start_sample <= x <= end_sample:
            return True
    return False


def convert_time_to_sample(time):
    hour = get_hour(time)
    min = get_min(time)
    sec = get_sec(time)
    msec = get_msec(time)
    return SAMPLES_PER_SECOND * (3600 * hour + 60 * min + sec + (msec / 1000))

def get_hour(time):
    return int(time[0:2])

def get_min(time):
    return int(time[3:5])

def get_sec(time):
    return int(time[6:8])

def get_msec(time):
    return int(time[9:12])

if __name__ == '__main__':
    # freeze_support()
    # main.download_all_files()
    start_time = timer()
    total_R_annotations = 0
    total_positive = 0
    total_false_positive = 0
    total_false_negative = 0
    with open('resultsParams.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['File', 'Prev_weight', 'Y_threshold', 'R_annotations', 'Positive', '%ofPositive', 'FalsePositive', '%OfFalsePositive',
                         'FalseNegative', '%OfFalseNegative'])
        weights = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
        y_thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        for weight in weights:
            for threshold in y_thresholds:
                print('weight: ' + str(weight) + ', threshold:' + str(threshold))
                for file in FILES:
                    anno, t_pos, f_pos, f_neg = test_detection(file, weight, threshold)
                    # write_stats_row(writer, file, weight, threshold, anno, t_pos, f_pos, f_neg)
                    total_R_annotations += anno
                    total_positive += t_pos
                    total_false_positive += f_pos
                    total_false_negative += f_neg
                write_stats_row(writer, 'TOTAL', weight, threshold, total_R_annotations, total_positive, total_false_positive,
                                total_false_negative)
    end_time = timer()
    print('Testing time: ' + str(end_time - start_time) + ' seconds')
