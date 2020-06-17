from multiprocessing.dummy import freeze_support
import wfdb
import numpy as np
import main
import json
import functools
import time

SAMPLES_PER_SECOND = 360
# FILES = ['233']
FILES = ['100', '107', '108', '200', '203', '207', '222', '233']
FILES_FRAGMENTS = {'100': [('00:00:00.000', '00:30:05.600')],
                   '207': [('00:00:00.000', '00:29:08.500')]}


def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start_time = time.perf_counter()
        value = func(*args, **kwargs)
        end_time = time.perf_counter()
        run_time = end_time - start_time
        print(f"{func.__name__!r} took {run_time:.4f} secs")
        return value
    return wrapper_timer

def get_r_peaks_from(anno):
    return main.get_r_samples(anno)

# @timer
def get_found_r_peaks(record_signal_ch0, weight, threshold):
    # return main.find_R_peaks(record_signal_ch0)
    return main.find_R_peaks_weights2(record_signal_ch0, weight, threshold)

# @timer
def calculate_stats_for_tests_bitmap(anno_r_peaks_x, found_r_peaks_x):
     return main.calculate_stats_for_tests_bitmap(anno_r_peaks_x, found_r_peaks_x)

def test_detection(file_number, record_signal_ch0, anno_r_peaks_x, weight, threshold):
    found_r_peaks = get_found_r_peaks(record_signal_ch0, weight, threshold)

    found_r_peaks = filter_r_peaks(found_r_peaks, file_number)
    found_r_peaks_x = list(map(lambda x: x[0], found_r_peaks))

    t_pos, f_pos, f_neg = calculate_stats_for_tests_bitmap(anno_r_peaks_x, found_r_peaks_x)
    return t_pos, f_pos, f_neg

def filter_r_peaks(real_r_peaks, filenumber):
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

def create_stats(weight, threshold, t_pos, f_pos, f_neg):
    return {
        'weight': weight,
        'threshold': threshold,
        'stats': {
            'true_positive': t_pos,
            'false_positive': f_pos,
            'false_negative': f_neg
        }
    }

@timer
def test_all(weights, y_thresholds):
    # Map<Filename, List<Result>>
    statistics = { filename : [] for filename in FILES }

    for file in FILES:
        filename = 'db/' + file
        record = wfdb.rdrecord(filename)
        record_signal_ch0 = list(map(lambda x: x[0], record.p_signal))

        anno = wfdb.rdann(filename, 'atr')
        anno_r_peaks = get_r_peaks_from(anno)
        anno_r_peaks = filter_r_peaks(anno_r_peaks, file)
        anno_r_peaks_x = list(map(lambda x: x[0], anno_r_peaks))
        annotation_count = len(anno_r_peaks)

        file_stats = {
            'filename': file,
            'annotation_count': annotation_count,
            'measurements': []
        }

        for weight in weights:
            for threshold in y_thresholds:
                print('weight: ' + str(weight) + ', threshold:' + str(threshold))
                t_pos, f_pos, f_neg = test_detection(file, record_signal_ch0, anno_r_peaks_x, weight, threshold)
                file_stats['measurements'].append(create_stats(weight, threshold, t_pos, f_pos, f_neg))

        statistics[file] = file_stats

    return statistics

if __name__ == '__main__':
    # freeze_support()
    # main.download_all_files()
    
    # weights = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]
    # y_thresholds = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    weights = [0.2]
    y_thresholds = [0.1]

    stats = test_all(weights, y_thresholds)

    f = open('result.json', 'w')
    f.write(json.dumps(stats))
    f.close()
