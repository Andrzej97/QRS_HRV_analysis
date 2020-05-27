from multiprocessing.dummy import freeze_support
import wfdb
import numpy as np
import csv
import main

SAMPLES_PER_SECOND = 360
FILES = ['100', '107', '108', '200', '203', '207', '222', '233']
FILES_FRAGMENTS = {'100': [('00:00:00.000', '00:30:05.600')],
                   '207': [('00:00:00.000', '00:29:08.500')]}


def get_real_r_peaks(filename):
    anno = wfdb.rdann(filename, 'atr')
    return main.get_r_samples(anno)


def get_found_r_peaks(filename):
    record = wfdb.rdrecord(filename)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    ecg = np.array(signal_ch0)
    return main.find_R_peaks(ecg)


def test_detection(file_number):
    filename = 'db/' + file_number
    found_r_peaks = get_found_r_peaks(filename)
    found_r_peaks_x = list(map(lambda x: x[0], found_r_peaks))
    real_r_peaks = get_real_r_peaks(filename)
    real_r_peaks = filter_real_r_peaks(real_r_peaks, file_number)
    real_r_peaks_x = list(map(lambda x: x[0], real_r_peaks))
    print('file: ' + filename)
    record = wfdb.rdrecord(filename)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    t_pos, f_pos, f_neg = main.calculate_stats(signal_ch0, real_r_peaks_x, found_r_peaks_x)
    return len(real_r_peaks), len(t_pos), len(f_pos), len(f_neg)


def write_stats_row(writer, file, anno, t_pos, f_pos, f_neg):
    writer.writerow([file, anno, t_pos, "{:.4f}".format(t_pos / anno), f_pos, "{:.4f}".format(f_pos / anno),
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
    total_R_annotations = 0
    total_positive = 0
    total_false_positive = 0
    total_false_negative = 0
    with open('results.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=',')
        writer.writerow(['File', 'R_annotations', 'Positive', '%ofPositive', 'FalsePositive', '%OfFalsePositive',
                         'FalseNegative', '%OfFalseNegative'])
        for file in FILES:
            anno, t_pos, f_pos, f_neg = test_detection(file)
            write_stats_row(writer, file, anno, t_pos, f_pos, f_neg)
            total_R_annotations += anno
            total_positive += t_pos
            total_false_positive += f_pos
            total_false_negative += f_neg
        write_stats_row(writer, 'TOTAL', total_R_annotations, total_positive, total_false_positive,
                        total_false_negative)
