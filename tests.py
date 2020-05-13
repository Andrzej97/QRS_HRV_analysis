from multiprocessing.dummy import freeze_support
import wfdb
import numpy as np
import csv
import main

FILES = ['100', '107', '108', '200', '203', '207', '222', '233']


def get_real_r_peaks(filename):
    anno = wfdb.rdann(filename, 'atr')
    return main.get_r_samples(anno)


def get_found_r_peaks(filename):
    record = wfdb.rdrecord(filename)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    ecg = np.array(signal_ch0)
    return main.find_R_peaks(ecg)


def test_detection(file_numeber):
    filename = 'db/' + file_numeber
    found_r_peaks = get_found_r_peaks(filename)
    found_r_peaks_x = list(map(lambda x: x[0], found_r_peaks))
    real_r_peaks = get_real_r_peaks(filename)
    real_r_peaks_x = list(map(lambda x: x[0], real_r_peaks))
    print('file: ' + filename)
    record = wfdb.rdrecord(filename)
    signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    t_pos, f_pos, f_neg = main.calculate_stats(signal_ch0, real_r_peaks_x, found_r_peaks_x)
    return len(real_r_peaks), len(t_pos), len(f_pos), len(f_neg)


def write_stats_row(writer, file, anno, t_pos, f_pos, f_neg):
    writer.writerow([file, anno, t_pos, "{:.4f}".format(t_pos / anno), f_pos, "{:.4f}".format(f_pos / anno),
                     f_neg, "{:.4f}".format(f_neg / anno)])


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
        write_stats_row(writer, 'TOTAL', total_R_annotations, total_positive, total_false_positive, total_false_negative)