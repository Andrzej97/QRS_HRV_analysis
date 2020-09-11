import wfdb
import tests
from os import listdir
from os.path import isfile, join

FILES_TO_SKIP = []
FILES = list({f.split('.')[0] for f in listdir('./db') if isfile(join('./db', f))} - set(FILES_TO_SKIP) - {'.'})
MAX_BREAK = 540 #1.5s

def fragmentize_file(file):
    filename = 'db/' + file
    record = wfdb.rdrecord(filename)
    record_signal_ch0 = list(map(lambda x: x[0], record.p_signal))
    anno = wfdb.rdann(filename, 'atr')
    anno_r_peaks = tests.get_r_peaks_from(anno)
    # anno_r_peaks = filter_r_peaks(anno_r_peaks, file)
    anno_r_peaks_x = list(map(lambda x: x[0], anno_r_peaks))
    # print('file: ' + file + str(anno_r_peaks_x))
    print("\'" + file + "\'" + ': [', end='')
    is_first = True
    for i in range(0, len(anno_r_peaks_x) - 1):
        if anno_r_peaks_x[i+1] - anno_r_peaks_x[i] > MAX_BREAK:
            if not is_first:
                print(', ', end='')
            print('(' + str(int(anno_r_peaks_x[i]) + MAX_BREAK - 1) + ', ' + str(anno_r_peaks_x[i+1] - 53) + ')', end='') # 53 - range of detection
            is_first = False
    print('],')

def main():
    for file in FILES:
        fragmentize_file(file)

main()
