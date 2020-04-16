from multiprocessing.dummy import freeze_support
import wfdb
import os

DESTINATION_PATH = os.getcwd() + '/db'
DB = 'mitdb'

def download_all_files():
    wfdb.dl_database(DB, DESTINATION_PATH, records=['100', '107', '108', '200', '203', '207', '222', '233'])

if __name__ == '__main__':
    freeze_support()
    download_all_files()
    record = wfdb.rdrecord('db/100', sampto=3000)
    wfdb.plot_items(signal=record.p_signal, title='MITDB Record 100')


