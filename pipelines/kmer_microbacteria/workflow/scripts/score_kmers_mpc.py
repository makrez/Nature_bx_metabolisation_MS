#!/usr/bin/python3
import pandas as pd
import os
import fire
import time
import numpy as np
from pathlib import Path
import logging
import math
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import as_completed
import queue
import multiprocessing
from functools import partial

from itertools import zip_longest

# def grouper(iterable, n, fillvalue=None):
#     args = [iter(iterable)] * n
#     return zip_longest(*args, fillvalue=fillvalue)
#
# class ProcessPoolExecutorWithQueueSizeLimit(ProcessPoolExecutor):
#     def __init__(self, maxsize=50, *args, **kwargs):
#         super(ProcessPoolExecutorWithQueueSizeLimit, self).__init__(*args, **kwargs)
#         self._work_queue = queue.Queue(maxsize=maxsize)

def calculate_score(row):
    if row['Score'] == 1:
        if row['isPresent'] == 1:
            val = 1
        if row['isPresent'] == 0:
            val = -1
    if row['Score'] == 0:
        #val = 0
        if row['isPresent'] == 1:
             val = -1
        if row['isPresent'] == 0:
             val = 0
    if pd.isnull(row['Score']):
        val = 0
    return val


def score_by_line(csv_file, line_number, phenotype_df):
    count = csv_file.iloc[line_number].to_frame()
    kmer = count.iat[0,0]
    count = count[1:]
    count.index.name = 'Sample'
    count.reset_index(inplace=True)
    count.columns = ['Sample', 'count']
    count = count.replace({'Sample': r'_kmers'}, {'Sample': r''}, regex=True)
    count = count.merge(phenotype_df, how='left', on= 'Sample')
    count['isPresent'] = np.where(count['count'] > 0, 1, 0)
    count['finalScore'] = count.apply(calculate_score, axis=1)
    Score = count['finalScore'].sum()
    data = {'kmer': [kmer], 'Score': [Score]}
    df_to_append = pd.DataFrame(data)
    return df_to_append

def get_scores_dataframe(phenotype_file, csv):
    logging.info(f'Start scoring file {csv}')
    start = time.perf_counter()
    file = pd.read_csv(csv)
    filename = os.path.basename(csv)
    result = pd.DataFrame(columns=['kmer', 'Score'])
    counter = 0
    for line in range(0, len(file)):
        to_append = score_by_line(file, line, phenotype_file)
        logging.info(f'{csv}: {counter}')
        result = result.append(to_append, ignore_index=True)
        counter += 1
    finish = time.perf_counter()
    logging.info(f'Finished scoring {csv}. The operation took {round(finish - start,2)} second(s)')
    #return result
    return (result, csv)

def score_kmer(csv_file_path, phenotype_file,  output_path, log_file):
    """
    """
    Path(output_path).mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.DEBUG, filename=log_file,
                        format='%(asctime)s:%(message)s')
    start_of_script = time.perf_counter()

    file_list = []
    for file in os.listdir(csv_file_path):
        full_path = os.path.join(csv_file_path, file)
        if os.path.isfile(full_path):
            file_list.append(full_path)

    file_list = [file for file in file_list if 'split' in file]
    phenotype = pd.read_csv(phenotype_file)


    POOL_SIZE = 68 #multiprocessing.cpu_count() - 1 # leave 1 for main process


    def compute_chunksize(iterable_size):
        if iterable_size == 0:
            return 0
        chunksize, extra = divmod(iterable_size, POOL_SIZE * 4)
        if extra:
            chunksize += 1
        return chunksize


    with multiprocessing.Pool(POOL_SIZE) as pool:
        chunksize = 15 #compute_chunksize(len(file_list))
        worker = partial(get_scores_dataframe, phenotype)
        # it is assumed that start_processing returns a tuple: (data frame, original filename argument)
        for df, file in pool.imap_unordered(worker, file_list, chunksize):
            csv_filename = os.path.basename(file)
            df.to_csv(os.path.join(f'{output_path}', f'{csv_filename}.csv'), index = False)


if __name__ == '__main__':
    fire.Fire(score_kmer)
