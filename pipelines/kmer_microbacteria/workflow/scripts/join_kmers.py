#!/usr/bin/python3
import pandas as pd
import os
import fire
import logging

def join_kmers(file_path, log_file):
    """
    """
    logging.basicConfig(level=logging.DEBUG, filename=log_file,
                        format='%(asctime)s:%(message)s')

    file_list = []
    for file in os.listdir(file_path):
        full_path = os.path.join(file_path, file)
        if os.path.isfile(full_path):
            file_list.append(full_path)

    counter = 0


    for file in file_list:
        sample = os.path.splitext(os.path.basename(file))[0]
        df_new = pd.read_csv(file, sep='\t', names = ['kmer', sample])
        if counter == 0:
            df = df_new
        else:
            df = df.merge(df_new,how='outer', on='kmer')
        counter += 1
        logging.info(f'Processed file {file}')
    print(df.to_csv(index=False, na_rep='NA'), end='')


if __name__ == '__main__':
    fire.Fire(join_kmers)
