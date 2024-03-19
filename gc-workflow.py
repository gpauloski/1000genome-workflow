
import csv
import os
import sys
from time import time_ns

from globus_compute_sdk import Executor

from bin.individuals import processing_2
from bin.individuals import processing_chrom_parts
from bin.individuals_merge import merging
from bin.sifting import sifting
from bin.mutation_overlap import run_moverlap
from bin.frequency import run_frequency


def run_workflow(endpoint_id):
    datafile = 'data.csv'
    dataset = '20130502'
    ind_jobs = 1
    execution_site = 'local'
    columns = 'columns.txt'

    populations = []
    for base_file in os.listdir(os.getcwd() + '/data/populations/'):
        f_pop = base_file
        populations.append(f_pop)

    c_nums = []
    individuals_files = []
    sifted_files = []
    sifted_jobs = []
    individuals_merge_jobs = []

    with Executor(endpoint_id=endpoint_id) as gce:
        with open(datafile, 'r') as f:
            for row in csv.reader(f):
                base_file = row[0]
                threshold = int(row[1])
                # To ensure we do not create too many individuals jobs
                ind_jobs = min(ind_jobs, threshold)
                step = threshold // ind_jobs
                rest = threshold % ind_jobs
                if rest != 0:
                    sys.exit("ERROR: for file {}: required individuals jobs {} does not divide the number of rows {}.".format(
                        base_file, ind_jobs, threshold))

                counter = 1

                individuals_jobs = []
                output_files = []

                # Individuals Jobs
                f_individuals = base_file
                base_file_path = os.path.join(os.getcwd(), 'data', dataset, f_individuals)

                c_num = base_file[base_file.find('chr')+3:]
                c_num = c_num[0:c_num.find('.')]
                c_nums.append(c_num)

                individuals_fns = {}
                
                while counter < threshold:
                    stop = counter + step

                    out_name = os.path.join('chr%sn-%s-%s.tar.gz' % (c_num, counter, stop))
                    output_files.append(out_name)
                    chrn_dfs = []
                    
                    
                    f_chrn = out_name

                    f_chrom_parts = gce.submit(
                        processing_chrom_parts,
                        inputfile=base_file_path,
                        columnfile=os.path.join(os.getcwd(), 'data', dataset, columns),
                        c=c_num,
                        counter=counter,
                        stop=stop,
                        total=threshold
                    )
                    
                    chrom_parts = f_chrom_parts.result()
                    
                    #print(f_chrom_parts)
                    print('Done f_chrome_parts')
                    future_chrn_df = []
                    for chrp in chrom_parts[2]:
                        f_chrn_df = gce.submit(
                            processing_2,
                            chrp_element=chrp,
                            data=chrom_parts[0],
                            ndir=chrom_parts[1],
                            c=c_num
                        )
                        future_chrn_df.append(f_chrn_df)
                    
                    for fut in future_chrn_df:
                        chrn_dfs.append(fut.result())

                    print(chrn_dfs)

                    counter = counter + step

                # merge job
                individuals_filename = 'chr%sn.tar.gz' % c_num

                individuals_fns[individuals_filename] = chrn_dfs
                individuals_files.append(individuals_fns)
                
                # Sifting Job
                f_sifting = row[2]
                f_sifting = os.path.join(os.getcwd(), 'data', dataset, 'sifting', f_sifting)

                f_sifted = gce.submit(
                    sifting,
                    inputfile=f_sifting,
                    c=c_num
                )

                sifted_files.append(f_sifted)
                
        sifted_files = [s.result() for s in sifted_files]
        print(f'{sifted_files=}')

        # Analyses jobs
        mutations = []
        frequencies = []
        for i in range(len(individuals_files)):
            for f_pop in populations:

                mutation_res = gce.submit(
                    run_moverlap,
                    input_dir=individuals_files[i],
                    siftfile=sifted_files[i],
                    c=c_nums[i],
                    columns=columns,
                    pop=f_pop
                )
                mutations.append(mutation_res)
                
                frequency_res = gce.submit(
                    run_frequency,
                    input_dir=individuals_files[i],
                    siftfile=sifted_files[i],
                    c=c_nums[i],
                    columns=columns,
                    pop=f_pop
                )
                frequencies.append(frequency_res)

        mutations = [m.result() for m in mutations]
        print(f'{mutations=}')

        frequencies = [freq.result() for freq in frequencies]
        print(f'{frequencies=}')

if __name__ == "__main__":
    start = time_ns()
    run_workflow(endpoint_id=sys.argv[1])
    end = time_ns()
    print(f'Workflow executed with Globus Compute took: {start=}, {end=}, time_ns={end-start}')