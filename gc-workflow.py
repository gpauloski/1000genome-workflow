
import csv
import os
import sys
from time import time_ns

from globus_compute_sdk import Executor

from bin.individuals import processing_2
from bin.individuals import processing_chrom_parts
from bin.sifting import sifting
from bin.mutation_overlap import run_moverlap
from bin.frequency import run_frequency

import pandas as pd
import matplotlib.pyplot as plt

def create_gantt(benchmarks, name):
    df = pd.DataFrame(benchmarks)
    colours = {'processing_chrom_parts': 'b', 'processing_2': 'r', 'sifting': 'm', 'mutation_overlap': 'o', 'frequency': 'g'}
    df['colour'] = df['task'].map(colours)
    df['norm_start'] = df['start'] - df['start'].min()
    df.to_pickle(f'{name}.pkl')
    plt.barh(y=df['thread_id'], width=df['duration'], left=df['norm_start'], color=df['color'])
    plt.savefig(f'{name}.pdf')


# class Workflow:
#     datafile: str
#     dataset: str
#     ind_jobs: int
#     columns: str
#     populations: list[str]
#     populations_dir: str

def run_workflow(endpoint_id, debug=False):
    datafile = 'data.csv'
    dataset = '20130502'
    ind_jobs = 1
    execution_site = 'local'
    columns = 'columns.txt'

    populations = []
    for base_file in os.listdir(os.getcwd() + '/data/populations/'):
        f_pop = base_file
        populations.append(f_pop)

    benchmarks = []
    c_nums = []
    individuals_files = []
    sifted_files = []
    output_fns = {}

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

                output_files = []

                # Individuals Jobs
                f_individuals = base_file
                base_file_path = os.path.join(os.getcwd(), 'data', dataset, f_individuals)

                c_num = base_file[base_file.find('chr')+3:]
                c_num = c_num[0:c_num.find('.')]
                c_nums.append(c_num)

                if debug: 
                    c_nums = [1, 2]
                    threshold = 2

                individuals_fns = {}
                
                individuals_filename = 'chr%sn' % c_num

                while counter < threshold:
                    stop = counter + step

                    out_name = os.path.join('chr%sn-%s-%s.tar.gz' % (c_num, counter, stop))
                    output_files.append(out_name)
                    
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
                    chrom_bench = chrom_parts[0]
                    benchmarks.append(chrom_bench)
                    
                    #print(f_chrom_parts)
                    print('Done obtaining chromosome parts')
                    future_chrn_df = []
                    f_chrn_bench = []
                    chrp_elements = chrom_parts[1][2]

                    if debug:
                        chrp_elements = chrp_elements[0:5]
                    
                    futures = []
                    for chrp in chrp_elements:
                        f_chrn_df = gce.submit(
                            processing_2,
                            chrp_element=chrp,
                            data=chrom_parts[1][0],
                            ndir=chrom_parts[1][1],
                            c=c_num,
                        )

                        futures.append(f_chrn_df)

                    for fut in futures:
                        res = fut.result()
                        future_chrn_df.append(res[1])
                        f_chrn_bench.append(res[0])

                    benchmarks.extend(f_chrn_bench)
                    
                    for fut in future_chrn_df:
                        if individuals_filename in output_fns:
                            output_fns[individuals_filename].append(fut)
                        else:
                            output_fns[individuals_filename] = [fut] 
                        
                    counter = counter + step
                
                # Sifting Job
                f_sifting = row[2]
                f_sifting = os.path.join(os.getcwd(), 'data', dataset, 'sifting', f_sifting)

                f_sifted = gce.submit(
                    sifting,
                    inputfile=f_sifting,
                    c=c_num
                )

                sifted_files.append(f_sifted)
                
            # merge task
            individuals_files = {}
            for key, val in output_fns.items():
                for data in val:
                    if key in individuals_files:
                        if data[0] in individuals_files[key]:
                            individuals_files[key][data[0]].append(data[1])
                        else:
                            individuals_files[key][data[0]] = [data[1]]
                    else:
                        individuals_files[key] = { data[0] : [data[1]] }

                        
        sifted_files = [s.result() for s in sifted_files]
        benchmarks.extend([s[0] for s in sifted_files])
        print(f'Sifting completed')

        # Analyses jobs
        mutations = []
        frequencies = []

        if debug:
            individuals_files = list(individuals_files.items())[0:2]
        else:
            individuals_files = individuals_files.items()

        for i, inf in enumerate(individuals_files):
            for f_pop in populations:

                mutation_res = gce.submit(
                    run_moverlap,
                    input_dir=inf[1],
                    siftfile=sifted_files[i][1],
                    c=c_nums[i],
                    columns=columns,
                    pop=f_pop,
                    debug=debug
                )
                mutations.append(mutation_res)
                
                frequency_res = gce.submit(
                    run_frequency,
                    input_dir=inf[1],
                    siftfile=sifted_files[i][1],
                    c=c_nums[i],
                    columns=columns,
                    pop=f_pop,
                    debug=debug
                )
                frequencies.append(frequency_res)

        benchmarks.extend([m.result() for m in mutations])
        benchmarks.extend([freq.result() for freq in frequencies])
        return benchmarks

if __name__ == "__main__":
    start = time_ns()
    benchmarks = run_workflow(endpoint_id=sys.argv[1], debug=bool(int(sys.argv[2])))
    end = time_ns()
    print(f'Workflow executed with Globus Compute took: {start=}, {end=}, time_ns={end-start}')

    df = create_gantt(benchmarks)

    