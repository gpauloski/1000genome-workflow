
import csv
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import dask

from dask.distributed import Client
from time import perf_counter

from proxystore.connectors.redis import RedisConnector
from proxystore.store import Store
from proxystore.store.future import ProxyFuture

from globus_compute_sdk import Executor
from globus_compute_sdk.sdk.executor import ComputeFuture

from bin.individuals import processing_2
from bin.individuals import processing_chrom_parts
from bin.sifting import sifting
from bin.mutation_overlap import run_moverlap
from bin.frequency import run_frequency
from bin.utils import Bench


def create_gantt(benchmarks, name):
    df = pd.DataFrame(benchmarks)
    print(df)
    colours = {'processing_chrom_parts': 'b', 'processing_2': 'r', 'sifting': 'm', 'mutation_overlap': 'y', 'frequency': 'g'}
    df['colour'] = df['task'].map(colours)
    df['norm_start'] = df['start'] - df['start'].min()
    df.to_pickle(f'{name}.pkl')
    plt.barh(y=df['thread_id'], width=df['duration'], left=df['norm_start'], color=df['colour'])
    plt.savefig(f'{name}.pdf')

class Workflow:
    populations_dir: str = '/data/populations/'
    datafile: str = 'data.csv'
    dataset: str = '20130502'
    ind_jobs: int = 40
    columns: str = 'columns.txt'
    populations: list[str] = []
    benchmarks: list = []
    c_nums: list = []
    individuals_files: list = []
    sifted_files: list =[]
    output_fns: dict = {}

    def __init__(self) -> None:
        for base_file in os.listdir(os.getcwd() + self.populations_dir):
            f_pop = base_file
            self.populations.append(f_pop)

    def run_gc_workflow(self, gce, debug=False):
        with open(self.datafile, 'r') as f:
            for row in csv.reader(f):
                base_file = row[0]
                threshold = int(row[1])
                # To ensure we do not create too many individuals jobs
                self.ind_jobs = min(self.ind_jobs, threshold)
                step = threshold // self.ind_jobs
                rest = threshold % self.ind_jobs
                if rest != 0:
                    sys.exit("ERROR: for file {}: required individuals jobs {} does not divide the number of rows {}.".format(
                        base_file, self.ind_jobs, threshold))

                counter = 1

                output_files = []

                # Individuals Jobs
                f_individuals = base_file
                base_file_path = os.path.join(os.getcwd(), 'data', self.dataset, f_individuals)

                c_num = base_file[base_file.find('chr')+3:]
                c_num = c_num[0:c_num.find('.')]
                self.c_nums.append(c_num)

                if debug: 
                    self.c_nums = [1, 2]
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
                        columnfile=os.path.join(os.getcwd(), 'data', self.dataset, self.columns),
                        c=c_num,
                        counter=counter,
                        stop=stop,
                        total=threshold
                    )
                    
                    chrom_parts = f_chrom_parts.result()
                    chrom_bench = chrom_parts[0]
                    self.benchmarks.append(chrom_bench)
                    
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

                    self.benchmarks.extend(f_chrn_bench)
                    
                    for fut in future_chrn_df:
                        if individuals_filename in self.output_fns:
                            self.output_fns[individuals_filename].append(fut)
                        else:
                            self.output_fns[individuals_filename] = [fut] 
                        
                    counter = counter + step

                # Sifting Job
                f_sifting = row[2]
                f_sifting = os.path.join(os.getcwd(), 'data', self.dataset, 'sifting', f_sifting)

                f_sifted = gce.submit(
                    sifting,
                    inputfile=f_sifting,
                    c=c_num
                )

                self.sifted_files.append(f_sifted)
                
            # merge task
            individuals_files = {}
            for key, val in self.output_fns.items():
                for data in val:
                    if key in individuals_files:
                        if data[0] in individuals_files[key]:
                            individuals_files[key][data[0]].append(data[1])
                        else:
                            individuals_files[key][data[0]] = [data[1]]
                    else:
                        individuals_files[key] = { data[0] : [data[1]] }

                        
        self.sifted_files = [s.result() for s in self.sifted_files]
        self.benchmarks.extend([s[0] for s in self.sifted_files])
        print(f'Sifting completed')

        # Analyses jobs
        mutations = []
        frequencies = []

        if debug:
            individuals_files = list(individuals_files.items())[0:2]
        else:
            individuals_files = individuals_files.items()

        for i, inf in enumerate(individuals_files):
            for f_pop in self.populations:

                mutation_res = gce.submit(
                    run_moverlap,
                    input_dir=inf[1],
                    siftfile=self.sifted_files[i][1],
                    c=self.c_nums[i],
                    columns=self.columns,
                    pop=f_pop,
                    debug=debug
                )
                mutations.append(mutation_res)
                
                frequency_res = gce.submit(
                    run_frequency,
                    input_dir=inf[1],
                    siftfile=self.sifted_files[i][1],
                    c=self.c_nums[i],
                    columns=self.columns,
                    pop=f_pop,
                    debug=debug
                )
                frequencies.append(frequency_res)

        self.benchmarks.extend([m.result() for m in mutations])
        self.benchmarks.extend([freq.result() for freq in frequencies])
        return self.benchmarks


    def run_proxy_workflow(self, gce, store, debug=False):
        with open(self.datafile, 'r') as f:
            for row in csv.reader(f):
                base_file = row[0]
                threshold = int(row[1])
                # To ensure we do not create too many individuals jobs
                self.ind_jobs = min(self.ind_jobs, threshold)
                step = threshold // self.ind_jobs
                rest = threshold % self.ind_jobs
                if rest != 0:
                    sys.exit("ERROR: for file {}: required individuals jobs {} does not divide the number of rows {}.".format(
                        base_file, self.ind_jobs, threshold))

                counter = 1

                output_files = []

                # Individuals Jobs
                f_individuals = base_file
                base_file_path = os.path.join(os.getcwd(), 'data', self.dataset, f_individuals)

                c_num = base_file[base_file.find('chr')+3:]
                c_num = c_num[0:c_num.find('.')]
                self.c_nums.append(c_num)

                if debug: 
                    self.c_nums = [1, 2]
                    threshold = 2

                individuals_fns = {}
                
                individuals_filename = 'chr%sn' % c_num

                while counter < threshold:
                    stop = counter + step

                    out_name = os.path.join('chr%sn-%s-%s.tar.gz' % (c_num, counter, stop))
                    output_files.append(out_name)
                    chrp_data_future: ProxyFuture[tuple] = store.future(polling_interval=0.001)
                    data_future: ProxyFuture[...] = store.future(polling_interval=0.001)

                    chrom_bench = gce.submit(
                        processing_chrom_parts,
                        inputfile=base_file_path,
                        columnfile=os.path.join(os.getcwd(), 'data', self.dataset, self.columns),
                        c=c_num,
                        counter=counter,
                        stop=stop,
                        total=threshold,
                        chrp_data_future=chrp_data_future,
                        data_future=data_future
                    )

                    ndir, chrp_elements = chrp_data_future.result()
                    self.benchmarks.append(chrom_bench.result())
                    
                    print('Done obtaining chromosome parts')
                    future_chrn_df = []
                    f_chrn_bench = []

                    if debug:
                        chrp_elements = chrp_elements[0:5]
                    

                    futures = [store.future(polling_interval=0.001) for i in range(len(chrp_elements))]
                    #pf_chrn_df: ProxyFuture[(pd.DataFrame, str)] = store.future()

                    for i,chrp in enumerate(chrp_elements):
                        f_chrn_bench.append(
                            gce.submit(
                                processing_2,
                                chrp_element=chrp,
                                data=data_future.proxy(),
                                ndir=ndir,
                                c=c_num,
                                pfuture=futures[i],
                            )
                        )

                        name = f'chr{c_num}.{chrp["name"]}'
                        if individuals_filename not in self.output_fns:
                            self.output_fns[individuals_filename] = []
                        self.output_fns[individuals_filename].append(
                            (name, futures[i].proxy())
                        )

                    self.benchmarks.extend(f_chrn_bench)

                    counter = counter + step

                # Sifting Job
                f_sifting = row[2]
                f_sifting = os.path.join(os.getcwd(), 'data', self.dataset, 'sifting', f_sifting)

                pf_sifting = store.future(polling_interval=0.001)
                f_sifted = gce.submit(
                    sifting,
                    inputfile=f_sifting,
                    c=c_num,
                    pfuture=pf_sifting
                )

                self.sifted_files.append(pf_sifting.proxy())
                self.benchmarks.append(f_sifted)

            # merge task
            print('Merging results')
            individuals_files = {}
            for key, val in self.output_fns.items():
                for name, data_proxy in val:
                    data_proxy = data_proxy.result()
                    if key in individuals_files:
                        if name in individuals_files[key]:
                            individuals_files[key][name].append(data_proxy)
                        else:
                            individuals_files[key][name] = [data_proxy]
                    else:
                        individuals_files[key] = { name : [data_proxy] }

        print(f'Sifting completed')

        # Analyses jobs
        mutations = []
        frequencies = []

        if debug:
            individuals_files = list(individuals_files.items())[0:2]
        else:
            individuals_files = individuals_files.items()

        for i, inf in enumerate(individuals_files):
            for f_pop in self.populations:

                mutation_res = gce.submit(
                    run_moverlap,
                    input_dir=inf[1],
                    siftfile=self.sifted_files[i],
                    c=self.c_nums[i],
                    columns=self.columns,
                    pop=f_pop,
                    debug=debug
                )
                # mutation_res = run_moverlap(
                #     input_dir=inf[1],
                #     siftfile=self.sifted_files[i],
                #     c=self.c_nums[i],
                #     columns=self.columns,
                #     pop=f_pop,
                #     debug=debug
                # )
                mutations.append(mutation_res)
                
                frequency_res = gce.submit(
                    run_frequency,
                    input_dir=inf[1],
                    siftfile=self.sifted_files[i],
                    c=self.c_nums[i],
                    columns=self.columns,
                    pop=f_pop,
                    debug=debug
                )
                frequencies.append(frequency_res)

        self.benchmarks = [b.result() if isinstance(b, ComputeFuture) else b for b in self.benchmarks]
        self.benchmarks.extend([m.result() for m in mutations])
        print('Collected mutations')
        self.benchmarks.extend([freq.result() for freq in frequencies])
        print('Collected frequencies')
        return self.benchmarks

    def run_dask_delayed_wf(self, debug=False):
        client = Client()
        with open(self.datafile, 'r') as f:
            for row in csv.reader(f):
                base_file = row[0]
                threshold = int(row[1])
                # To ensure we do not create too many individuals jobs
                self.ind_jobs = min(self.ind_jobs, threshold)
                step = threshold // self.ind_jobs
                rest = threshold % self.ind_jobs
                if rest != 0:
                    sys.exit("ERROR: for file {}: required individuals jobs {} does not divide the number of rows {}.".format(
                        base_file, self.ind_jobs, threshold))

                counter = 1

                output_files = []

                # Individuals Jobs
                f_individuals = base_file
                base_file_path = os.path.join(os.getcwd(), 'data', self.dataset, f_individuals)

                c_num = base_file[base_file.find('chr')+3:]
                c_num = c_num[0:c_num.find('.')]
                self.c_nums.append(c_num)

                if debug: 
                    self.c_nums = [1, 2]
                    threshold = 2

                individuals_filename = 'chr%sn' % c_num

                while counter < threshold:
                    stop = counter + step

                    out_name = os.path.join('chr%sn-%s-%s.tar.gz' % (c_num, counter, stop))
                    output_files.append(out_name)
                    
                    f_chrom_parts = dask.delayed(processing_chrom_parts)(
                        inputfile=base_file_path,
                        columnfile=os.path.join(os.getcwd(), 'data', self.dataset, self.columns),
                        c=c_num,
                        counter=counter,
                        stop=stop,
                        total=threshold,
                        dask=True
                    ).persist()
                    
                    #chrom_parts = f_chrom_parts.compute()
                    chrom_bench = f_chrom_parts[0]
                    self.benchmarks.append(chrom_bench)
                    
                    #print(f_chrom_parts)
                    print('Done obtaining chromosome parts')
                    future_chrn_df = []
                    f_chrn_bench = []
                    chrom_parts = f_chrom_parts[1].compute()
                    chrp_elements = chrom_parts[2]
                    if debug:
                        chrp_elements = chrp_elements[0:5]
                    
                    futures = []
                    for chrp in chrp_elements:
                        f_chrn_df = dask.delayed(processing_2)(
                            chrp_element=chrp,
                            data=chrom_parts[0],
                            ndir=chrom_parts[1],
                            c=c_num,
                            dask=True
                        )

                        futures.append(f_chrn_df)

                    for fut in futures:
                        future_chrn_df.append(fut[1].persist())
                        f_chrn_bench.append(fut[0])

                    self.benchmarks.extend(f_chrn_bench)
                    
                    for fut in future_chrn_df:
                        if individuals_filename in self.output_fns:
                            self.output_fns[individuals_filename].append(fut)
                        else:
                            self.output_fns[individuals_filename] = [fut] 
                        
                    counter = counter + step
                
                # Sifting Job
                f_sifting = row[2]
                f_sifting = os.path.join(os.getcwd(), 'data', self.dataset, 'sifting', f_sifting)

                f_sifted = dask.delayed(sifting)(
                    inputfile=f_sifting,
                    c=c_num,
                    dask=True
                )

                self.sifted_files.append(f_sifted.persist())

            # merge task
            individuals_files = {}
            for key, val in self.output_fns.items():
                for data in val:
                    if key in individuals_files:
                        if data[0] in individuals_files[key]:
                            individuals_files[key][data[0]].append(data[1])
                        else:
                            individuals_files[key][data[0]] = [data[1]]
                    else:
                        individuals_files[key] = { data[0] : [data[1]] }

                        
        #self.sifted_files = [s.result() for s in self.sifted_files]
        self.benchmarks.extend([s[0] for s in self.sifted_files])
        print(f'Sifting completed')

        # Analyses jobs
        mutations = []
        frequencies = []

        if debug:
            individuals_files = list(individuals_files.items())[0:2]
        else:
            individuals_files = individuals_files.items()

        for i, inf in enumerate(individuals_files):
            for f_pop in self.populations:

                mutation_res = dask.delayed(run_moverlap)(
                    input_dir=inf[1],
                    siftfile=self.sifted_files[i][1],
                    c=self.c_nums[i],
                    columns=self.columns,
                    pop=f_pop,
                    debug=debug,
                    
                )
                mutations.append(mutation_res)
                
                frequency_res = dask.delayed(run_frequency)(
                    input_dir=inf[1],
                    siftfile=self.sifted_files[i][1],
                    c=self.c_nums[i],
                    columns=self.columns,
                    pop=f_pop,
                    debug=debug
                )
                frequencies.append(frequency_res)

        # self.benchmarks = client.gather(self.benchmarks)
        self.benchmarks.extend(mutations)
        self.benchmarks.extend(frequencies)
        self.benchmarks = dask.compute(*self.benchmarks)
        client.close()
        return self.benchmarks
        pass

if __name__ == "__main__":
    w = Workflow()
    endpoint_id = sys.argv[1]
    debug = bool(int(sys.argv[2]))
    option = sys.argv[3]

    tic = perf_counter()
    if option == 'dask':
        with dask.config.set(scheduler='processes'):
            benchmarks = w.run_dask_delayed_wf(debug=debug)
            duration = perf_counter() - tic
            print(f'Workflow executed with Dask Delayed took: {duration=}s')
            df = create_gantt(benchmarks=benchmarks, name='dask-trial')
    elif option == 'proxy':
        with Executor(endpoint_id=endpoint_id) as gce:
            with Store('genome-store', RedisConnector('localhost', 6379)) as store:
                benchmarks = w.run_proxy_workflow(gce=gce, store=store, debug=debug)
                duration = perf_counter() - tic
                print(f'Workflow executed with Globus Compute and Proxy Futures took: {duration=}s')
                df = create_gantt(benchmarks=benchmarks, name='proxy-bench-trial')
    else:
        with Executor(endpoint_id=endpoint_id) as gce:
            benchmarks = w.run_gc_workflow(gce=gce, debug=bool(int(sys.argv[2])))
            duration = perf_counter() - tic
            print(f'Workflow executed with Globus Compute took: {duration=}s')
            df = create_gantt(benchmarks=benchmarks, name='gc-trial')

    # end = time_ns()
    # print(f'Workflow executed with Globus Compute took: {start=}, {end=}, time_ns={end-start}')

    # df = create_gantt(benchmarks, 'gc-bench')

    