
import csv
import os
import sys

from bin.individuals import processing_2
from bin.individuals_merge import merging
from bin.sifting import sifting
from bin.mutation_overlap import run_moverlap
from bin.frequency import run_frequency

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

        while counter < threshold:
            stop = counter + step

            out_name = os.path.join('chr%sn-%s-%s.tar.gz' % (c_num, counter, stop))
            output_files.append(out_name)
            f_chrn = out_name

            f_chrn = processing_2(
                inputfile=base_file_path,
                columnfile=os.path.join(os.getcwd(), 'data', dataset, columns),
                c=c_num,
                counter=counter,
                stop=stop,
                total=threshold
            )

#             # j_individuals = (
#             #     Job('individuals')
#             #         .add_args(f_individuals, c_num, str(counter), str(stop), str(threshold))
#             #         .add_inputs(f_individuals, self.columns)
#             #         .add_outputs(f_chrn, stage_out=False, register_replica=False)
#             # )

#             counter = counter + step

#         # merge job
#         tar_files = [out_name for out_name in output_files]
#         individuals_filename = 'chr%sn.tar.gz' % c_num
#         f_chrn_merged = individuals_filename
#         f_chrn_merged = merging(c=c_num, tar_files=tar_files)
#         # j_individuals_merge = Job('individuals_merge').add_args(c_num)

#         # for out_name in output_files:
#         #     f_chrn = File(out_name)
#         #     j_individuals_merge.add_inputs(f_chrn)
#         #     j_individuals_merge.add_args(f_chrn)

#         # individuals_files.append(f_chrn_merged)
#         # j_individuals_merge.add_outputs(f_chrn_merged, stage_out=False, register_replica=False)
#         # if self.use_decaf or self.use_pmc:
#         #     j_individuals_merge.add_profiles(Namespace.PEGASUS, key="label", value="cluster1")

#         # self.wf.add_jobs(j_individuals_merge)
#         # individuals_merge_jobs.append(j_individuals_merge)
        
#         # Sifting Job
#         f_sifting = row[2]
#         f_sifting_path = os.path.join(os.getcwd(), 'data', dataset, 'sifting', f_sifting)
#         # self.rc.add_replica(site=self.file_site, lfn=f_sifting, pfn=self.src_path +
#         #                     '/data/' + self.dataset + '/sifting/' + f_sifting.lfn)

#         f_sifted = 'sifted.SIFT.chr%s.txt' % c_num
#         sifted_files.append(f_sifted)

#         f_sifted = sifting(inputfile=f_sifting, c=c_num)

#         # j_sifting = (
#         #     Job('sifting')
#         #         .add_inputs(f_sifting)
#         #         .add_outputs(f_sifted, stage_out=False, register_replica=False)
#         #         .add_args(f_sifting, c_num)
#         # )

#         # self.wf.add_jobs(j_sifting)
#         # sifted_jobs.append(j_sifting)

# # Analyses jobs
# for i in range(len(individuals_files)):
#     for f_pop in populations:
#         # Mutation Overlap Job
#         f_mut_out = 'chr%s-%s.tar.gz' % (c_nums[i], f_pop.lfn)

#         run_moverlap(input_dir=individuals_files[i], siftfile=sifted_files[i], c=c_nums[i], columns=columns, pop=f_pop)
#         # j_mutation = (
#         #     Job('mutation_overlap')
#         #         .add_args('-c', c_nums[i], '-pop', f_pop)
#         #         .add_inputs(individuals_files[i], sifted_files[i], c_nums[i], f_pop)
#         #         .add_outputs(f_mut_out, stage_out=True, register_replica=False)
#         # )
#         # Frequency Mutations Overlap Job
#         f_freq_out = 'chr%s-%s-freq.tar.gz' % (c_nums[i], f_pop.lfn)
        
#         run_frequency(input_dir=individuals_files[i], siftfile=sifted_files[i], c=c_nums[i], columns=columns, pop=f_pop)
#         # j_freq = (
#         #     Job('frequency')
#         #         .add_args('-c', c_nums[i], '-pop', f_pop)
#         #         .add_inputs(individuals_files[i], sifted_files[i], f_pop, columns)
#         #         .add_outputs(f_freq_out, stage_out=True, register_replica=False)
#         # )
#         # self.wf.add_jobs(j_mutation, j_freq)