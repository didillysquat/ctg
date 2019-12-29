"""
This script will produce the in_groups.csv and in_libs.csv files that are
used by the PrepareAllPathsInputs.pl script to prepare the input data for the allpathslg algorithm
from fastq files.
You can look at the allpathslg documentation for further details on the format these files should take

pseudo code
For each full path go through a set of conditionals to work out which files is which and fill out the settings
accordingly
The different datasets will be
pe Lane1 merged
pe Lane2 merged
pe Lane1 unmerged
pe Lane2 unmerged
For each of the matepair
mp 3.5-4.5 mp (genuine mate pair reads)
mp 3.5-4.5 pe (look at the stdev for the insert size for this one from mapping) (paired end)
mp 3.5-4.5 se (single end)
NB I think that for the time being we're probably best staying clear of the unknown.
Because these are a mix of types of sequences, I don't see how we can use them.
mp 3.5-4.5 unknown (these are supposed to be a mixture of mate pairs and paired ends according to the nxtrim documentation)
"""

import os
import pandas as pd
import statistics
import sys
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from itertools import chain

class GenerateAllPathsLGDataInputCSVs:
    def __init__(self):
        self.seq_data_files_abs_paths_list = self._populate_library_list()
        self.in_groups_stream = open(os.path.join('.', 'in_groups.csv'), 'w')
        self.in_libs_stream = open(os.path.join('.', 'in_libs.csv'), 'w')
        self._write_headers_to_csv_files()
        self.f, self.axarr = plt.subplots(4,5, figsize=(10, 10))
        self.axarr = list(chain.from_iterable(self.axarr))
        self.ax_count = 0

    def _write_headers_to_csv_files(self):
        self.in_groups_stream.write('group_name,library_name,file_name\n')
        self.in_libs_stream.write('library_name,project_name,organism_name,type,paired,frag_size,frag_stddev,insert_size,insert_stddev,read_orientation,genomic_start,genomic_end\n')

    def _populate_library_list(self):
        """ Return a list of the fastq files. From this list we will
        be able to get the full paths to the *.unmapped.len_info.txt and *.log.txt files
        from which we will get the length and stdev of the sequences for unpaired seq types
        and we can get the insert av length and stdev, respectively
        """
        files = list(os.walk('.'))[0][2]
        return [os.path.abspath(os.path.join('.', file_name)) for file_name in files if 'fastq' in file_name]

    def write_csv_files(self):
        for seq_file in self.seq_data_files_abs_paths_list:
            if 'M_17' in seq_file:
                # Then we are working with pe
                pe_mp = 'pe'
                if 'Lane1' in seq_file:
                    lane1_lane2 = 'lane1'
                else:
                    lane1_lane2 = 'lane2'
                # Then we are working with lane 1
                if '.unmerged.' in seq_file:
                    merged_unmerged = 'unmerged'
                    paired = '1'
                else:
                    merged_unmerged = 'merged'
                    paired = '0'
                    
                # Then we are working with the merged file
                group_name = f'{pe_mp}_{lane1_lane2}_{merged_unmerged}'
                library_name = f'{seq_file.split("/")[-1].replace(".fastq", "")}'
                
                # get the mean size of the fragment and the standard deviation from the len_infofile
                with open(f'{seq_file.replace(".fastq", ".len_info.txt")}', 'r') as f:
                    info = [line.rstrip() for line in f]
                mean_fragment_len = int(float(info[0].split()[2].split('=')[1]))
                stdev_fragment_len = int(float(info[0].split()[3].split('=')[1]))

                self.in_groups_stream.write(f'{group_name},{library_name},{seq_file}\n')
                self.in_libs_stream.write(f'{library_name},ct_genome,C.thermophilum,fragment,{paired},{mean_fragment_len},{stdev_fragment_len},,,inward,,\n')
                if '.unmerged.' in seq_file:
                    # This is purely out of interest. I want to investigate what the insert lengths of the other paired
                    # libraries look like
                    if os.path.isfile(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}'):
                        with open(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}', 'r') as f:
                            lines = [line.rstrip().split() for line in f]
                        # Use this df later for plotting a histogram
                        df = pd.DataFrame(lines[6:], columns=lines[5]).astype(int)
                        # Here we have a pandas hopefully with two columns
                        # Bin in the first and frequency in the second column.

                        fyi_mean_insert_len = int(float(lines[0][1]))
                        fyi_stdev_insert_len = int(float(lines[3][1]))
                        fyi_mode_insert_len = int(float(lines[2][1]))
                        print(
                            f'{pe_mp}_{lane1_lane2}_{merged_unmerged}: mean_insert: {fyi_mean_insert_len}; stdev: {fyi_stdev_insert_len}')
                        self.axarr[self.ax_count].plot(df.iloc[:, 0], df.iloc[:, 1])
                        self.axarr[self.ax_count].set_title(
                            f'{pe_mp}_{lane1_lane2}_{merged_unmerged}\nmean: {fyi_mean_insert_len}; mode: {fyi_mode_insert_len}',
                            fontsize='x-small')
                        self._set_xlim(self.axarr[self.ax_count], fyi_mode_insert_len)
                        self.ax_count += 1
            else:
                pe_mp = 'mp'
                if '35-45' in seq_file:
                    size_range = '35-45'
                elif '5-7' in seq_file:
                    size_range = '5-7'
                elif '8-11' in seq_file:
                    size_range = '8-11'
                else:
                    raise RuntimeError("expected size range was not found")
                if 'Lane1' in seq_file:
                    lane1_lane2 = 'lane1'
                else:
                    lane1_lane2 = 'lane2'
                if '.se.' in seq_file:
                    lib_type = 'se'
                    paired = '0'
                elif '.pe.' in seq_file:
                    lib_type = 'pe'
                    paired = '1'
                elif '.unknown.' in seq_file:
                    # skip this as we won't be using the unknown seq files for the time being
                    # This is purely out of interest. I want to investigate what the insert lengths of the other paired
                    # libraries look like
                    if os.path.isfile(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}'):
                        with open(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}', 'r') as f:
                            lines = [line.rstrip().split() for line in f]
                        # Use this df later for plotting a histogram
                        df = pd.DataFrame(lines[6:], columns=lines[5]).astype(int)
                        # Here we have a pandas hopefully with two columns
                        # Bin in the first and frequency in the second column.

                        fyi_mean_insert_len = int(float(lines[0][1]))
                        fyi_stdev_insert_len = int(float(lines[3][1]))
                        fyi_mode_insert_len = int(float(lines[2][1]))
                        print(
                            f'{pe_mp}_{lane1_lane2}_{size_range}_unknown: mean_insert: {fyi_mean_insert_len}; stdev: {fyi_stdev_insert_len}')
                        self.axarr[self.ax_count].plot(df.iloc[:,0], df.iloc[:,1])
                        self.axarr[self.ax_count].set_title(f'{pe_mp}_{lane1_lane2}_{size_range}_unknown\nmean: {fyi_mean_insert_len}; mode: {fyi_mode_insert_len}', fontsize='x-small')
                        self._set_xlim(self.axarr[self.ax_count], fyi_mode_insert_len)
                        self.ax_count += 1
                    continue
                elif '.mp.' in seq_file:
                    lib_type = 'mp'
                    paired = '1'
                else:
                    raise RuntimeError('seq type not recognised')
                
                # Now set the insert of fragment length and stdev info.
                if 'mp' in seq_file:
                    seq_type = 'jumping'
                    mean_fragment_len = ''
                    stdev_fragment_len = ''
                    # get the mean and stdev insert length from the log file
                    with open(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}', 'r') as f:
                        lines = [line.rstrip().split() for line in f]
                    # Use this df later for plotting a histogram
                    df = pd.DataFrame(lines[6:], columns=lines[5]).astype(int)
                    # Here we have a pandas hopefully with two columns
                    # Bin in the first and frequency in the second column.

                    mean_insert_len = int(float(lines[0][1]))
                    stdev_insert_len = int(float(lines[3][1]))
                    mode_insert_len = int(float(lines[2][1]))
                    print(f'{pe_mp}_{lane1_lane2}_{size_range}_{lib_type}: mean_insert: {mean_insert_len}; stdev: {stdev_insert_len}')
                    self.axarr[self.ax_count].plot(df.iloc[:, 0], df.iloc[:, 1])
                    self.axarr[self.ax_count].set_title(
                        f'{pe_mp}_{lane1_lane2}_{size_range}_{lib_type}\nmean: {mean_insert_len}; mode: {mode_insert_len}',
                        fontsize='x-small')
                    self._set_xlim(self.axarr[self.ax_count], mode_insert_len)
                    self.ax_count += 1

                else:
                    # No need to calculate insert lengths or stdevs.
                    seq_type = 'fragment'
                    # get the mean size of the fragment and the standard deviation from the len_infofile
                    with open(f'{seq_file.replace(".fastq", ".len_info.txt")}', 'r') as f:
                        info = [line.rstrip() for line in f]
                    mean_fragment_len = int(float(info[0].split()[2].split('=')[1]))
                    stdev_fragment_len = int(float(info[0].split()[3].split('=')[1]))
                    mean_insert_len = ''
                    stdev_insert_len = ''

                    # This is purely out of interest. I want to investigate what the insert lengths of the other paired
                    # libraries look like
                    if os.path.isfile(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}'):
                        with open(f'{seq_file.replace(".unmapped.fastq", ".mapped.ihist.txt")}', 'r') as f:
                            lines = [line.rstrip().split() for line in f]
                        # Use this df later for plotting a histogram
                        df = pd.DataFrame(lines[6:], columns=lines[5]).astype(int)
                        # Here we have a pandas hopefully with two columns
                        # Bin in the first and frequency in the second column.

                        fyi_mean_insert_len = int(float(lines[0][1]))
                        fyi_stdev_insert_len = int(float(lines[3][1]))
                        fyi_mode_insert_len = int(float(lines[2][1]))
                        print(
                            f'{pe_mp}_{lane1_lane2}_{size_range}_{lib_type}: mean_insert: {fyi_mean_insert_len}; stdev: {fyi_stdev_insert_len}')
                        self.axarr[self.ax_count].plot(df.iloc[:,0], df.iloc[:,1])
                        self.axarr[self.ax_count].set_title(f'{pe_mp}_{lane1_lane2}_{size_range}_{lib_type}\nmean: {fyi_mean_insert_len}; mode: {fyi_mode_insert_len}', fontsize='x-small')
                        self._set_xlim(self.axarr[self.ax_count], fyi_mode_insert_len)
                        self.ax_count += 1
                    else:
                        foo = 'bar'
                
                group_name = f'{pe_mp}_{lane1_lane2}_{size_range}_{lib_type}'
                library_name = f'{seq_file.split("/")[-1].replace(".fastq", "")}'
                self.in_groups_stream.write(f'{group_name},{library_name},{seq_file}\n')
                self.in_libs_stream.write(f'{library_name},ct_genome,C.thermophilum,{seq_type},{paired},{mean_fragment_len},{stdev_fragment_len},{mean_insert_len},{stdev_insert_len},inward,,\n')
        self.f.suptitle("Insert length distribution infered from P. lutea genome mapping", fontsize="small")
        plt.tight_layout()
        self.f.subplots_adjust(top=0.88)
        self.in_groups_stream.close()
        self.in_libs_stream.close()
        plt.savefig("insert_len_distribution_from_lutea_mapping.svg")
        plt.savefig("insert_len_distribution_from_lutea_mapping.png")

    def _set_xlim(self, ax, mode):
        ax.set_xlim(0, 2 * mode)

gap = GenerateAllPathsLGDataInputCSVs()
gap.write_csv_files()