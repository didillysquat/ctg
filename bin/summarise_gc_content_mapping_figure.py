import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import math
import pandas as pd
import os
from itertools import chain
"""
This will produce a figure that is a plot for every library.
Each plot will plot the various levels of quality control for a bbmap
against the lutea genome, plotting percentage of mapping and gc content
of the mapped and unmapped fastqs.
We have run the mapping after doing the bbmerge for the paired end reads
and after doing nxtrim for the mate pair reads.
As such we have a looooot of files to go through and a lot of plots.
We whould also bear in mind that some of the mappings are paired and some are unpaired (merged (bbmerge) and single end reads (nxtrim))

It's gonna be awesome.
The aim will be to find the sweetspot for the quality control threshold.
We were originally doing this before doing the bbmerge and the nxtrim. The problem was that some of the matepairs reads are not
genuine matepair reads and are single end reads and paired end reads or unknowns.
Becuase we will be putting these different sequence types into the alignment algorithms seperately
it is important to know the size and stdev for each seperately (we will generate this info during the full lutea mapping).
It is also important to know if there is a specific set of minratio settings for the mapping for each of these types
We will find those out here and use them in the full lutea mapping.
We can keep lowering the quality control threshold and more sequences
will continue to map to the genome. However, presumably the number of false positives
will incease due to Symbiodnium DNA mapping to the host. If we are seeing that
mapping percentage continue to increase, but the GC content does not continue to rise
then it is likely that we are not removing any additional host.
This is because the GC content is higher in the zooxs than in the host.

This has worked fairly well. The only problem is that the GC content output doesn't
look very precise, despite having 1000 bins they are still only being populated every 10,
or so so our data has big jumps in it.
"""

class gc_map_figure:
    def __init__(self):
        # First thing to do is to work out how many samples we have and get the names
        # Do this by looking for the *.gchist.txt file
        self.lib_list = self._populate_library_list()
        # We want a plot for each library
        self.rows, self.cols = self._get_dims_for_plot()
        self.f, self.axarr = plt.subplots(self.rows, self.cols, figsize=(10, 10))
        self.axarr = list(chain.from_iterable(self.axarr))

    def _populate_library_list(self):
        files = list(os.walk('.'))[0][2]
        #files = list(os.walk('.'))[0][2]
        return [file_name.replace(".gchist.txt", "") for file_name in files if '.gchist.txt' in file_name]

    def _get_dims_for_plot(self):
        nlib = len(self.lib_list)
        if abs(nlib - int(math.sqrt(nlib))) > abs(nlib - ( int(math.sqrt(nlib) + 1 ))):
            # Then we go with the larger dim and this can just be a square plot
            nrows = int(math.sqrt(nlib)) + 1
            ncols = int(math.sqrt(nlib)) + 1
        else:
            # Then we go with the larger dim for the number of columns, but need
            # to work out the number of rows still
            ncols = int(math.sqrt(nlib))
            nrows = int(nlib / ncols) + 1
        return nrows, ncols

    def plot(self):
        """
        We want to plot 3 things here:
        1 - The percentage of mapping
        2 - The GC content of the mapped
        3 - The GC content of the unmapped
        4 - The GC content of the whole library set
        """
        for lib_name, ax in zip(self.lib_list, self.axarr):
            gc_plotter = GCPlotter(lib_name, ax)
            gc_plotter.plot_unmapped_gc_levels()
            gc_plotter.plot_mapped_gc_levels()
            gc_plotter.plot_mapping_pct()
            gc_plotter.annotate()
            gc_plotter.plot_vline()
        plt.tight_layout()
        plt.savefig("gc_to_mapping_pct.svg")
        plt.savefig("gc_to_mapping_pct.png", dpi=1200)
            
class GCPlotter:
    def __init__(self, lib_name, ax):
        self.ax = ax
        self.ax2 = self.ax.twinx()
        self.lib_name = lib_name
        self.qual_vals = ['1', '0.95', '0.90', '0.85', '0.8', '0.75', '0.7', '0.65', '0.60', '0.55', '0.50', '0.45', '0.40', '0.35', '0.30']
        self.gchist_file_path = os.path.join('.', f'{lib_name}.gchist.txt')
        self.unmapped_gc_path_list = [os.path.join('.', f'{lib_name}.{qual_val}.unmapped.awk_gc.txt') for qual_val in self.qual_vals]
        self.mapped_gc_path_list = [os.path.join('.', f'{lib_name}.{qual_val}.mapped.awk_gc.txt') for qual_val in self.qual_vals]
        self.logged_mapping_path_list = [os.path.join('.', f'{lib_name}.{qual_val}.log.txt') for qual_val in self.qual_vals]
        self.raw_gc_content = self._get_raw_gc_content()
        self.insert_len = None
        self.title = self._generate_title()

    def plot_vline(self):
        self.ax.vlines(x=0.6, ymin=self.ax.get_ylim()[0], ymax=self.ax.get_ylim()[1], color='green')

    def _generate_title(self):
        if 'M_17' in self.lib_name:
            pe_mp = 'pe'
            if 'Lane1' in self.lib_name:
                lane1_lane2 = 'lane1'
            else:
                lane1_lane2 = 'lane2'
            if 'unmerged' in self.lib_name:
                merged_unmerged = 'unmerged'
            else:
                merged_unmerged = 'merged'
            return f'{pe_mp}_{lane1_lane2}_{merged_unmerged}'

        else:
            pe_mp = 'mp'
            if '35-45' in self.lib_name:
                size_range = '35-45'
            elif '5-7' in self.lib_name:
                size_range = '5-7'
            elif '8-11' in self.lib_name:
                size_range = '8-11'
            else:
                raise RuntimeError("expected size range was not found")
            if 'Lane1' in self.lib_name:
                lane1_lane2 = 'lane1'
            else:
                lane1_lane2 = 'lane2'
            if 'se' in self.lib_name:
                lib_type = 'se'
            elif 'pe' in self.lib_name:
                lib_type = 'pe'
            elif 'unknown' in self.lib_name:
                lib_type = 'unknown'
            elif 'mp' in self.lib_name:
                lib_type = 'mp'
            else:
                raise RuntimeError('seq type not recognised')
            return f'{pe_mp}_{size_range}_{lib_type}_{lane1_lane2}'

    def annotate(self):
        if not self.insert_len:
            self.ax.set_title(
                f'{self.title}\nGC content:{self.raw_gc_content:.2f}\nAv. insert len: N/A',
                fontsize='x-small')
        else:
            self.ax.set_title(f'{self.title}\nGC content:{self.raw_gc_content:.2f}\nAv. insert len: {self.insert_len:.2f}', fontsize='x-small')

    def _get_raw_gc_content(self):
        with open(self.gchist_file_path, 'r') as f:
            df = pd.read_csv(f, delimiter='\t', skiprows=range(4))
        return df.iat[df['Count'].idxmax(), 0]


    # def plot_unmapped_gc_levels(self):
    #     max_unmapped_gc_contents = []
    #     for unmapped_path in self.unmapped_gc_path_list:
    #         with open(unmapped_path, 'r') as f:
    #             df = pd.read_csv(f, delimiter='\t')
    #         max_unmapped_gc_contents.append(df.iat[df['scaffolds'].idxmax(), 0]*100)
    #     self.ax.plot([float(_) for _ in self.qual_vals], max_unmapped_gc_contents, '-b')

    def plot_unmapped_gc_levels(self):
        max_unmapped_gc_contents = []
        for unmapped_path in self.unmapped_gc_path_list:
            with open(unmapped_path, 'r') as f:
                val = float(f.read().rstrip())
            max_unmapped_gc_contents.append(val*100)
        self.ax.plot([float(_) for _ in self.qual_vals], max_unmapped_gc_contents, '-b')

    # def plot_mapped_gc_levels(self):
    #     max_mapped_gc_contents = []
    #     for mapped_path in self.mapped_gc_path_list:
    #         with open(mapped_path, 'r') as f:
    #             df = pd.read_csv(f, delimiter='\t')
    #         max_mapped_gc_contents.append(df.iat[df['scaffolds'].idxmax(), 0] * 100)
    #     self.ax.plot([float(_) for _ in self.qual_vals], max_mapped_gc_contents, '-r')

    def plot_mapped_gc_levels(self):
        max_mapped_gc_contents = []
        for mapped_path in self.mapped_gc_path_list:
            with open(mapped_path, 'r') as f:
                val = float(f.read().rstrip())
            max_mapped_gc_contents.append(val*100)
        self.ax.plot([float(_) for _ in self.qual_vals], max_mapped_gc_contents, '-r')

    def plot_mapping_pct(self):
        av_mapped_list = []

        for map_path in self.logged_mapping_path_list:
            with open(map_path, 'r') as f:
                lines = [line.rstrip() for line in f]
            mapped_total = 0
            mapped_val_count = 0
            for line in lines:
                if line.startswith("mapped:"):
                    mapped_val = float(line.split()[1].replace('%',''))
                    mapped_total += mapped_val
                    mapped_val_count += 1
                if line.startswith("insert size avg:"):
                    self.insert_len = float(line.split()[3])
            av_mapped_list.append(mapped_total/mapped_val_count)

        x = [float(_) for _ in self.qual_vals]
        self.ax2.plot(x, av_mapped_list, c='black')

gcmf = gc_map_figure()
gcmf.plot()