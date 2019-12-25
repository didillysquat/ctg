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
It's gonna be awesome.
The aim will be to find the sweetspot for the quality control threshold
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
        return [file_name.replace(".gchist.txt", "") for file_name in files if '.gchist.txt' in file_name]

    def _get_dims_for_plot(self):
        nrows = 0
        ncols = 0
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
        plt.tight_layout()
        foo = 'bar'
            
class GCPlotter:
    def __init__(self, lib_name, ax):
        self.ax = ax
        self.ax2 = self.ax.twinx()
        self.lib_name = lib_name
        self.qual_vals = ['1', '0.95', '0.90', '0.85', '0.8', '0.75', '0.7', '0.65', '0.60', '0.55', '0.50']
        self.gchist_file_path = os.path.join('.', f'{lib_name}.gchist.txt')
        self.unmapped_gc_path_list = [os.path.join('.', f'{lib_name}.{qual_val}.gc_content_hist_unmapped.txt') for qual_val in self.qual_vals]
        self.mapped_gc_path_list = [os.path.join('.', f'{lib_name}.{qual_val}.gc_content_hist_mapped.txt') for qual_val in self.qual_vals]
        self.logged_mapping_path_list = [os.path.join('.', f'{lib_name}.{qual_val}.log.txt') for qual_val in self.qual_vals]
        self.raw_gc_content = self._get_raw_gc_content()
        self.insert_len = None

    def annotate(self):
        self.ax.set_title(f'{self.lib_name}\nGC content:{self.raw_gc_content:.2f}\nAv. insert len: {self.insert_len:.2f}', fontsize='x-small')

    def _get_raw_gc_content(self):
        with open(self.gchist_file_path, 'r') as f:
            df = pd.read_csv(f, delimiter='\t', skiprows=range(4))
        return df.iat[df['Count'].idxmax(), 0]


    def plot_unmapped_gc_levels(self):
        max_unmapped_gc_contents = []
        for unmapped_path in self.unmapped_gc_path_list:
            with open(unmapped_path, 'r') as f:
                df = pd.read_csv(f, delimiter='\t')
            max_unmapped_gc_contents.append(df.iat[df['scaffolds'].idxmax(), 0]*100)
        self.ax.scatter(x = [float(_) for _ in self.qual_vals], y = max_unmapped_gc_contents)

    def plot_mapped_gc_levels(self):
        max_mapped_gc_contents = []
        for mapped_path in self.mapped_gc_path_list:
            with open(mapped_path, 'r') as f:
                df = pd.read_csv(f, delimiter='\t')
            max_mapped_gc_contents.append(df.iat[df['scaffolds'].idxmax(), 0] * 100)
        self.ax.scatter(x=[float(_) for _ in self.qual_vals], y=max_mapped_gc_contents)

    def plot_mapping_pct(self):
        av_mapped_list = []
        for map_path in self.logged_mapping_path_list:
            with open(map_path, 'r') as f:
                lines = [line.rstrip() for line in f]
            mapped_total = 0
            for line in lines:
                if line.startswith("mapped:"):
                    mapped_val = float(line.split()[1].replace('%',''))
                    mapped_total += mapped_val
                if line.startswith("insert size avg:"):
                    self.insert_len = float(line.split()[3])
            av_mapped_list.append(mapped_total/2)
        x = [float(_) for _ in self.qual_vals]
        self.ax2.plot(x, av_mapped_list, c='black')

gcmf = gc_map_figure()
gcmf.plot()