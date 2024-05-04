import pysam as ps
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

depths_tumor_file = "data/results/bwa/tu_depth.txt"
depths_wildtype_file = "data/results/bwa/wt_depth.txt"

depths_tumor = pd.read_csv(depths_tumor_file, sep='\t', names=['chr', 'pos', 'depth'])
depths_wildtype = pd.read_csv(depths_wildtype_file, sep='\t', names=['chr', 'pos', 'depth'])


def makePlot(data):
    plt.plot(data['window'], data['log2'])
    plt.xlabel('window')
    plt.ylabel('log2')
    plt.title('Read depth plot')
    plt.grid(True)
    plt.savefig("basic_plot.png")

def countReadsInFractions(readFile, windowSize):
    start = readFile.iloc[0]['pos'] // windowSize
    end = readFile.iloc[-1]['pos'] // windowSize
    read_depths = np.zeros((end - start + 1))
    for line in readFile.iloc:
        index = (line['pos'] // windowSize) - start
        read_depths[index] += line['depth']
    return read_depths

tu_read_depths = countReadsInFractions(depths_tumor, 1000)
wt_read_depths = countReadsInFractions(depths_wildtype, 1000)

zipped_depths = zip(tu_read_depths, wt_read_depths)

depths = []
for(tu_read, wt_read) in zipped_depths:
    depths.append(np.log2(1 if wt_read == 0 else tu_read/wt_read))

# saving and immediate retrieval is for the sake of division of the data generation and plot construction
depths_to_save = pd.DataFrame(depths)
depths_to_save.to_csv("depths_to_plot")

depths_to_plot = pd.read_csv("depths_to_plot", names=['window', 'log2'])
makePlot(depths_to_plot)






