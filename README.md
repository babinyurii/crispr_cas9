`crispr_count_indels.py` counts exact number of insertions and deletions in fasta file which is the result of reads mapping onto the reference sequence.  The script creates xlsx spreadsheets with the raw count.

`crispr_create_plots.py` uses the output of the `count_indels.py` script and constructs plots and spreadsheets containg the percentage of eash indel and its position according to the reference sequence

The example of output file. Here are all the deletions found in mapping irrespective of their length: 
![bars](example_output/dels_bars.png)

And heatmap which shows both percent of indels and their lengths:
![heatmap](example_output/dels_heatmap.png)
