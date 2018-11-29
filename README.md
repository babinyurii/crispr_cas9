The scripts `crispr_count_indels.py` and `crispr_create_plots.py` are used to process the deep sequencing data.

`crispr_count_indels.py` counts exact number of insertions and deletions in fasta file which is the result of reads mapping onto the reference sequence. It creates xlsx spreadsheets with the raw count.

`crispr_create_plots.py` uses the output of the 'count_indels.py' script and makes plots and spreadsheets containg the percentage of eash indel and its position according to the reference sequence

The example of output file. Here are all the deletions found in mapping irrespective of their length: 
![bars](example_output/dels_bars.png)

And heatmap which shows both percent of indels and their lengths:
![heatmap](example_output/dels_heatmap.png)
