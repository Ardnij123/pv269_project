## ModDotPlot

### Where to get ModDotPlot

There are multiple ways to obtain ModDotPlot. Some easy ones are:

- Install it yourself

    - You can get a guide at [https://github.com/marbl/ModDotPlot](github.com/marbl/ModDotPlot)

- Run it as Docker image

```bash
docker pull jindmen/moddotplot
docker run --rm -it --entrypoint bash -v .:<directory you want to use> jindmen/moddotplot:latest
```

### Running ModDotPlot

ModDotPlot has many arguments. First off, there are 2 modes of operation
-- static and interactive. The static mode computes similiarity dotplot
and outputs `.bedpe` file which is effectively an allignment heatmap.
It also allows for exporting this heatmap into various other file formats
for viewing such as `png` or `svg`.
The interactive mode allows for viewing the heatmap in a web browser.

For a static mode, parameters I found out to be very usefull are:

- `--load BEDPE` allows for loading `.bedpe` files outputted in previous runs.
    This allows for instance to generate visual files from a run.

- `--fasta FILE` specifies files with sequences to compare.

- `--output-dir DIR` outputs files into directory other than current.

- `--resolution RESOLUTION` specifies number of 'pixels' to use as resolution
    of heatmap.

- `--window WINDOW` specifies size of a single pixel, e.g. `--window=1000`
    means each pixel instantiates allignment of two 1000bp long subsequences.

- `--identity PERC` specifies minimal identity percentage for a pixel
    to be shown.

- `--compare`, `--compare-only` specifies how to compare sequences. By default,
    moddotplot only compares each sequence with itself. Using `--compare-only`
    it only compares each sequence with other ones, and `--compare` each to each.

- `--no-bedpe`, `--no-plot`, `--no-hist` disables output of some files.
    Usually I have used options `--no-plot` and `--no-hist` and only added
    graphical representation to the bedpe file in a subsequent run using `--load`.
    I do not recommend using `--no-bedpe` option ever.

- `--width W`, `--dpi DPI`, `--deraster` controls size and appearance
    of outputted plot and histogram.

- `--axes-number N` specifies the number of ticks on the axes. This may be useful
    when generating plots for a large number of sequences.

- `--grid` allows for plotting the allignment as a grid. The tool hovewer seems
    not to allow gridding more than 4 sequences into a single graph.

### Running ModDotPlot pt.2

In following text I write some suggestions about running ModDotPlot.

First off, it is often nice to in first run create the BEDPE file
without plotting images, and only then plot the images. Creation of images
chews a lot of RAM.

Time needed for run of ModDotPlot depends on resolution, not window size.
It may be good to first create allignment of the full sequence set in low
resolution, then choose promissing regions/sequences and create a plot in
high resolution of this subset.

When examining the resulting heatmaps, it may be usefull to use a tool such
as `event_counter.py` in this directory. Specifically, this tool allows
to aggregate the heatmap by larger chunks.

Option `--grid` does seem to not allow for allignment of more than 4 sequences.
To overcome this issue, there is a script `merge_split.sh` in this directory
that chains multiple sequences splitting them using strings of `N`'s.

Full run of ModDotPlot might look like following script.

```bash
# Merge sequences
merge_split.sh input.fa -o merged.fa

# Run ModDotPlot in low resolution
moddotplot static -f merged.fa -o w5000 -w 5000 --no-plot --no-hist

# Visualise the bedpe file
event_counter.py w5000/merged.bedpe -w 550000 -s merged.fa.seqs

# Select promissing sequences and plot in high res
seqkit grep input.fa -r -p "chr1_" "chr4_" "chr7_" > group.fa
merge_split.sh group.fa -o group.merged.fa
moddotplot static -f group.merged.fa -o w1000 -w 1000
```
