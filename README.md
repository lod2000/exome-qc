# exome-qc

## Install necessary libraries

This requires Python 3, so if you run an operating system that uses Python 2 by default, make sure to explicitly run pip3.

`$ pip install pandas pymongo matplotlib adjustText`

Also install `python-tk` or `python3-tk` if it isn't already.

## Data structure

Before running, make sure your `data` folder looks like this:

```
data/
    ground_truth/
        ground_truth.csv
    panel/
        <small-panel-coverage>.bed
    samples/
        sample1/.../ann.target_transcripts.tab
        sample2/.../ann.target_transcripts.tab
        ...
```

## Run

`$ python src/main.py` (or `python3 src/main.py`)

Alternately, to generate only the combined data file with no analysis or joint callers:

`$ python src/sample_parser.py <ground-truth.xlsx> <small-panel.bed> <tab-archive>`

## Output

A tab-delimited file of data from all samples combined will be generated at `output/combined.tab` This will be used if the script is run again, to avoid generating the same large DataFrame multiple times. Also, a CSV with basic analysis of the variant callers will be generated at `output/analysis.csv`. The script also generates a plot of the different callers. When it pops up, adjust the plot and drag around the labels to the heart's content, then save it as a png.

## Additional callers

Joint callers can be easily added by creating a script and dropping it in the `callers` directory. Make sure the script follows the format below (see modules already provided for more details), and everything should run smoothly.

### Format

Every caller module should look something like this:

```
NAME = 'JOINT_EXAMPLE'

def get_roc(df): # Only necessary for callers with different possible cutoff values
    ...
    return {'TP': true_pos, 'FP': false_pos, 'TN': true_negs, 'FN': false_negs}

def add_caller(df, training): # df is the full sample DataFrame, training is a training subset of df
    ...
    df[NAME] = <list of all variant calls> # Formatted such that positives are marked True and negatives are './.'
```
