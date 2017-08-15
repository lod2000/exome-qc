# exome-qc

## Install necessary libraries

`sudo pip install pandas pymongo matplotlib adjustText`

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

On Mac or Linux:

`python src/main.py`

On Windows:

`src\main.py`

Alternately, to generate only the combined data file with no analysis or joint callers:

`python src/sample_parser.py <ground-truth.xlsx> <small-panel.bed> <tab-archive>`

## Output

A tab-delimited file of data from all samples combined will be generated at `data/samples/combined.tab` This will be used if the script is run again, to avoid generating the same large DataFrame multiple times. Also, a CSV with basic analysis of the variant callers will be generated at `data/analysis.csv`.
