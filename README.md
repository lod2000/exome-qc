# exome-qc

## Install necessary libraries

`sudo pip install pandas xlrd`

## Data structure

Before running, make sure your `data` folder looks like this:

```
data/
    ground_truth/
        <ground-truth>.xlsx
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

Alternately, to generate only the combined data file with no analysis:

`python src/sample_parser.py <ground-truth.xlsx> <small-panel.bed> <tab-archive>`

## Output

A tab-delimited file of data from all samples combined will be generated at `data/samples/combined.tab`. Also, a CSV with basic analysis of the variant callers will be generated at `data/analysis.csv`.
