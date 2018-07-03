# CRISPR spacer aligner
Create strong ordering of CRISPR spacers using graphical topology sorting


## Requirements
This script requires:
* `python` 3.6 or above
* `cd-hit` 4.6.1 or above


## Installation
```bash
pip3 install --user git+https://github.com/scwatts/crispr_alignment.git
```


## Usage
Results in `json` format from `CRISPRCasFinder` are used as input:
```bash
mkdir output/
crispr_spacer_alignment.py --input_fps input/*.json --output_prefix output/kp_gc23
```

## TODO
Features
- [x] use CRISPRCasFinder output (json format)
- [ ] collapsed contiguous, harmonious regions in resolve spacer graph

Output files, data:
- [x] spacer content matrix
    - crispr spacer in order of resolved spacer graph
- [ ] node mapping
    - columns: collapsed node, spacer identifier, spacer nucleotide sequence
- [ ] observed spacer order
    - true spacer order per isolate


## Contributors
Stephen Watts, Alex Tokolyi, Kat Holt
