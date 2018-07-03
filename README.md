# CRISPR spacer aligner
[![Build Status](https://travis-ci.org/scwatts/crispr_alignment.svg?branch=master)](https://travis-ci.org/scwatts/crispr_alignment)
[![Code Coverage](https://codecov.io/gh/scwatts/crispr_alignment/branch/master/graph/badge.svg)](https://codecov.io/gh/scwatts/crispr_alignment)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Create strong ordering of CRISPR spacers using graphical topology sorting


## Table of contents
* [Requirements](#requirements)
* [Installation](#installation)
* [Usage](#usage)
* [Todo](#todo)
* [Contributors](#contributors)
* [License](#license)


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
- [ ] allow optional use of CD-HIT
- [ ] provide option to search in RC space

Output files, data:
- [x] spacer content matrix
    - crispr spacer in order of resolved spacer graph
- [ ] node mapping
    - columns: collapsed node, spacer identifier, spacer nucleotide sequence
- [ ] observed spacer order
    - true spacer order per isolate


## Contributors
Stephen Watts, Alex Tokolyi, Kat Holt


## License
[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html)
