# CCB
NRP Comprehensive Combinatorial Biosynthesis (CCB)

NRP-CCB is a tool that allow the generation of large libraries of valid non-ribosomal peptides (NRPs). An NRP is considered valid, if and only if, there is an existing IML for every pair of monomers within that generated peptide,

## Installation
This script uses Python 3.7.x. If you don't have Python, I would recommend downloading it from [Anaconda](https://www.continuum.io/downloads).

Copy or clone this package from Github.

Open the Terminal/Command Line and navigate to where you copied the package:

    $ cd path/to/copied/directory

### Linux and MacOS

Install the dependencies by entering:

    $ pip install -r requirements.txt

## Usage

To run from the command-line, just do:

    $ python CCB_prd.py
Example: Running tool with replacement and with a particular genus
    $ python CCB_prd.py -in path_to/IML_genus_db.csv -o path_to_output/CCB/ -l 3 -r 1 -g Bacillus


## Questions and Comments

Feel free to direct any questions or comments to the Issues page of the repository.

## License

See the [LICENSE](LICENSE.md) file for license rights and limitations (MIT).
