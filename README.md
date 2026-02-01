# Compact U(1) lattice gauge theory

C++ code and analysis scripts to run compact U(1) lattice gauge theory. Prepared for the Swansea Univsersity Physics Journal Club.

The C++ code in `code/` uses a Metropolis followed by an overelaxation step to sample the gauge theory.

## Usage

To compile the lattice executable:

```bash
cd code
make
```

To run at a particular value of beta, e.g. 0.92:

```bash
./exec 0.92
```

The results are printed in two columns, 
the first being the average value of the plaquette
and the second the magnetic monopole density.

To run the full analysis, after compiling
the C++ code, from the project home directory:

```bash
cd analysis
./runscript.sh
python analysis.py
```

Generating the data from the lattice using `analysis/runscript.sh` may take some time.
