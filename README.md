# TDALP
Time-dependent-asymmetric-linear-parsimony parametric reconstruction

The time-dependent-asymmetric-linear parsimony is an ancestral state reconstruction method which extends the standard linear parsimony (a.k.a. Wagner parsimony) approach by taking into account both branch lengths and asymmetric evolutionary costs for reconstructing quantitative characters (asymmetric costs rely Asymmetric costs rely on a single parameter; they allow to ``model'' some evolutionary trend).

The time-dependent-asymmetric-linear parsimony infers states which are all taken among the known states, except for some degenerate cases corresponding to special values of the asymmetry parameter (this holds in particular for Wagner parsimony).

The parametric reconstruction of node is the given of all the possible reconstructed states for its node, associated to the corresponding ranges of the asymmetry parameter.

The software 'tdalp' computes the time-dependent-asymmetric-linear-parsimony parametric reconstruction of all the nodes of phylogenetic tree given as input (see below).


Directory "src" contains the C sources of the software.

'example_tree.newick' is a file containing a phylogenetic tree in newick format
'example_states.csv' is a file containing a list of known states in tex/'.csv' format (i.e. lines '<node ident> <state value>')

If the software 'tdalp' is in the same directory as example_tree.newick and example_states.csv, typing in a console with this current directory :

./tdalp -f 1 -s ./example_states.csv ./example_tree.newick example_out.txt

will output two files:

	'example_out.txt': a text file which contains the result of the parametric reconstruction of 'example_tree.newick'/'example_states.csv' in text format. For all nodes <node ident> of 'example_tree.newick', it displays a line:
	'<node ident>'	1.0 >1./3.< 2.0 >2./3.< 3.0 >3./2.< 4.0
	which means that for an asymmetry parameter between 0 and 1./3., the state reconstructed for the node <node ident> is 1.0. It is 2.0 for an asymmetry parameter between 1./3. and 2./3., 3.0 for an asymmetry parameter between 2./3. and 3./2., and 4.0 for an asymmetry parameter greater than 3./2.

	'example_out.pdf': is a pdf document resulting from option '-f 1' and displaying the parametric reconstruction 'example_tree.newick'/'example_states.csv'

A complete description of the options allowed is given below.

---------
| TDALP |
---------

--------------------------
REQUIREMENT

	The software needs the Cairo library.

--------------------------
COMPILING

	Just type
	> make tdalp
	to build the binary.

--------------------------
DESCRIPTION

	'tdalp' computes time-dependent-asymmetric-linear-parsimony parametric reconstruction.


--------------------------
--------------------------
MANUAL
--------------------------
--------------------------


--------------------------

NAME
	tdalp - compute time-dependent-asymmetric-linear-parsimony parametric reconstruction
	
SYNOPSIS
	tdalp [OPTIONS] <inputFile name> <outputFile name>

DESCRIPTION
	Compute the time-dependent-asymmetric-linear-parsimony parametric reconstruction of the tree contained in <inputFile> (it must be in Newick format). It outputs the parametric reconstruction in text format in the file <outputFile> and, when on sets option '-f', in graphic format in the file <outputFile> with the corresponding extension (i.e .pdf, .png etc.)

	Options are
	-s <file name>
		load the file <file name> which must contains the character states in '.csv' format.
			-> a series of lines '<name of the node> <state value>'
	-t <type>
		set the function applied on the branch length of the tree: 
			-t i -> inverse
			-t u -> identity
	-f <number>
		set the graphic format of the output (option is required if on wants a graphic output)
			-f 1 -> pdf
			-f 2 -> postscript
			-f 3 -> png
			-f 4 -> svg
			-f 5 -> pstricks
	-x <file name>
		output the detail of the parametric reconstruction of the root in the file <file name>
	-c <r1> <g1> <b1> <r2> <g2> <b2>
		set the color scale to go from the color (r1,g1,b1) to the color (r2,g2,b2) (rgb codes of colors have their components in [0,1])
	-h
		display help

--------------------------
