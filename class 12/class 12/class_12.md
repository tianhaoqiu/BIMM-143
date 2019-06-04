class 12 - Drug Design
================
Tianhao Qiu
5/9/2019

``` r
library(bio3d)
file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

``` r
hiv <- read.pdb(file.name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

> Q1: What is the name of the two non protein resid values in this structure? What does resid correspond to and how would you get a listing of all reside values in this structure?

``` r
# HOH, MK1
```

``` r
# Create ligand + protein file (two ways)
prot1 <- trim.pdb(hiv, "protein")
lig1  <- trim.pdb(hiv, "ligand")
prot <- atom.select(hiv, "protein", value = T)
lig <-  atom.select(hiv, "ligand", value = T)
write.pdb(prot, file="1hsg_protein.pdb")
write.pdb (lig, file="1hsg_ligand.pdb")
```

Read all.pdbqt
--------------

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi=TRUE)
write.pdb(res, "results.pdb")
```

Nome Mode Analysis(NMA)
-----------------------

``` r
pdb <- read.pdb("1hel")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma( pdb )
```

    ##  Building Hessian...     Done in 0.015 seconds.
    ##  Diagonalizing Hessian...    Done in 0.091 seconds.

``` r
m7 <- mktrj(modes, mode=7, file="mode_7.pdb")

library("bio3d.view")
view(m7, col=vec2color(rmsf(m7)))
```

    ## Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()
