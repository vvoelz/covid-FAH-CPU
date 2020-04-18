## Directory contents

Here are some examples of FEP simulations being readied for production on FAH:

```
100_ligands                  -- the 100 test ligands from Tim Dudgeon, prepared OpenMM --> gmx 5.0.4 by Matt
100_ligands_FEPready         -- FEP-ready projects 

100_ligands_noreceptor       -- the 100 test ligands, alone in a small box of solvent
100_ligands_noreceptor_FEPready   -- FEP-ready versions

7_hits
7_hits_FEPready -- b
```

Note: To scp these projects to FAH, you can make a compressed tarball using, e.g:
`tar -zcvf 100_ligands_FEPready.tar.gz 100_ligands_FEPready`

Also:

* `example-from-tim` - the original example of fragment hit and docked compounds from Tim Dudgeon (InformaticsMatters)
* `example_FEP` - initial noodling scripts and testing for nutlin (1MQ) - MDM2 FEP calculations.  Can we get ee to work?
* `hits` -  prep of fragment hits from DiamondMX
* `test_set` - prep (by Dylan?) of the 100 example compounds from Tim



