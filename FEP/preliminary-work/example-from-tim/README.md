# Explanation

Using the 7 fragment *hits* from the xtal screening, 

https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem.html

... the 


From John's forwarded email Mar 10:

```
---------- Forwarded message ---------
From: Tim Dudgeon <tdudgeon.ml@gmail.com>
Date: Tue, Mar 10, 2020 at 4:50 AM
Subject: Re: free energy calcs - Mpro candidates
To: Frank von Delft <frank.von-delft@diamond.ac.uk>, John Chodera <john.chodera@choderalab.org>
Cc: Clyde, Austin Robert <aclyde@anl.gov>


John, Austin,

Attached is some typical data that we are generating. Hopefully it will let you work out what you can do with it.

You'll find:

1. results.sdf.gz - the results of the docking and scoring runs. 100 examples are present. In total there are currently ~20,000 different molecules for each target. More below on the different fields.

2. receptor.mol2 - the protein structure of one of the fragment screening hits (currently we have 7 but more are coming) that was used for the docking. This has had fairly minimal preparation - just removal of all waters and protonation by OpenBabel at pH 7.4

3. The corresponding 7 ligand structures that were used to collect the candidates we tested. 

So #1 contains the poses that generated by docking into #2 and were then scored, and you will want to screen in some way. The candidate molecules were enumerated for charge variants. The original SMILES is found as the Name field of the record (also as the first line) and the CTAB block will be related to that SMILES but may have a different charge state.

The SDF fields that may be of interest are:

Name: the SMILES of the candidate (see above about charges)
SCORE: the rDock docking score
SCORE.norm: the normalised (by HAC) rDock docking score
XChemDeepScore: the predicted score using Jack Scantlebury's deep leaning algorithm (now termed TransFS). 0 means bad, 1 means good.
SuCOS_Max_Score: the best shape/feature overlay score against the 7 fragment screening hit ligands
SuCOS_Max_FeatureMap_Score: ditto but the Feature Map score
SuCOS_Max_Protrude_Score: ditto but the Protrude score
SuCOS_Max_Index: the index of the hit that was the max score (the name of this mol is also supposed to be added as a field but I forgot that)
SuCOS_Cum_Score: the sum of the individual SuCOS scores against all 7 hits.
SuCOS_Cum_Score: ditto for the Feature Map score
SuCOS_Cum_Protrude_Score: ditto but the Protrude score

Let me know if anything doesn't make sense.

Tim 



```



## Email from John Mar 10

```
Greg,

We're posting all the input structures online here:
https://github.com/choderalab/coronavirus

It's a rapidly-evolving work-in-progress, so please open an issue if some information needs to be clariifed!

The DiamondMX fragment screening data is rapidly being posted publicly here:
https://www.diamond.ac.uk/covid-19/for-scientists/Main-protease-structure-and-XChem.html
The DiamondMX docking structures will likely come Friday directly to Diamond to Vince/Matt.

If we want to coordinate a bit better, we are also on an open Slack set up by Alan Aspuru-Guzik to link scientists working on COVID-19. I'll invite you all in case there is some interest in using that as a venue to coordinate SARS-CoV-2 projects (since other interested scientists monitor this too) instead of the FAH slack.

Best,

John

See More from Bowman, Greg

-- 
John D. Chodera
Associate Member, Computational and Systems Biology Program
Memorial Sloan-Kettering Cancer Center
BIH Einstein Visiting Fellow, Charité Universitätsmedizin, Berlin
email: john.chodera@choderalab.org
office: 646.888.3400
fax: 646.888.3105
mobile: 415.867.7384
url: http://www.choderalab.org
```


