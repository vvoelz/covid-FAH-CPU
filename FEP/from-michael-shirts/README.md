# Examples and advice about expanded ensemble simulations from Michael Shirts

## Contents

## Email correspondence March 12, 2020

Hi, Vince-
 
 Amidst a HUGE spike in Folding@home participation due to people wanting to help with coronavirus, my lab is going to try doing some massively parallel FEP of drug screening hits.
 
That’s great! 
 
 How confident are you that  expanded-ensemble Wang Landau works "out of the box" in gmx 5.0.4 (the version in the GRO_A7 core)? It think it might be the right tool for the job.  Thoughts?
 
I have a number of example files that were run with 5.0.4, so it should work correctly. I’m attaching files that are adapted for running expanded ensemble in the SAMPLing paper where we used expanded ensemble for host-guest simulation (https://www.biorxiv.org/content/10.1101/795005v1), where we got very comparable answers to other methods.   There might be a couple of vague things – I’m happy to help clarify anything.
 
The one thing that is problematic, and the reason I never got around to publishing, is the convergence.  This showed up in the harder If one is running a bunch of replicates with _independently_ determined weights, it should work, though; you can average weights (or use consensus weights).  We can talk about the best analysis moving forward, but analysis doesn’t necessarily need to be run between generations.
 
What I can’t remember in the current core is how much of the wang-landau information is saved for restarts via checkpointing (though I’m not sure that is how restarts for each generation are running!).  What you would want to do is have the histogram saved, as well as the wang-landau delta and the weights.  I can’t remember what of that is put into the checkpoints.  If not saved, it can slow down convergence.  We may need to do a bit of tweaking.   I worked on the core at some point to try to get this in, so we should be able to re-find the solution relatively quickly.
 
I’m happy to help with any troubleshooting on getting these started, and on the analysis once things are running!
 
Best,
~~~~~~~~~~~~~~~~
Michael Shirts
Associate Professor
michael.shirts@colorado.edu
http://www.colorado.edu/lab/shirtsgroup/
Phone: (303) 735-7860
Office: JSCBB C123
Department of Chemical and Biological Engineering
University of Colorado Boulder
```


