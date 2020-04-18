# Note from Yunhui about continuing an expanded ensemble run with `convert-tpr`

This is the line in the project.xml file I mentioned during our meeting:

```
/home/server/server2/projects/Gromacs/p14135/extend-tpr-with-fep-state.py $jobdir/frame$prev-gen.tpr $results/frame$prev-gen.trr $jobdir/frame$gen.tpr $results/PLCpep7_dhdl.xvg 5000
```

The original fils is here (on VAV4): /home/server/server2/projects/Gromacs/p14135/project.xml
This is the key file used:

```
#!/home/server/anaconda2/bin/python
import os, sys
usage = """Usage:
$ extend-tpr-with-fep-state.py [previous-gen.tpr] [previous-gen.trr]
                [output-this-gen.tpr] [previous-gen_dhdl.xvg] [number of ps to extend the tpr]
"""
if len(sys.argv) < 6:
  print usage
  sys.exit(1)

prev_gen_tpr_filename = sys.argv[1]
prev_gen_trr_filename = sys.argv[2]
this_gen_tpr_filename = sys.argv[3]
prev_gen_dhdl_filename = sys.argv[4]
extend_in_ps = int(sys.argv[5])

print 'prev_gen_tpr_filename', prev_gen_tpr_filename
print 'prev_gen_trr_filename', prev_gen_trr_filename
print 'this_gen_tpr_filename', this_gen_tpr_filename
print 'prev_gen_dhdl_filename', prev_gen_dhdl_filename
print 'extend_in_ps', extend_in_ps

# First, let's find the last fep state index sampled
fin = open(prev_gen_dhdl_filename, 'r')
lines = fin.readlines()
fin.close()
fep_state_index = int(lines[-1].split()[1])
print 'The last fep_state_index:', fep_state_index
Testing = False
# Next, let's use this index to build the next tpr file
cmd = '/usr/local/bin/gromacs-5.0.4/bin/gmx convert-tpr -s %s -f %s -o %s -init_fep_state %d -extend %d\n'%(prev_gen_tpr_filename,prev_gen_trr_filename, this_gen_tpr_filename, fep_state_index, extend_in_ps)
print 'cmd: >>>', cmd

if not Testing:
  os.system(cmd)
```


