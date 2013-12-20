#!/usr/bin/perl -w


   relax.linuxgccrelease  -database <rosetta database> -in:file:fullatom -out:file:silent_struct_type binary -relax:membrane -membrane:Membed_init -score:weights membrane_highres.wts -in:file:spanfile $basename.span -nstruct 10 -relax_script <dir_ProQ_scripts>/resampling/repack.script -in:file:s *.pdb -out:file:silent ProQM.repacked.silent

Score each of these with ProQM:
    score.linuxgccrelease -database <rosetta database> -in:file:fullatom -ProQ:membrane -score:weights ProQM -in:file:silent ProQM.repacked.silent -out:file:scorefile ProQM.repacked.silent.score
