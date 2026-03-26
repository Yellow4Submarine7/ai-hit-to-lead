# WWP2-C2 Target Provenance

The public target provenance for this repository is `C2_shifted` / WWP2-C2.

## Investigation summary

During migration, a stale document in the private HPC workspace labeled one docking status summary as `6M0J (ACE2-RBD Interface)`.

The execution-side evidence did not support that label:

- the 9M Uni-Dock PBS job used `C2_shifted.pdbqt`
- the GNINA rescoring PBS job used `receptor/C2_shifted.pdbqt`
- the active-learning scripts used `wwp2_c2.pdbqt`
- the checked receptor files matched byte-for-byte across the relevant HPC directories

## Public repository rule

The public repo treats `6M0J` as stale documentation metadata and does not use it in runnable examples, configs, or provenance summaries.
