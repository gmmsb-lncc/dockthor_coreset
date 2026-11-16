"""Library of the qualification-set scripts (scripts/grid_toml_generator.py, scripts/qualification_set.py,
scripts/summarize_results.py, scripts/qualify_coreset.py).

Modules:
    tomlfile   the <pdb>.toml files: TOML literals, the file layout
    dataset    the set on disk: one <pdb>/ per complex (receptors/, ligands/, cofactors/), the root of a TOML
               folder, cloning from git
    molecules  the ligand topologies DockThor reads (.top, MMFFLigand's .cif/.mmcif) and their coordinates
    gridbox    a docking box around a ligand and DockThor's grid-point count, in float32 as DockThor
    docking    running DockThor on the TOMLs
    summary    qualification_summary.csv, success_rate.txt, the two modes of summarize_results.py
    prepare    building a set from CASF-2013 with MMFFLigand and PDBThorBox (qualify_coreset.py)
    affinity   CASF's affinity index converted to index/affinity.csv (prepare)
"""
