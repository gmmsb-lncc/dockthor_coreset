"""Library of the qualification-set scripts (scripts/qualification_set.py, scripts/grid_toml_generator.py).

Modules:
    tomlfile   the <pdb>.toml files: TOML literals, the --set option, the file layout
    dataset    the set on disk: systems, config/*.txt boxes, cofactors, cloning from git
    molecules  the ligand topologies DockThor reads (.top, MMFFLigand's .cif/.mmcif) and their coordinates
    gridbox    a docking box around a ligand and DockThor's grid-point count, in float32 as DockThor
    docking    running DockThor on the TOMLs
    summary    qualification_summary.csv, success_rate.txt, --collect and --generate-summary
"""
