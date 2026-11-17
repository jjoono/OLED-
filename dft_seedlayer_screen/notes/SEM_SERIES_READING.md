# SEM full-resolution reading: HATCN(2-30)/Ag(3-25) series (2026-08-15)

Nine 1024x768 originals reviewed (70 kX, 5 kV). Key finding first:

## Detector confound
Every Ag3 image = InLens; every Ag5 image = SE2. InLens resolves the
Ag/void labyrinth starkly; SE2 flattens it. The apparent "Ag5 is smooth"
is NOT demonstrated by this set - the detectors differ. ACTION: re-image
the Ag5 films with InLens at identical WD/EHT.

## Per-panel
- Ag3 series (InLens): semicontinuous labyrinth network everywhere -
  bright Ag domains 15-30 nm, dark void channels 5-10 nm wide. Visual
  void fraction 15-25 % (InLens-contrast upper bound). HATCN6/30 add
  bright agglomerated patches (heterogeneity grows with HATCN thickness).
- Ag5 series (SE2): weak granularity + sparse pinholes; smooth-looking
  but not diagnostic (see confound).
- HATCN30_Ag5: dark crack-like linear defects, 50-150 nm long, sparse
  (<1 % area) - consistent with the HATCN 30 nm crack prediction
  (runs/hatcn30_crack.json).
- HATCN2_Ag25: normal grain texture of a thick film.

## Consistency with optics (scripts/72, eta = 0.96)
- Ag3 labyrinth <-> A rising to red (plasmonic tail, 14->17.5 %): MATCH.
- Ag5 absorption 10-13 % despite smooth SE2 look: texture below SE2
  sensitivity, no contradiction.
- Visual void fraction range matches the EMA f_void 0.12-0.24 band.
The "SEM smooth vs optics absorbing" contradiction is RESOLVED as a
detector-choice artefact.

## Open
- Pixel-level quantification (35_sem_texture.py pipeline: autocorrelation
  length, PSD scale, void segmentation) pending original TIFF files
  uploaded as file attachments.
