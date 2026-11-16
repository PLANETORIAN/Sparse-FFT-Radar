# Sparse FFT on Real Radar

Is real HD automotive radar data (RADIal) sparse enough for the Sparse
Fourier Transform, and what would a Sparse FFT cost on an FPGA?

Deadline: 29 Nov 2026. Detailed notes: `notes/log.md`.

## Progress (updated 16 Nov 2026, evening)

- [x] Phase 0, setup and novelty check (2-3 Nov)
- [x] Phase 1, baseline radar chain on RADIal (4-10 Nov)
  - [x] Raw ADC to range-Doppler map
  - [x] CA-CFAR, threshold frozen at 4 dB
  - [x] Label check with chance-level control
  - [x] Camera pairing and Figure 1
  - [x] Simple MTI tested (little effect; noted as a limitation)
- [x] Phase 2, sparsity study (11-14 Nov)
  - [x] City, countryside and highway sequences checked (Figure 1 per scene)
  - [x] Whole-map K90 with a noise control (ranks scenes backwards; explained by the noise floor share)
  - [x] Peaks per range row over every 5th frame, Figure 2
  - [x] 6 dB check: scene ranking unchanged
  - [x] K in rows with a labelled vehicle: about 36-41 (12 DDM copies), never under 8
- [x] Code moved to src/sfftradar (radial, rd, sparsity); Phase 2 numbers unchanged
- [x] Phase 3, Sparse FFT along Doppler (15-16 Nov)
  - [x] Algorithm: STOC 2012 exactly-sparse (phase method), src/sfftradar/sfft.py
  - [x] Synthetic tests and operation count: cheaper than a dense FFT only for about 1-2 tones per row
  - [x] Multi-channel version with a noise-relative threshold for real rows
  - [x] Real RADIal rows: peaks kept vs cost, against the dense FFT and an energy-based row skip
  - [x] Checks: detection threshold c, realistic row skip
  - [x] Figure 3: the sparse FFT keeps fewer peaks than both alternatives at every cost
- [ ] Open: label frame index vs radar frame (dropped blocks; highway labels offset about 1.8 m)
- [ ] Phase 4, FPGA cost (mostly cut; decide 20 Nov)
- [ ] Phase 5, write-up and submission (report draft 17-19 Nov, final 26-28 Nov)
