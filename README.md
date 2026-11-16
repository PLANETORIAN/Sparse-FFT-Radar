# Sparse FFT on Real Radar

Is real HD automotive radar data (RADIal) sparse enough for the Sparse
Fourier Transform, and what would a Sparse FFT cost on an FPGA?

Deadline: 29 Nov 2026. Detailed notes: `notes/log.md`.

## Progress (updated 16 Nov 2026)

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
- [ ] Phase 3, Sparse FFT along Doppler (15-19 Nov) - in progress
  - [x] Algorithm chosen: STOC 2012 exactly-sparse (phase method)
  - [x] Permutation, flat window, hash to bins, locate (src/sfftradar/sfft.py)
  - [ ] Full loop and operation counter
  - [ ] Synthetic tests: buckets vs K, off-grid tones, noise
  - [ ] Real RADIal rows: detections kept vs dense FFT and skip-empty-rows
  - [ ] Figure 3, recall against computation
- [ ] Open: label frame index vs radar frame (dropped blocks; highway labels offset about 1.8 m)
- [ ] Phase 4, FPGA cost (mostly cut; decide 20 Nov)
- [ ] Phase 5, write-up and submission (26-29 Nov)
