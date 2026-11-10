# Sparse FFT on Real Radar

Is real HD automotive radar data (RADIal) sparse enough for the Sparse
Fourier Transform, and what would a Sparse FFT cost on an FPGA?

Deadline: 29 Nov 2026. Detailed notes: `notes/log.md`.

## Progress (updated 11 Nov 2026)

- [x] Phase 0, setup and novelty check (2-3 Nov)
- [x] Phase 1, baseline radar chain on RADIal (4-10 Nov)
  - [x] Raw ADC to range-Doppler map
  - [x] CA-CFAR, threshold frozen at 4 dB
  - [x] Label check with chance-level control
  - [x] Camera pairing and Figure 1
  - [x] Simple MTI tested (little effect; noted as a limitation)
- [ ] Phase 2, sparsity study (11-14 Nov) - in progress
  - [x] City and highway sequences downloaded
  - [ ] Sanity check on the new sequences
  - [ ] K90 on the unfolded map, with a noise control
  - [ ] Figure 2, one result per scene type
- [ ] Phase 3, Sparse FFT along Doppler (15-20 Nov)
  - [ ] Read Hassanieh et al., SODA 2012
  - [ ] Implement and test
  - [ ] Figure 3, recall against computation
  - [ ] Optional: compare phase shifts
- [ ] Phase 4, FPGA cost (mostly cut; decide 20 Nov)
- [ ] Phase 5, write-up and submission (26-29 Nov)
