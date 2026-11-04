# Novelty check (3 Nov 2026)

Question: has anyone run a Sparse FFT on real HD automotive radar data?
Short answer: I couldn't find anyone who has. Not a claim that nobody has.

| Paper | Real data? | Method | Hardware? | Note |
|---|---|---|---|---|
| Wang, Patel, Petropulu: RFPS-SFT | No, simulation | Projection-slice SFT | No | 3 targets, N0=512, N1=256, 16 antennas. About 3% of the samples, but needs higher SNR than the FFT. No TDM/DDM. |
| Wang, Patel, Petropulu: RSFT | No, simulation | Bucketization with windowing | No | Savings are op counts, not runtime. Fails when K is too big for B. |
| Singh et al., IEEE TIM 2024 | Abstract only so far | Zero-skipping + thresholds, not hashing | Yes | Closest hardware work. Read the full paper before Phase 4. |
| Agarwal et al., FPL | No | Hashing-type Sparse FFT | Yes, FPGA | Million-point transform, not radar. |
| Prediction-Induced Sparsity for Tracked FMCW Radar (arXiv) | Public FMCW dataset | Tracker-based sparsity | Yes, Artix-7 | Different method, but real radar plus FPGA. Read before Phase 4. |

## What exists
- Sparse FFT for automotive radar is shown mostly in simulation (Rutgers group).
- FPGA sparse-FFT hardware exists, but not hashing-type on real radar data.
- Their radar examples use about 3 targets. Real scenes have far more.

## What looks open
- Hashing-type Sparse FFT tested on real RADIal frames.
- Measured sparsity (K, K90, K99) by scene type: city, highway, countryside.
- RADIal uses Doppler-division multiplexing, so every target shows up 12 times
  along Doppler. None of the SFT radar papers I read mention this.
- FPGA cost of the kernel (latency, LUT, DSP, BRAM) under Vitis HLS vs Dynamatic.

## Limits of this check
- Done with web search, not Google Scholar yet. Still to run the three Scholar queries.
- Did not read the Rutgers thesis or the full TIM 2024 paper.
- My expectation: real frames will break the K << B assumption. If so, the
  result is "where it fails", and that is fine to report.

## Wording to use
"We found no prior evaluation of hashing-type Sparse FFT on measured HD
imaging radar data." Never write "first".