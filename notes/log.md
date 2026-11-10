# Project log: Sparse FFT on Real Radar

## Nov 4-5 (late night): data and the log parser

Got labels_CVPR.csv and one RADIal sequence (RECORD@2020-11-22_12.54.38) into data/. Unzipping nested the folder twice, so the real files sit in data/RECORD@.../RECORD@.../. I keep DATA (inner folder) and SEQ (the name) as two variables in the notebook.

Wrote the parser for the events_log.rec file. A radar line splits into 17 pieces. Sensor is piece 10, sample is 7, offset is 13, datasize is 16. Shorter lines (gps) get skipped. I store radar[chip][sample] = (offset, datasize).

Counts in the log: ch0 463, ch1 459, ch2 462, ch3 462. The .bin files have more blocks than that (469, 465, 468, 468), so the log doesn't list the last 6 of each. Lesson: always read at the offset the log gives, never "block i of the file".

One block is 2,097,152 bytes. That's 1,048,576 int16, which is 512 samples x 256 chirps x 4 Rx as complex numbers (I and Q alternate).

Pairing the chips into frames: take a ch3 sample n, then ch0 is n+1, ch2 is n+2, ch1 is n+3 (this is how Valeo's SyncReader does it). The offsets are different for each chip. For frame 13457, ch3 and ch0 start at 0 but ch2 and ch1 start at 2,097,152.

Result: 462 ch3 blocks, 452 complete frames, 10 incomplete. First frame is 13457, last is 15300. Five sample numbers are just missing from the log (13900, 14139, 14835, 14910, 15005).

Mistake along the way: I first grouped by (n - 13455) // 4 and got rotating chip orders and repeated chips. The recording starts in the middle of a cycle, so the ch3 anchor is the right way.

## Nov 6: cube, range FFT, range-Doppler map

Building the cube for one frame: turn each chip's int16 into complex (x[0::2] + 1j*x[1::2]), reshape to (512, 4, 256) with order="F" because the file goes sample, then Rx, then chirp. Then transpose to (512, 256, 4). Stack the four chips along the last axis in the order ch3, ch0, ch1, ch2. Final cube is (512, 256, 16), complex128, about 33 MB.

Bugs I hit: f.seek() with nothing in it, count given in bytes instead of values, and datasize/2 instead of datasize // 2.

Range FFT plot (Rx 0): first I plotted X and R directly, so I got negative values and a mess of 256 lines. Should have plotted mag_db and avg_db. After fixing: peaks near range bins 85, 160, 198, 217, 240. At 0.2 m per bin that's about 17, 32, 40, 43, 48 m. Noise floor around 63 dB, biggest peak about 96 dB. Small bump at bin 0 is the DC offset, and the floor drops after bin 470. Same peaks show up in the single chirp and the average, so the reshape is fine.

Pipeline, in the same order as Valeo's rpl.py: subtract the DC mean (one per Rx, over samples and chirps, so axis=(0, 1)), Hamming window and range FFT along axis 0, Hamming window and Doppler FFT along axis 1. rd is (512, 256, 16). The map I plotted is only Rx 0. The all-antenna version is the power sum over axis 2, and I'll use that before CFAR. (Check that my saved figure was made with axis=(0, 1) for the mean and not axis=0.)

The DDM comb is real. In the Rx 0 map there are 12 columns at Doppler bins about 66, 82, ... 242, spaced 16, and an empty band at about 2-65 (4 slots of 16). In the row at range bin 198, find_peaks gives 70, 86, 102, 118, 134, 150, 166, 182, 198, 214, 230, 246, every gap exactly 16. Peaks are 112-129 dB, the empty band only has bumps around 93 dB. My guess for Tx0 here is bin 246 (slot just before the empty band), but I haven't checked that, and the Tx1-Tx11 numbering after it is also unchecked.

Taking the 12 largest cells with argsort was wrong, it picked neighbours of the same peak (70/71, 86/87...) and missed the weak ones. Local maxima fixed that.

Things I understood today:
- The FFT acts as if the 512 samples repeat forever. A tone that isn't exactly on a bin leaves a jump at the seam, and that leaks into every bin. Hamming fades the ends so the leakage drops (toy example, tone at bin 10.5, N=64: about -19 dB without a window, about -42 dB with it) but peaks get wider. The window changes how many bins a target covers, so it changes any K I count. Must say which window I used in the write-up.
- One np.fft.fft call with axis= runs lots of 1D FFTs. Range: 256 x 16 = 4,096 FFTs of length 512. Doppler: 512 x 16 = 8,192 of length 256. Using N/2 log2 N, the Doppler stage is about 47% of the FFT multiplications. I had assumed range dominated, which was too strong. Rough count, ignores windows and memory.
- Range migration: the pipeline doesn't correct for a target moving across range bins during a frame. My estimate assumes a 77 GHz carrier and the paper's 0.1 m/s Doppler resolution (about 19.5 ms for 256 chirps). Targets faster than roughly 10 m/s relative smear over a few bins. Just a limitation for the write-up, compensation methods are in the study-next list.

Still open:
- Find Tx0 automatically (4 empty slots in a row, Tx0 is the slot just before) and fold to the 16-bin map, then compare one frame with the SDK.
- Check that the ch3-first grouping is physically right (clean peaks and smooth phase across the 16 channels).
- Plan step 3 says subtract the zero-Doppler mean. That doesn't carry over, because the car moves and the data is DDM, so the static clutter isn't at bin 0 (it's near bin 242 in this frame). Decide the clutter method on Nov 8, before CFAR.
- Haven't looked at the other 15 channels one by one or at any other sequence.

Schedule: started the day about a day behind. The Doppler FFT, map and comb were done this evening, so I'm roughly back on plan. Nov 10 is still the buffer day.

Next (Nov 7): Tx0 finding and the 16-bin fold, check against the SDK on one frame.


## Nov 7-8: fold, Tx0, CA-CFAR

Fold: the all-antenna map is the sum of magnitudes over the 16 channels (np.sum(np.abs(rd), axis=2), the SDK does magnitudes, not powers). Then reshape (512, 16, 16) and sum over the middle axis to get (512, 16). A Doppler bin l is 16a + b, and the 12 transmitters' copies of one target sit 16 bins apart, so they all share column b and the sum stacks them. The fold width has to divide 16 and 16 is the biggest. Gain is about 3x (about 9.5 dB) if I assume the copies are independent. What it costs: two targets whose Doppler differs by a multiple of 16 land in the same cell, and I lose which slot each copy came from, which is why Tx0 has to be found afterwards.

Tx0 (SDK's __find_TX0_position): take the 16 slot values of column b at one range bin, add the first 4 again at the end so it wraps, cumsum, average every 4 in a row, and the smallest window is the empty band. Tx0 bin = 16*j + b. For range bin 198 I expect b = 6, j = 15, Tx0 = 246. I haven't written down the printed result, so check that the cell really prints 6 15 246.

CA-CFAR, copied from Valeo's rpl.py: window (9, 9, 3, 3) so a 19 x 19 mask with a 7 x 7 block removed in the middle (guard cells and the cell under test), which leaves 312 training cells. Work on power (square the map), get the training sums with convolve2d(..., mode="same"), divide by the counts from convolving a map of ones (so the edges are fair), and call it a hit if power / noise > 10^(thr/10). The window does not wrap on the Doppler axis, so columns 0 and 15 have fewer neighbours.

Why a fixed threshold when the noise estimate moves: the false-alarm rate depends on the ratio and the number of training cells, not on the noise level. Valeo uses 2 dB (ratio 1.585). That only works because every folded cell is a sum of 256 magnitudes, so its power barely wobbles (about 6.5% if independent, my simulation gave 6.55%). For a single FFT cell, power is exponential and the textbook formula Pfa = (1 + a/N)^(-N) would give about 21% false alarms at a = 1.585 with N = 312, so that formula does not describe our map. I choose the threshold empirically.

Results on frame 13457 (thr_db 2 / 4 / 6): 593 / 349 / 213 detections. No flat part in the counts, so the threshold alone can't separate targets from clutter. Looking at the plot: columns 1-4 over range bins 30-350 are the static world (the static clutter was near Doppler bin 242, and 242 mod 16 = 2), the moving car is at range bin 198 in columns 5-7, there's a ridge at near range, and range bins 0-3 and 500-511 are edge artefacts to mask. P_red is in linear units (about 3.3e5 to 3.3e8, roughly 110-170 dB). The 0-1 colorbar I got was a plotting quirk, fixed by keeping the imshow handle and using vmin=120. Saved: figures/cfar_folded_thr6.png and the version with the label line.

## Nov 9 (after midnight): labels and the first hit rate

labels_CVPR.csv has 24,179 rows for all sequences. My recording has 556 rows over 405 labelled frames (index 0 to 448). Columns I use: dataset, index (frame number), radar_R_m, radar_A_deg, Annotation. radar_D_mps looks wrong (52 m/s for a car at 18 m), so no Doppler matching until I understand its units.

Range scale: the SDK does RangeBin / 512 * 103, so 103/512 = 0.2012 m per bin. I first wrote 103/150, which is wrong. Frame 0: label at 18.36 m is range bin 91.3. The dashed line lands a few bins from the bright blob at 85-88, so weak evidence. Frame 22: label at 36.7 m is bin 182.6, and the line runs through the brightest return in the whole frame (about 170 dB), in the middle-right columns. Two frames with two ranges both agree, so frames[k] matches label index k as far as I can tell.

Wrote frame_to_pred(n) so any frame gives its folded map. Mistakes: I swapped the sample numbers for ch1 and ch2 (it is ch2 = n+2, ch1 = n+3), and twice put the wrong axis on a sum or FFT (channel sum is axis 2, Doppler FFT is axis 1).

Matching rule: a label is found if some detection is within 5 range bins (about 1 m). Because 200+ detections cover a lot of the range axis, I also measured the chance level: the fraction of empty range bins (more than 10 bins from any label) that are within 5 bins of a detection. I first did this with random fake labels, then switched to the exact coverage after pointing out that a random fake can land on a real vehicle.

Full run, all 556 labels, thr_db 6, gate 5 bins: hit rate h = 0.896, chance c = 0.540, lift (h - c) / (1 - c) = 0.773. Took 2 min 50 s. The first 37 labels gave 0.946 and 0.568, so the small sample was too optimistic. This is a baseline with no clutter removal and no angle processing, so I should not call it a precise detector.

Decisions: no Vivado now (Phase 4 is the stretch and mostly cut, FPGA cost will be an operation count, say so in the report). Clutter removal (plan step 3) is cut for now.

Where I am: about half a day behind, but Nov 10 can be a buffer again if the next day goes well. Left in Phase 1: sweep the threshold at 2, 4, 6, 8, 10 (build all the maps once into a dict so each sweep takes seconds) and freeze one value for Phase 3, get the camera image (the .mjpg index mapping is not verified, give it two hours at most), make Figure 1. SDK cross-check is optional. Rule stays: Figure 1 by Nov 10, otherwise use the SDK's map and move on.

To do: write down my exam dates and how many hours a day I can work, so Phase 3 (Nov 15-20) can be checked.

Next: threshold sweep, then camera, then Figure 1.

## Nov 9 (evening): threshold frozen, camera, Figure 1, clutter tests

Threshold sweep, all 556 labels, gate 5 bins, using cached maps (P_all, one dict entry per labelled frame, about 26 MB, built once in 2.5 min; each threshold then takes seconds). hit rate / chance / lift:
2 dB: 0.978 / 0.890 / 0.804
4 dB: 0.951 / 0.710 / 0.832
6 dB: 0.896 / 0.540 / 0.773
8 dB: 0.865 / 0.391 / 0.779
10 dB: 0.815 / 0.285 / 0.741
Frozen at 4 dB for everything from here on, including the Sparse FFT comparison. The lift beats 6 and 8 dB by only about 0.05 and the frames are neighbours, so it is a modest preference. Picked on the same labels I measure on, fine for a baseline with one knob, but say so in the write-up. Mistakes: I typed dist <= 10 in the chance line instead of dist <= gate (the 10 is only the exclusion zone around real vehicles), and I once left the label loop outside the for k loop, so it ran on one frame.

Camera: the .mjpg is not a video file, it is about 2,700 JPEGs one after another, found through the same events log (offset and datasize, the bytes at the offset start with ff d8 ff e0). I first tried cv2.VideoCapture, which is wrong for this file. Pairing with the radar: Valeo's SDK matches by time. Using their timestamp mode with the -180 ms radar delay gave camera picture 125 for frame 22 and the label box sat beside the car. Matching on timeofissue (their other mode, no delay) gave picture 131 and the box sits on the car. Second labelled frame (ks[200]) also fits with timeofissue. So the camera is paired by timeofissue: camera_for(n) = argmin |cam_tio - ch3_tio[n]|. I do not know why the SDK delay does not apply to the labels, so I describe it as the pairing under which the dataset's boxes fit. Label index is the dataset_index (frame index in the sequence), and the README says the labelled point is the middle of the vehicle's back or front face, which explains the few-bin offsets in range.

Figure 1 is saved as figures/fig1_camera_rd_cfar.png: camera with the green label box on the left, folded range-Doppler map with CFAR detections (4 dB, 357 detections) and the dashed label range line on the right, frame 22, label at 36.7 m (bin 182.6). The vehicle shows up as the brightest streak in the frame in folded columns 9-11. Caption has to say: Doppler is folded (mod 16), detections at the edge columns are mostly static clutter, no clutter removal.

Clutter removal: decided not to build it. Test on 45 labelled frames (every 9th), with a separate script whose CFAR reproduced my counts 213, 349 and 357 exactly (I have not re-run this script myself yet, so these numbers are provisional until I do):
- Clutter columns are only clearly peaked in 8 of 45 frames; the 3 busiest folded columns hold a median 32% of detections (19% if spread evenly).
- Masking those 3 columns would lose 38% of the labelled vehicles (chance 19%) while removing 32% of detections, so it is worse than random. Not doing it.
- The plan's step 3 (subtract the mean over chirps per range bin) changes almost nothing: K 271 to 269, K90 4323 to 4316, lift 0.795 to 0.799. On frame 22: 357 detections against 354. frame_to_pred(n, mti=True) does this, default is off.
- Why it is hard: a static scatterer's Doppler is about v_ego * cos(azimuth), so real removal needs the ego speed, the azimuth (angle FFT and calibration matrix) and Tx0 per cell. Days of work, so it goes in the study-next doc. In Figure 2 I keep two series (DC removal only, plus step 3) and write one sentence saying no ego-motion-compensated suppression was applied, so static returns remain and K/N overstates the moving-target content.

Two warnings for Phase 2:
- The definition of K matters a lot. On the folded map, CFAR detections give K = 271 of 8192 cells (3.3%), but K90 (cells holding 90% of the energy) is 4323 (53%) and K99 is 7059 (86%). First look at 45 frames, on the folded sum-of-magnitudes map, which is not the stage where a Sparse FFT would run. I have to decide where K is measured before quoting a number.
- I have one sequence (looks like countryside). Figure 2 needs city and highway too.

Where I am: Phase 1 is done, one day before the Nov 10 buffer rule. Skipped on purpose: SDK cross-check (optional), clutter removal.

To do: download two more sequences (city and highway; radar .bin files, events log and camera only, labels file already covers all), write down my exam dates and daily hours, commit (git status first, data/ must not show).

Next (Nov 10-11): define K properly, then K/N on this sequence, then the other sequences.

## Nov 10: two more sequences, and what I studied about sparsity

Data: downloaded and unzipped a city sequence (RECORD@2020-11-21_13.46.26) and a highway one (RECORD@2020-11-22_12.25.47), next to the first one (RECORD@2020-11-22_12.54.38, countryside). I picked them from labels_CVPR.csv (vehicles per frame and range as a rough guide to the scene type) and confirmed the scene types by watching the preview videos. The labels file has no scene-type column. Disk is 96% full; the two zips (about 10.5 GB) can be deleted.
Not done yet: Step 0 sanity check on both new sequences (change SEQ, run Figure 1 on the first labelled frame, box on the vehicle and a red dot near the dashed line). Pairing, ch3-first grouping and the range scale were only checked on the first sequence.

What I went through (not yet used in code):
- Sparse means few bins carry most of the energy. K90 = smallest number of bins holding 90% of the energy, found by sorting the bins from largest and adding up. The reference is noise: a pure-noise spectrum needs 149 of 256 bins (58%) for 90%, the 3-tone example needs 2. A map counts as sparse only if K90/N is well below the noise value, so every K90 needs a noise control run through the same pipeline.
- Where to measure: Phase 2 sparsity must be on the unfolded 512 x 256 power map (power = sum over the 16 channels of |.|^2), because the Sparse FFT would replace the 256-point Doppler FFT. The folded 512 x 16 map is the wrong stage. With DDM a single target already has about 12 peaks (one per Tx slot, 16 bins apart) in the unfolded 256 bins, so K cannot go below about 12 per target.
- On the grid and off the grid: bins are the frequencies that fit a whole number of cycles in the window (spacing fs/N); a tone between them leaks into many bins. Hamming window, N = 256 tone: on the grid K90 is 1 without a window and 3 with Hamming; between bins it is 4 without a window (66 bins still within 40 dB of the peak) and 2 with Hamming. Real targets are off the grid, so Hamming gives the lower and more stable K. Also look at bins within X dB of the peak, not only K90.
- Bucketing: keep every (N/B)-th sample and take a B-point FFT; bin f lands in bucket f mod B. With N = 256 and B = 16, spikes at 37, 100 and 201 land in buckets 5, 4 and 9. The bucket value gives the spike size.
- Finding the bin: repeat with the samples shifted by tau (tau < N/B), take the ratio of shifted to unshifted bucket value, and the angle gives f = (N / 2 pi tau) x angle. Then snap to the nearest of the candidates in the bucket (for bucket 5: 5, 21, 37, ...). With tau = 1, N = 256 and B = 16 the candidates are 22.5 degrees apart, so the phase error must stay below 11.25 degrees. A tiny N = 8 example was worked by hand.
- Collisions and operation counts: see the Nov 11 entry below.

Plan for the paper: superseded, see the Nov 11 entry below.

Optional, only if time permits in Phase 3: test several shifts tau (1, 3, and two shifts together). A larger tau gives a finer estimate of an off-grid frequency (frequency error = N / (2 pi tau) x phase error, so 1 degree is 0.71 bins at tau = 1 and 0.24 at tau = 3) but must be used together with tau = 1, never alone (see Nov 11). Start with tau = 1.

Next (Nov 11): Step 0 on both new sequences, then K90 on one labelled frame per scene type, with a noise control.

## Nov 11 (after midnight): buckets, shifts, collisions, which paper

Understood:
- One B-point FFT of a sample group (x0, x16, x32, ... for B = 16) gives all 16 buckets at once. Bucket b is the complete sum of every bin f with f mod 16 = b; a tone at bin 17 is present in every sample, so group 0 sees it even though it does not contain sample 17.
- Group tau = samples starting at tau (x1, x17, ... is group 1). The regular 256-point FFT is 16 group FFTs plus 16 more across the groups (32 FFTs of size 16). The sparse version computes only groups 0 and 1: with one spike in a bucket there are 2 unknowns (size, position), group 0 gives the size, the angle of group 1 / group 0 gives the position, and groups 2-15 are predictable from those two (checked: they match to 1e-15).
- tau rule: use odd tau. Even tau makes candidates in the same bucket share an angle (tau = 2: bins 5 and 133 both at 14.1 degrees). Any odd tau keeps the 16 candidates 22.5 degrees apart, so on-grid tones are found exactly with tau = 1 or 15.
- Off grid, an offset of d bins turns the angle by 1.41 x tau x d degrees; snapping holds only if that is under 11.25 degrees, i.e. d < 8/tau. tau = 1 is always safe (d <= 0.5). tau = 15 fails: 37.6 snaps to bin 21. But once the integer bin is known, a large odd tau reads the fraction finely (37.4 -> 8.4 degrees -> 37.40; 1 degree costs 0.05 bin instead of 0.71). So: tau = 1 for the bin, a large odd tau only as an optional refinement.
- Collisions: bins 37 (1.0) and 53 (0.7) share bucket 5. The bucket reads 1.7, the ratio has angle 61.3 degrees (-> 43.6) and size 0.981. The size test is weak; two shifts that disagree are a better test (tau = 1 -> 37, tau = 3 -> 213).
- Important for RADIal: the 12 DDM copies of a target are exactly 16 bins apart, so with simple every-16th-sample buckets all 12 land in the same bucket (the same thing as the 512 x 16 fold). The papers fix this by permuting (read x at sigma*m + tau, sigma odd, so bin f moves to sigma*f) and using a flat window so each bucket is a block of neighbouring permuted bins (hash = round(sigma*f*B/N)). With sigma = 5 the 12 copies land in 12 different buckets. So Phase 3 must use the permutation, not the simplified bucketing.
- Rough collision chance for one spike: (K-1)/B, so B has to be several times K.

Paper correction: the SODA 2012 paper (Simple and Practical) locates frequencies by voting across rounds, not by phase. The two-sample phase method comes from Nearly Optimal Sparse Fourier Transform (same authors, STOC 2012, arXiv 1201.2501), exactly-sparse algorithm in Section 3. Proposal: implement that one in Phase 3 (simpler, matches what I learned; real data is not exactly sparse, which is the point to test). Decision pending. Reading: STOC 2012 Sections 1-3, skip proofs; SODA 2012 Section 1 only.

Deferred: operation counts (lesson 5) until the end of Phase 3. One point kept: with about 12 bins per DDM target, the saving at N = 256 may be small or none; Phase 2 measures K.

Estimate: Phases 2 and 3 done by about Nov 16 at 5-6 hours a day. Still need to write down exam dates and daily hours.

Next (Nov 11): Step 0 on both new sequences, then frame_to_rd, k_frac and the noise control on one frame per scene.
