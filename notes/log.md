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

## Nov 11 (evening): new sequences checked, real highway, Phase 2 started

Data:
- The "highway" sequence RECORD@2020-11-22_12.25.47 turned out to be a two-lane rural road (Figure 1 at label frame 51). Kept only as a spare rural sequence.
- Downloaded RECORD@2020-11-22_12.38.47 instead. Confirmed as a divided highway from preview frames at 10-90% of the recording and from the GPS file: speed 107-117 km/h (median 111 km/h, about 31 m/s). 296 labelled frames, one vehicle each, all beyond 40 m.
- Scenes now: city = 2020-11-21_13.46.26, countryside = 2020-11-22_12.54.38, highway = 2020-11-22_12.38.47. Zips deleted, about 64 GB free.

Step 0 (notebooks/02_check_sequences.ipynb, Phase 1 cells only, SEQ and TAG at the top, figures saved as fig1_<TAG>_k<k>.png):
- City, label frame 1: all 5 boxes on vehicles; detections near the 6.8 m, 14.0 m and 32-34 m labels. Pass.
- Highway, label frame 0: box on the lead car; static clutter forms the curved band expected at speed. The nearest detection to the 70.0 m label is about 7 range bins away (read off the plot), so the 5-bin gate would count it as missed. Possible gate issue at long range; check later.
- Countryside, label frame 0: same as the original Figure 1.
- Not done: the same check at k = ks[100] on city and highway.

Open issue, label frame index: the code maps label index k to frames[k]. Every sequence has dropped radar blocks (10 countryside, 15 city, 7 rural), so after the first drop frames[k] and the k-th ch3 sample differ. Quick test on labelled frames where they differ: hit rate 0.885 vs 0.923 (countryside), 0.833 vs 0.897 (city), 0.682 vs 0.545 (rural); about 0.98 where they agree. Neither mapping is fully right, so the 0.951 hit rate may be slightly low. Needs a timestamp-based mapping before quoting hit rates again. Does not affect Phase 2 (K90 uses no labels). label numSample is not the radar sample number.

Phase 2 started (notebooks/02_sparsity.ipynb, cells 1-7 copied from 02_check_sequences):
- frame_to_rd(n): same as frame_to_pred up to the Doppler FFT, then power summed over the 16 channels, sum |X|^2 over axis 2, with no fold -> (512, 256). Differences from frame_to_pred: power instead of magnitude (K90 counts energy), and no fold (the fold adds a target's 12 DDM copies into one cell, which would make the data look about 12 times sparser than it is for a Doppler sparse FFT).
- Uses frames[0], not labelled frames; K90 needs no labels.

Next: test frame_to_rd (shape 512 x 256), then k_frac, the noise control, and K90 on one frame per scene.

## Nov 12-13: Phase 2 sparsity results, Figure 2

Notebook (notebooks/02_sparsity.ipynb):
- Cells 0-2 replaced by SEQS (tag -> sequence) and load_seq(tag), which sets the globals SEQ, TAG, DATA, radar and frames. read_block and frame_to_rd read these globals, so every later cell works on whichever scene was loaded last. Complete frames: city 576, countryside 452, highway 388 (old cell comments saying highway 483 were from the wrong download).
- k_frac(P, frac): sort all bin powers, running share of total energy, searchsorted for the first point reaching frac, +1 for a count. Test on [50, 30, 10, 6, 4] gives 3.
- Noise control: complex Gaussian (512, 256, 16) through the same windows, FFTs and channel power sum.

Whole-map results, frame 0 of each scene (share of the 131072 bins):

| Scene | K50 | K80 | K90 | K99 | Floor share | Bins > 6 dB | Bins > 10 dB |
|---|---|---|---|---|---|---|---|
| City | 0.12% | 1.65% | 5.29% | 59.8% | 0.055 | 19611 | 11723 (8.9%) |
| Countryside | 0.06% | 1.26% | 13.3% | 75.1% | 0.128 | 8167 | 4369 (3.3%) |
| Highway | 0.34% | 19.4% | 46.1% | 84.8% | 0.291 | 8065 | 3238 (2.5%) |
| Noise | | | 84.6% | | | 0 | 0 |

Floor share = median(P) x number of bins / total energy, a rough estimate of how much of the energy is noise floor (one median for the whole map, although the floor changes with range).

Findings:
- Noise K90 is 84.6%, not 60-70%. One channel of complex noise has exponential bin power (K90 about 59%); summing 16 independent channels makes the bins nearly equal (std/mean = 1/sqrt(16) = 0.25), so K90 approaches the flat-spectrum value of 90%. The windows do not change it.
- K90 ranks the scenes backwards (city sparsest, highway least sparse). K90 measures how the energy is split between peaks and floor, not how many peaks there are. Each curve stays small until the energy level reaches about 1 - floor share, then jumps into the flat floor: city near 94.5%, countryside near 87%, highway near 71%. The highway's 46% comes from its large floor share, not from many targets.
- Counting bins more than 10 dB above the median ranks them as expected: city > countryside > highway. Noise gives 0, so the counted bins are real returns.
- The counts are far above the number of labelled vehicles (highway: one vehicle, 3238 bins). Each target covers several range bins, Hamming leakage widens it, DDM makes 12 copies along Doppler, and static clutter (guardrails, road, signs) appears at Doppler values set by its angle because the radar is moving. K therefore includes clutter.
- So K90 should not be the headline number. Report peaks above the floor.

Per-row results (what the sparse FFT replaces: one 256-point Doppler FFT per range row). K = bins more than 10 dB above that row's own median (per-row median handles the range-dependent floor):

Frame 0: city median 8, p90 44, max 95, rows with K <= 8: 50%; countryside median 0, p90 38, max 79, 69%; highway median 0, p90 16, max 70, 85%.

Every 5th frame, all rows pooled (saved in results/k_rows.npz):

| Scene | Frames | Rows | Median K | p90 | Max | Rows with K = 0 | Rows with K <= 8 |
|---|---|---|---|---|---|---|---|
| City | 116 | 59392 | 2 | 36 | 108 | about 42% | 63.2% |
| Countryside | 91 | 46592 | 0 | 18 | 85 | about 60% | 81.9% |
| Highway | 78 | 39936 | 0 | 11 | 96 | about 73% | 88.5% |

- Frame 0 was busier than typical (city median 8 -> 2, countryside p90 38 -> 18), so one frame was not enough.
- The distribution is very skewed: most rows empty or nearly empty, a minority crowded (up to 108 of 256 bins). The realistic design is a hybrid: sparse FFT on sparse rows, dense FFT on crowded rows.
- DDM catch: a row that contains a real target has about 12 DDM copies, so K >= 12 or so. Rows with K <= 8 are mostly rows without a target. A simple test that skips empty rows already gets much of the saving (42-73% of rows have K = 0). Phase 3 must compare the sparse FFT with "skip empty rows" as well as with the dense FFT, and look at K in rows that contain targets.
- Idea for Phase 3 (not committed): the DDM copies are a known 16 bins apart, so a sparse FFT could recover a target from its copies together instead of treating them as 12 separate peaks.
- K = 8 break-even is still the rough, unmeasured figure; real operation counts come in Phase 3.

Figure 2 (figures/fig2_row_k_cdf.png): share of range rows with at most K peaks, per scene, step plot, dashed line at K = 8 marked illustrative. Check: at K = 8 the curves read 0.63, 0.82, 0.89, matching the table. Caption should say noise gives K = 0 in every row, and that K includes clutter and the 12 DDM copies.

Not done / optional: K in rows that contain labelled targets (needs the timestamp-based label-to-frame mapping); repeat at 6 dB to show the ranking does not depend on the threshold.

Phase 2 status: done except the optional items. About a day ahead of plan.

Next: Phase 3. Confirm the paper (STOC 2012 exactly-sparse algorithm, Section 3), read Sections 1-3 without proofs, then implement the sparse FFT on one range row (N = 256) with the permutation and flat window.

## Nov 13 (evening): the optional Phase 2 items, label offset

6 dB check (every 5th frame, same loop as Figure 2, not saved):

| Scene | Median K | p90 | Max | Rows with K = 0 | Rows with K <= 8 |
|---|---|---|---|---|---|
| City | 16 | 54 | 120 | 22.3% | 39.9% |
| Countryside | 6 | 39 | 106 | 31.2% | 56.6% |
| Highway | 1 | 29 | 108 | 46.3% | 71.4% |

- The ranking (highway sparsest, then countryside, then city) is the same at 6 dB and 10 dB.
- The absolute numbers move a lot: from 10 to 6 dB, rows with K <= 8 drop by 18-23 points. So always quote K with its threshold.
- The threshold should not be a free choice. CFAR detects at 4 dB above the local noise, so judging the sparse FFT by "does CFAR on its output keep the same detections as on the dense FFT" is cleaner than any fixed K threshold. Plan for Phase 3.

K in rows that contain a labelled car:
- Only labelled frames before the first dropped radar block are used (frames[k] == ch3[k]), so label k surely points at radar frame k.
- First run (rows r-1 to r+1 around the label row, r = radar_R_m / 0.201): city 5 cars, median K 38; countryside 91 cars, median 36; highway 32 cars, median 0. The highway 0 was a measurement problem: printing K for rows r-10 to r+10 showed the car's busy block 6-11 rows beyond the label row.
- Offset check: for each car, the busiest row within +-20 rows of the label row, offset = busiest row - label row. Countryside cars at 40-44 m: offset 0. Highway cars at 61-70 m: offset about +9 rows (about 1.8 m, radar sees the car further away). Many scattered points elsewhere are argmax picking clutter.
- Result, keeping cars within 5 rows of their scene's median offset: countryside 55 cars, median K 41, max 69; highway 32 cars, median 36, max 48; city 1 car (too few). K <= 8 in 0% of car rows in both scenes.
- Why about 36: a car shows as 12 DDM copies in the same row, each about 3 bins wide (peak plus Hamming leakage), 12 x 3 = 36.

Conclusion for Phase 2:
- Most range rows are sparse, but those are the rows without a target. A row with a car has K about 36-41, well above the rough break-even of 8.
- So a plain sparse FFT mostly saves work that simply skipping empty rows would also save. To help where the cars are, Phase 3 has to use the DDM structure: 12 copies exactly 16 bins apart are one target, not 12 peaks.

Label range offset (open, does not affect sparsity):
- Not proportional to range (countryside at 42 m would be at about +6 if the bin size were wrong; it is at 0), and not a fixed shift (that would also show in the countryside).
- Labels come from lidar and camera, converted to radar coordinates; they are not the radar's own detection. Likely cause: a time gap between the lidar/camera snapshot and the radar frame paired with it. On the highway the car's range changes from frame to frame (70.0, 69.7, 69.1 m ...), so a small time gap becomes a range gap; in the countryside the cars barely change range.
- Possible sources: separate sensor clocks, label k vs radar frame k off by a frame or two, or the unresolved Valeo delay constants (radar -180 ms, lidar -40 ms).
- Explains the Step 0 highway gate miss at 70 m. Matters only for Phase 1 hit rates. Test later: compare highway labels against radar frames k-10 ... k+10.

Phase 2 done.

Decisions:
- Phase 3 implements the STOC 2012 exactly-sparse algorithm (Nearly Optimal Sparse Fourier Transform, Section 3), the phase method already learned.
- Moving the Phase 1-2 functions to src/ waits until cleanup (about Nov 24-26). The sparse FFT itself goes in src/sfft.py from its first line.

Next: read STOC 2012 Sections 1-3 (skip proofs; questions to answer: what breaks off-grid, why sigma is odd, why the flat window is needed, how large B must be for K = 36 by Claim 3.1, why subtracting found coefficients helps), then the flat window lesson, then src/sfft.py on a synthetic row.

Still to write here: exam end date and daily hours until Nov 19.

## Nov 14: Phase 3 theory, permutation and flat window

Studied (before reading STOC 2012 Sections 1-3):

Collisions with DDM:
- With plain buckets (take every 16th sample, 16-point FFT), bin f goes to bucket f mod 16. A car's 12 DDM copies are exactly 16 bins apart, so all 12 land in one bucket.
- The phase trick needs one tone per bucket. With 12 tones the bucket is a sum of 12 arrows; shifting by tau turns each arrow by its own angle, so the ratio's angle is a mix of 12 angles and points at no real copy (example: copies at 3, 19, ..., 179 gave angle 136 degrees -> bin 96.7, |ratio| 0.64 instead of 1). 12 tones = 24 unknowns, 2 measurements. Using all 16 groups would solve it but is the full FFT again.

Two ways to group bins into buckets (write the 256 bins as a 16 x 16 grid):
- Mod buckets = columns (bucket = f mod 16): every 16th bin. What subsampling gives. Fails on evenly spaced tones (DDM copies: one column).
- Contiguous buckets = rows (bucket = f // 16): 16 neighbouring bins, like a band-pass filter / one bar of a music equalizer. Fails on neighbouring tones (two cars at nearly the same speed, one copy's own leakage).
- The 12 DDM copies (bins 3 ... 179) are in one column but in 12 different rows, so contiguous buckets already separate them.

Permutation (sigma odd):
- Read the time samples in the order x[sigma*m mod N]. This moves every tone from bin f to bin sigma*f mod N. No FFT needed, only a different reading order (memory addressing on an FPGA). Only this kind of reordering keeps a tone a single tone; a random order would smear it.
- Undo after finding a tone: multiply the found position by sigma^-1 mod N (sigma = 5, N = 256: 205; sigma = 3, N = 16: 11, e.g. 15 * 11 mod 16 = 5).
- sigma must be odd: with N a power of 2, an even sigma maps two bins to one place (2 * 128 = 0) and cannot be undone.
- Permutation alone cannot fix mod buckets: (sigma*f) mod 16 depends only on f mod 16, so a column stays a column. It is used with contiguous buckets, where it scatters neighbours so two particular bins rarely share a row (paper's Lemma 2.4).
- Summary: mod buckets break on evenly spaced tones, contiguous buckets break on neighbours, permutation + contiguous buckets break on neither (most of the time).
- Time/frequency pairs in the paper's P(sigma, a, b): index x sigma <-> bins move to sigma*f; shift by a (= my tau) <-> each bin turns by its own angle; multiply by a rotating arrow (b) <-> all bins slide.

Flat window (how contiguous buckets are made from time samples):
- Recipe: multiply the samples by a window w, fold (cut into pieces of B samples and add them), one B-point FFT.
- Result 1: fold + B-point FFT gives exactly bins 0, L, 2L, ... of the big N-point FFT (L = N/B), because the piece number drops out of the exponent. Without a window a tone between those bins is missed.
- Result 2: multiplying by w in time smears the spectrum: Z[m] = (1/N) sum_f X[f] W[m-f].
- Together: bucket j = sum_f X[f] G[f - jL], with G[d] = W[-d]/N the band-pass filter shape. A flat G (1 over L neighbouring bins, 0 elsewhere) makes bucket j the sum of band j.
- A window is a recipe for which neighbouring bins get added: each rotating term in w adds a copy of the spectrum shifted by that many bins. Example N = 8, B = 4: w[n] = 1 + e^(-2 pi i n/8) adds a copy shifted down one bin, so Z[m] = X[m] + X[m+1] and bucket j = bins 2j and 2j+1.
- Worked example (N = 8, B = 4, tones at bin 2 size 1 and bin 5 size 0.5): true X = [0,0,8,0,0,4,0,0]; without window buckets [0,8,0,0] (bin 5 missed); with window buckets [0,8,4,0]; shift tau = 1: angles 90 deg -> bin 2, 225 deg -> bin 5.
- The phase trick survives the window: the same G multiplies both groups and cancels in the ratio, so a tone at a band edge is still located exactly (only its amplitude reads low, corrected by dividing by G).
- Cost: a perfectly flat band needs a window over all N samples. Rounding the band edges slightly (paper Definition 2.3, delta error) makes w short: for N = 256, B = 16 about 137 samples, not 16. So the saving at N = 256 is small; the paper's speed-up is about large N. With K about 36 per car row the saving may be small or none, and measuring that is the Phase 3 result.
- 4-point DFT note: e^(-2 pi i jm/4) = (-i)^(jm), since e^(-2 pi i/4) = -i.

Not fully clear yet, come back to it: the flat window derivation (Results 1 and 2 and how W is chosen from G).

Next: go over the flat window once more, then read STOC 2012 Sections 1-3 with the 5 questions, then src/sfft.py on a synthetic row.

## Nov 15-16: code moved to src/, sparse FFT Parts 1-4 built

Moved to src/sfftradar (no more notebook globals):
- radial.py: SEQS, Seq dataclass, load_seq(tag) returns a Seq (tag, name, data, radar, frames); read_block(seq, ...), chip_cube, frame_cube(seq, n) (the 4-chip read written once), load_labels().
- rd.py: W_RANGE, W_DOPPLER (Hamming), range_doppler(cube), frame_to_rd(seq, n), frame_to_pred(seq, n, mti), ca_cfar.
- sparsity.py: k_frac, sparsity_report, n_peaks, row_peaks, noise_map (noise through the same range_doppler).
- Notebooks start with sys.path.insert(0, "../src") and %autoreload 2. 02_sparsity.ipynb left as the Phase 2 record.
- Check after the move: city frame 0 K90 = 6929, noise K90 = 110919, per-row median 8.0 and K <= 8 share 0.502, all unchanged.

Decision: start coding before finishing the paper; understand the algorithm by testing each piece. Read STOC 2012 Section 3 afterwards.

src/sfftradar/sfft.py so far (tests in notebooks/03_sfft.ipynb, N = 256, B = 16):
1. tone(N, f, amp) and permute(x, sigma): read x[sigma*m mod N]. Tone at bin 37 moves to 185 = 5 x 37; sigma^-1 = pow(5, -1, 256) = 205, 185 x 205 mod 256 = 37.
2. flat_window(N, B, edge=1.5, thr=1e-3): design the band shape G in frequency (box of width 16 with erf-rounded edges), w = N * ifft(G), drop samples below 0.1% of the peak. Result: 141 samples used, flat part 0.977-1.000, leakage outside 0.0004. Trade-off: edge 1.0 -> 199 samples, flatness 0.998; edge 1.5 -> 141, 0.977; edge 2.0 -> 113, 0.933 (window length vs flatness; the paper's delta). Worth a sentence in the report.
3. hash_to_bins(x, sigma, tau, w, B): read x[sigma*(m + tau)], multiply by w, fold into B (reshape(N//B, B).sum(axis=0)), one B-point FFT. Matches the formula u[j] = sum_p Y[p] G[p - jL] to 6e-14. One car (12 DDM copies 3, 19, ..., 179): one copy per bucket in 12 buckets with sigma = 5 and also with sigma = 1, so the band-pass (contiguous) buckets do the separating; plain every-16th-sample buckets would put all 12 in one.
4. locate(x, sigma, w, B, thresh, ratio_tol=0.1): hash with tau = 0 and 1, r = u1/u0, permuted bin p = round(angle * N / 2 pi), f = p * sigma^-1 mod N, value = u0 / band height at p. Collision checks: |r| far from 1, and p outside bucket j's band (also rejects leaks of a tone into the neighbouring bucket).

Walk-through on the test row (tones 27, 66, 76, 105, 211; sigma = 5): permuted to 135, 74, 124, 13, 31; buckets 8, 5, 8, 1, 2. Round 1 finds 105, 211, 66 with values within 0.1%; 27 and 76 collide in bucket 8 and produce a fake tone at 229 (|r| = 1.019 passed the check); leaks of 66 and 27 into buckets 4 and 9 are rejected by the band check. All 12 car copies found in one round.

Preview (container run of Part 5, the full loop; my code not written yet): subtract found tones and repeat with a new random sigma; the 5-tone row is then recovered exactly and the fake cancels. Recovery rate over 200 random rows with B = 16: K = 1, 5: 100%; K = 8: about 97%; K = 12: about 56-64%; K = 20: about 0.5%. So B = 16 handles about 8 tones; a car row (K about 36) needs more buckets, which raises the cost.

Next: Part 5 (sparse_fft loop), Part 6 (operation counter), then the B vs K sweep, off-grid and noise tests, then real RADIal rows.
