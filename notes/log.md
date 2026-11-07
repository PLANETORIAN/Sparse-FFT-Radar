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

