"""Reading raw RADIal recordings: radar block index, ADC blocks, frame cubes, labels."""
import os
from dataclasses import dataclass
import numpy as np
import pandas as pd

ROOT = r"C:\Users\daksh\Documents\SPARSE_FFT\data"
SEQS = {
    "city":        "RECORD@2020-11-21_13.46.26",
    "countryside": "RECORD@2020-11-22_12.54.38",
    "highway":     "RECORD@2020-11-22_12.38.47",
}
CHIPS = ("radar_ch0", "radar_ch1", "radar_ch2", "radar_ch3")


@dataclass
class Seq:
    tag: str        # "city", "countryside", "highway"
    name: str       # "RECORD@..."
    data: str       # folder with the .bin and .rec files
    radar: dict     # chip -> {sample number: (offset, size)}
    frames: list    # ch3 sample numbers with all 4 chips present


def load_seq(tag, root=ROOT):
    """Index one sequence's radar blocks and list its complete frames."""
    name = SEQS[tag]
    data = os.path.join(root, name, name)
    log_path = os.path.join(data, name + "_events_log.rec")

    radar = {chip: {} for chip in CHIPS}
    with open(log_path) as f:
        for line in f:
            pieces = line.split()
            if len(pieces) != 17 or pieces[10] not in radar:
                continue
            radar[pieces[10]][int(pieces[7])] = (int(pieces[13]), int(pieces[16]))

    frames = [n for n in radar["radar_ch3"]
              if (n + 1) in radar["radar_ch0"]
              and (n + 2) in radar["radar_ch2"]
              and (n + 3) in radar["radar_ch1"]]
    return Seq(tag, name, data, radar, frames)


def read_block(seq, chip, offset, datasize):
    """Read one chip's raw int16 ADC block from its .bin file."""
    path = os.path.join(seq.data, seq.name + f"_{chip}.bin")
    with open(path, "rb") as f:
        f.seek(offset)
        return np.fromfile(f, dtype=np.int16, count=datasize // 2)


def chip_cube(x):
    """Raw int16 block -> complex (512 samples, 256 chirps, 4 Rx)."""
    c = x[0::2] + 1j * x[1::2]
    c = np.reshape(c, (512, 4, 256), order="F")
    return c.transpose((0, 2, 1))


def frame_cube(seq, n):
    """Raw (512, 256, 16) cube for frame n: ch3 n, ch0 n+1, ch1 n+3, ch2 n+2."""
    keys = [("radar_ch3", n), ("radar_ch0", n + 1), ("radar_ch1", n + 3), ("radar_ch2", n + 2)]
    cubes = []
    for chip, s in keys:
        off, size = seq.radar[chip][s]
        cubes.append(chip_cube(read_block(seq, chip, off, size)))
    return np.concatenate(cubes, axis=2)


def load_labels(root=ROOT):
    """All RADIal labels (labels_CVPR.csv)."""
    return pd.read_csv(os.path.join(root, "labels_CVPR.csv"))