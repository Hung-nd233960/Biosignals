import mne
import matplotlib

matplotlib.use("TkAgg")
raw = mne.io.read_raw_edf(
    "project/data/sub-NORB00055_ses-1_task-EEG_eeg.edf", preload=True
)

print(raw.info)  # metadata
print(raw.ch_names)  # channel names
for ch in raw.info["chs"]:
    print(ch["ch_name"], ch["unit"], ch["kind"])
print(raw.annotations)
raw.plot(block=True)  # visualize signals
