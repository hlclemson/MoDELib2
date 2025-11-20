import os
import shutil
import subprocess
import tempfile
from pathlib import Path

# -------------------------------------------------
# Parameters you may want to tweak
# -------------------------------------------------
SRC_DIR   = Path.cwd()          # folder that contains img_*.png
DST_FILE  = Path("movie.mp4")    # output file
FPS       = 1                   # frames per second
CRF       = 20                   # quality (lower = larger file, 18-28 is sane)

# -------------------------------------------------
# 1. Collect and sort PNGs numerically
# -------------------------------------------------
pngs = sorted(
    SRC_DIR.glob("img_*.png"),
    key=lambda p: int(p.stem.split("_")[1])
)

if not pngs:
    raise RuntimeError("No img_*.png files found!")

# -------------------------------------------------
# 2. Copy into a temp dir with zero-padded names
#    (FFmpeg needs lexicographic order)
# -------------------------------------------------
tmp = tempfile.TemporaryDirectory()
tmp_path = Path(tmp.name)

for idx, src in enumerate(pngs, start=1):
    dst = tmp_path / f"{idx:04d}.png"
    shutil.copy2(src, dst)

# -------------------------------------------------
# 3. Call FFmpeg
# -------------------------------------------------
ffmpeg_cmd = [
    "ffmpeg",
    "-y",                      # overwrite output
    "-framerate", str(FPS),
    "-i", str(tmp_path / "%04d.png"),
    "-c:v", "libx264",
    "-pix_fmt", "yuv420p",
    "-crf", str(CRF),
    str(DST_FILE)
]

print("Running:", " ".join(ffmpeg_cmd))
subprocess.run(ffmpeg_cmd, check=True)

tmp.cleanup()
print(f"Done → {DST_FILE.absolute()}")

