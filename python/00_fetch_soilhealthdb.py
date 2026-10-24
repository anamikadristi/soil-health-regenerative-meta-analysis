"""Download SoilHealthDB_V1.xlsx (and the CSV copy) from the figshare record
cited in Jian, Du & Stewart (2020) Scientific Data 7:16.

Figshare record: SoilHealthDB, doi:10.6084/m9.figshare.8292176.v8 (CC0).
The record holds one 128 MB zip (SoilHealthDB.zip) that also contains large
global rasters we do not need. This script reads the zip's central directory
with HTTP range requests and extracts only the database files, so the
download is a few MB instead of 128 MB.

Standard library only. Run with:  python -I python/00_fetch_soilhealthdb.py data_raw
"""
import hashlib
import io
import json
import sys
import urllib.request
import zipfile
from pathlib import Path

ARTICLE_API = "https://api.figshare.com/v2/articles/8292176"
WANTED = ["SoilHealthDB/data/SoilHealthDB_V1.xlsx",
          "SoilHealthDB/data/SoilHealthDB_V1.csv",
          "SoilHealthDB/data/other_data/Koeppen-Geiger-ASCII.txt"]


class HTTPRangeFile(io.RawIOBase):
    """Minimal seekable read-only file object backed by HTTP range requests."""

    def __init__(self, url, size):
        self.url, self.size, self.pos = url, size, 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        if whence == 0:
            self.pos = offset
        elif whence == 1:
            self.pos += offset
        else:
            self.pos = self.size + offset
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        if n == 0 or self.pos >= self.size:
            return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        self.pos += len(data)
        return data

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(ARTICLE_API, timeout=60) as r:
        meta = json.load(r)
    f = [x for x in meta["files"] if x["name"] == "SoilHealthDB.zip"][0]
    log = [f"figshare record: {meta['title']} doi:{meta['doi']} license:{meta['license']['name']}",
           f"zip: {f['name']} size={f['size']} md5={f['computed_md5']} url={f['download_url']}"]
    # each range request goes through the figshare downloader (signed storage
    # URLs expire quickly); reads are buffered in 4 MB blocks
    zf = zipfile.ZipFile(io.BufferedReader(HTTPRangeFile(f["download_url"], f["size"]), buffer_size=4 * 2**20))
    for name in WANTED:
        info = zf.getinfo(name)
        data = zf.read(name)  # zipfile checks the CRC32 here
        dest = out / Path(name).name
        dest.write_bytes(data)
        log.append(f"extracted {name} ({info.file_size} bytes, CRC ok) sha256={hashlib.sha256(data).hexdigest()}")
    (out / "download_log.txt").write_text("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data_raw")
