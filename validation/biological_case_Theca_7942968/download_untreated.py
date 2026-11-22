"""Download and verify 10 untreated human theca-cell 10x archives from Zenodo."""
import concurrent.futures
import hashlib
import json
import urllib.request
from pathlib import Path

record_url = 'https://zenodo.org/api/records/7942968'
dest = Path('../RawData/Zenodo7942968_Theca')
dest.mkdir(parents=True, exist_ok=True)
with urllib.request.urlopen(record_url, timeout=30) as response:
    record = json.load(response)
files = sorted((x for x in record['files'] if x['key'].endswith('-C.tar')),
               key=lambda x: x['key'])
if len(files) != 10:
    raise RuntimeError(f'Expected 10 untreated archives, got {len(files)}')
(dest / 'zenodo_record.json').write_text(json.dumps(record, indent=2),encoding='utf-8')

def download(item):
    path = dest / item['key']
    checksum = item['checksum'].split(':',1)[-1]
    if path.exists():
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        if digest == checksum:
            return f'{path.name}: already verified'
    req = urllib.request.Request(item['links']['self'],
                                 headers={'User-Agent':'PCOS-CAM academic analysis'})
    with urllib.request.urlopen(req,timeout=120) as response, path.open('wb') as handle:
        while True:
            block = response.read(1024*1024)
            if not block:
                break
            handle.write(block)
    digest = hashlib.md5(path.read_bytes()).hexdigest()
    if digest != checksum:
        path.unlink(missing_ok=True)
        raise RuntimeError(f'{path.name}: checksum mismatch')
    return f'{path.name}: verified ({path.stat().st_size} bytes)'

with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for line in pool.map(download,files):
        print(line,flush=True)
