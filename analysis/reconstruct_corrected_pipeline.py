#!/usr/bin/env python3
from pathlib import Path
import base64
import hashlib
import re
import zlib

EXPECTED_PARTS = {
    "000": "b356c9855eda38a6ef64f268e53ce2f752a55476ba5656b032c132790f2d7120",
    "001": "cccd8e582a577f98b8566f45f631a4e7f998b24eb9d5bff3776d96d41cb5e2b8",
    "002": "41617a144f7e8e3803e44fea8358a7ab28033b99fb4ae0b6678752071cc55e1e",
    "003": "1dc6443439bc635e8e8a8690dc2f6201a1998f6c9aad00f1c7a7052d2ca8480b",
    "004": "d3dce91500498fe0171ca2da43780135274aaa72077638d50ac8c81e0f383ae1",
    "005": "8ce23d88c3f9f4f591baeaefaa2fdf3c4bddabe36ada8aa7f5b0fb8a25a03e14",
    "006": "604ada3e90cf54f607072fdab5977b5accd989a2f443bfd34b17f761ece4f173",
    "007": "1389b4a80194f05ff3d123ea496235c4672608cdc438a6f786d37fa972ee52c6",
    "008": "b6ef19d31937672430d0a20c2229e0c7872da4070f7d4d510822a27a1db3c0e7",
    "009": "8d6cbd82e477940c126e25a58073a34d918d32269bd5cb2ffe02a9a4ee8814f2",
    "010": "cb0033c4e89cedaef8a07fdb47cc761f4231d79fc60d0012679a87157e56b3e7",
    "011": "566fed810df079960afb9cd109895e6aad74c4eb8ee9bc48994c7b73a09e6b3b",
    "012": "dc3e0c184b90e5600cf7bb64913a7cac4abe139daee48a3c3d86666854f3c534",
    "013": "5e5692e4b5923ae6f15a3b68c37e7b787c09719ada301ed64c0d642eae7ecd4c",
}
EXPECTED_SCRIPT = "f586edb19690c27800907c74da8198c1cca10f13e1e0b50e16e5aea618bd8aa0"

chunks = []
for i in range(14):
    key = f"{i:03d}"
    path = Path(f"analysis/payload_parts/part.{key}")
    text = re.sub(r"\s+", "", path.read_text(encoding="utf-8"))
    observed = hashlib.sha256(text.encode()).hexdigest()
    if observed != EXPECTED_PARTS[key]:
        raise RuntimeError(f"Payload part {key} hash mismatch: {observed}")
    chunks.append(text)
raw = zlib.decompress(base64.b64decode("".join(chunks)))
observed_script = hashlib.sha256(raw).hexdigest()
if observed_script != EXPECTED_SCRIPT:
    raise RuntimeError(f"Reconstructed script hash mismatch: {observed_script}")
out = Path("analysis/corrected_resubmission_pipeline.py")
out.write_bytes(raw)
print(f"Reconstructed {out} ({len(raw)} bytes; SHA256 {observed_script})")
