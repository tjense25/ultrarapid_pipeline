#!/usr/bin/env bash
set -euo pipefail

python - <<'PY'
from pathlib import Path
import shutil

root = Path(shutil.which("clair3.py")).resolve().parent
dispatcher = root / "clair3" / "CallVariantsFromCffiGPU.py"
worker = root / "clair3" / "CallVariantsFromCffi.py"

d = dispatcher.read_text()
d = d.replace(
    '''        device_ids = [int(x) for x in args.device.split("cuda:")[-1].split(",")]
        os.environ["CUDA_VISIBLE_DEVICES"] = ','.join([str(item) for item in device_ids])''',
    '''        device_ids = [int(x) for x in args.device.split("cuda:")[-1].split(",")]'''
)
d = d.replace("    for device_id in all_device_ids:\n", "    for device_id in device_ids:\n", 1)
dispatcher.write_text(d)

w = worker.read_text()
old = '''    elif args.use_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu_id) if args.gpu_id else ""
        device = _select_device(True)'''
new = '''    elif args.use_gpu:
        if args.gpu_id is None:
            raise RuntimeError("GPU mode requested without --gpu_id")
        torch.cuda.set_device(args.gpu_id)
        device = torch.device(f"cuda:{args.gpu_id}")'''

if old not in w:
    raise RuntimeError("Expected Clair3 2.0.3 GPU code not found")

worker.write_text(w.replace(old, new))
print("Applied Clair3 2.0.3 GPU-selection patch")
PY
