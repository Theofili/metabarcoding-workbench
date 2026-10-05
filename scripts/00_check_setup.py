import importlib.util
import shutil
import sys

print("Python:", sys.version.split()[0])

packages = ["Bio", "sklearn", "matplotlib", "pandas"]
for name in packages:
    ok = importlib.util.find_spec(name) is not None
    print(f"{name}: {'OK' if ok else 'MISSING'}")

print("Rscript:", shutil.which("Rscript") or "NOT FOUND")

missing = [n for n in packages if importlib.util.find_spec(n) is None]
if missing:
    print("\nMissing Python packages:", ", ".join(missing))
    sys.exit(1)

print("\nPython setup looks good.")
