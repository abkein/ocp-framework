"""Exercise the executable and ADIOS2 reader against independent small-system sums."""

import csv
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def close(actual, expected, label):
    if not math.isclose(float(actual), expected, rel_tol=2e-12, abs_tol=2e-12):
        raise AssertionError(f"{label}: got {actual}, expected {expected}")


def observables(points, mean):
    # Ordered pairs with an explicit 1/2, independently of the C++ i<j loop.
    pair = math.fsum(
        1 / math.dist(a, b)
        for i, a in enumerate(points)
        for j, b in enumerate(points)
        if i != j
    ) / 2
    radial = math.fsum(v * v for p in points for v in p)
    uc = 0.6 * mean ** (2 / 3)
    return pair / mean - uc, (radial / mean - uc) / 2, uc


def check_virial(row, points, mean, density, gamma):
    # Independent reference: build every ion's vector force and dot it with r.
    # The production code instead reuses the scalar energy sums.
    total = 0.0
    for i, point in enumerate(points):
        force = [-v for v in point]  # Uniform-background force.
        for j, other in enumerate(points):
            if i == j:
                continue
            r = math.dist(point, other)
            for axis in range(3):
                force[axis] += (point[axis] - other[axis]) / r**3
        total += sum(v * f for v, f in zip(point, force))
    virial = total / (3 * mean)
    close(row["virial"], virial, "force virial ")
    close(row["virial"], (float(row["upex"]) - 2 * float(row["ubex"])) / 3,
          "virial-energy identity ")
    close(row["pex"], total * density / (3 * mean), "excess pressure ")
    if gamma is None:
        for key in ("p", "Z"):
            if not math.isnan(float(row[key])):
                raise AssertionError(f"{key} must be nan without --gamma")
    else:
        temperature = 1 / gamma
        # Construct pressure from the total force virial and the sphere volume.
        sphere_volume = mean / density
        pressure = (mean * temperature + total / 3) / sphere_volume
        close(row["p"], pressure, "total pressure ")
        close(row["Z"], pressure * sphere_volume / (mean * temperature), "compressibility ")


def check_frame(metadata, rdf, points, lengths, step, dr, gamma):
    cutoff = min(lengths) / 2
    close(metadata["timestep"], step, "metadata timestep")
    close(metadata["N"], len(points), "N")
    close(metadata["V"], math.prod(lengths), "volume")
    close(metadata["rmax"], cutoff, "rmax")
    close(metadata["nbins"], math.ceil(cutoff / dr) if rdf is not None else 0, "nbins")
    if gamma is None:
        if not math.isnan(float(metadata["T"])):
            raise AssertionError("Missing temperature must be reported as nan")
    else:
        close(metadata["T"], 1 / gamma, "T")
    if rdf is not None:
        bins = math.ceil(cutoff / dr)
        expected = [0.0] * bins
        for a, b in itertools.combinations(points, 2):
            distance = math.sqrt(sum(
                (d - length * math.floor(d / length + 0.5)) ** 2
                for d, length in zip((a[i] - b[i] for i in range(3)), lengths)
            ))
            if distance < cutoff:
                expected[int(distance / dr)] += 2 / len(points)
        close(rdf[0], step, "rdf timestep")
        if len(rdf) != bins + 1:
            raise AssertionError("Unexpected RDF column count")
        for actual, wanted in zip(rdf[1:], expected):
            close(actual, wanted, "rdf bin")


def check_sphere(row, points, lengths, step, rfrac, gamma):
    base_radius = min(lengths) / 2
    factor = math.sqrt(3) if rfrac == "c" else float(rfrac)
    radius = base_radius * factor
    radius_sq = 3 * base_radius**2 if rfrac == "c" else radius**2
    candidates = points
    if factor > 1:
        # Deliberately oversized lattice independent of production image bounds.
        candidates = [
            tuple(p[a] + shift[a] * lengths[a] for a in range(3))
            for p in points
            for shift in itertools.product(range(-3, 4), repeat=3)
        ]
    inside = [p for p in candidates if sum(v * v for v in p) < radius_sq]
    mean = len(points) * 4 * math.pi * radius**3 / (3 * math.prod(lengths))
    close(row["timestep"], step, "sphere timestep")
    close(row["r"], radius, "per-frame sphere radius")
    close(row["Ns"], mean, "mean population")
    close(row["Nrs"], len(inside), "selected population")
    for key, value in zip(("upex", "ubex", "uc"), observables(inside, mean)):
        close(row[key], value, key)
    check_virial(row, inside, mean, len(points) / math.prod(lengths), gamma)


def run_case(root, frames, lengths=(4, 4, 4), center=(0, 0, 0), scaled=False,
             flags=("--rdf", "0.5"), dr=0.5, gamma=None, metadata_out=None,
             radii=("0.5", "1", "c", "2", "2.9"), obs_out="observables.csv",
             spheres_out="spheres.json", frame_lengths=None):
    dump = root / "dump.bp"
    header = [int(scaled), *lengths, *center]
    lines = [" ".join(map(str, header))]
    for index, frame in enumerate(frames):
        lines.append(str(len(frame)))
        if frame_lengths is not None:
            lines.append(" ".join(map(str, frame_lengths[index])))
        lines.extend(" ".join(map(str, p)) for p in frame)
    writer_args = ["--vary-box"] if frame_lengths is not None else []
    subprocess.run([writer, str(dump), *writer_args], input="\n".join(lines), text=True, check=True)
    gamma_args = [] if gamma is None else ["--gamma", str(gamma)]
    metadata_args = [] if metadata_out is None else ["--meta-out", metadata_out]
    (root / obs_out).parent.mkdir(parents=True, exist_ok=True)
    (root / spheres_out).parent.mkdir(parents=True, exist_ok=True)
    # Place radii after other options to exercise interspersed option parsing.
    subprocess.run([program, str(dump), "--config", "", *flags, *gamma_args, *radii,
                    *metadata_args, "--obs-out", obs_out, "--spheres-out", spheres_out], cwd=root, check=True)
    manifest_path = root / spheres_out
    manifest = json.loads(manifest_path.read_text())
    if len(manifest) != len(radii):
        raise AssertionError("Unexpected sphere count in manifest")
    sphere_rows = []
    filenames = set()
    for index, (entry, rfrac) in enumerate(zip(manifest, radii)):
        if set(entry) != {"rfrac", "filename"}:
            raise AssertionError("Unexpected manifest fields")
        factor = math.sqrt(3) if rfrac == "c" else float(rfrac)
        close(entry["rfrac"], factor, "manifest multiplier")
        filename = (manifest_path.parent / entry["filename"]).resolve()
        base = root / obs_out
        expected_name = base.with_name(f"{base.stem}-{index + 1}{base.suffix or '.csv'}")
        if filename != expected_name.resolve() or filename in filenames:
            raise AssertionError("Incorrect or duplicate sphere filename")
        filenames.add(filename)
        with filename.open() as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != "timestep,r,Ns,Nrs,upex,ubex,uc,virial,pex,p,Z".split(","):
                raise AssertionError(f"Unexpected sphere columns: {reader.fieldnames}")
            rows = list(reader)
        if len(rows) != len(frames):
            raise AssertionError("Unexpected sphere frame count")
        sphere_rows.append(rows)
    if (root / obs_out).exists():
        raise AssertionError("Unexpected combined observables file")
    with (root / (metadata_out or "metadata.csv")).open() as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["timestep", "N", "V", "rmax", "nbins", "T"]:
            raise AssertionError(f"Unexpected metadata columns: {reader.fieldnames}")
        metadata_rows = list(reader)
    if len(metadata_rows) != len(frames):
        raise AssertionError("Unexpected metadata frame count")
    if metadata_out is not None and (root / "metadata.csv").exists():
        raise AssertionError("Default metadata file created despite --meta-out")
    rdf = None
    if "--rdf" in flags:
        with (root / "rdf.csv").open() as stream:
            rdf = list(csv.reader(stream))
        if len(rdf) != len(frames):
            raise AssertionError("Unexpected RDF frame count")
    for step, frame in enumerate(frames):
        current_lengths = frame_lengths[step] if frame_lengths is not None else lengths
        check_frame(metadata_rows[step], rdf[step] if rdf else None, frame, current_lengths, step, dr, gamma)
        for rfrac, rows in zip(radii, sphere_rows):
            check_sphere(rows[step], frame, current_lengths, step, rfrac, gamma)


program, writer = map(lambda s: str(Path(s).resolve()), sys.argv[1:])
rng = random.Random(17)
frames = [
    [(0, 0, 0)],  # For rfrac=2, the six images at radius L are excluded.
    [(1, 0, 0), (0, 0, 0)],  # rfrac=0.5 excludes the ion exactly at radius 1.
    [(-0.5, 0, 0), (0.5, 0, 0), (1.9, 1.9, 1.9), (-1.9, -1.9, -1.9)],
    [(1.9, 1.9, 1.9)],  # Empty inscribed sphere, but nonzero mean population.
    [(-2, -2, -2)],  # All eight corner images must be excluded for c.
    [(-1.9, -1.9, 0)],  # Edge images; pairs of images of the same ion.
    [(-2, 0, 0), (0, 0, 0)],  # Inscribed boundary and RDF's excluded upper edge.
    [(-0.75, -0.75, -0.75), (0.75, 0.75, 0.75)],  # Entirely beyond RDF cutoff.
    [(rng.uniform(-1.9, 1.9), rng.uniform(-1.9, 1.9), rng.uniform(-1.9, 1.9)) for _ in range(9)],
]
cases = [
    {},
    {"gamma": 5, "metadata_out": "frames.csv", "obs_out": 'outputs/sphere "quoted"\\name.csv',
     "spheres_out": "indices/spheres.json"},
    {"gamma": 0.25, "scaled": True, "center": (3, -2, 7)},
    {"scaled": True, "center": (3, -2, 7)},
    {"center": (3, -2, 7)},
    {"flags": (), "radii": ("c",), "gamma": 15},
    {"flags": (), "radii": ("1",), "gamma": 2},
    {"flags": ("--rdf", "0.6"), "dr": 0.6, "radii": ()},
    {"flags": ("--rdf", "0.5"), "gamma": 4},
    {"flags": ("--rdf", "3"), "dr": 3, "radii": (str(math.nextafter(3.0, 0.0)), "1", "2.9", "0.5", "c")},
]
with tempfile.TemporaryDirectory(prefix="ocp-regression-") as directory:
    for index, case in enumerate(cases):
        root = Path(directory) / str(index)
        root.mkdir()
        run_case(root, frames, **case)
    root = Path(directory) / "rectangular"
    root.mkdir()
    run_case(root, [[(0, 0, 0), (0.4, 1.8, 2.9)]], lengths=(1, 4, 6), gamma=3)
    for value in ("0", "-1", "nan", "inf", "-inf", "1oops", "abc", "1e999", "1e-999"):
        result = subprocess.run([program, "unused.bp", "1", "--gamma", value],
                                capture_output=True, text=True)
        if result.returncode == 0 or "--gamma must be finite and > 0" not in result.stderr:
            raise AssertionError(f"Invalid gamma accepted or misdiagnosed: {value}")
    result = subprocess.run([program, "unused.bp", "1", "--gamma"], capture_output=True, text=True)
    if result.returncode == 0 or "Missing value after --gamma" not in result.stderr:
        raise AssertionError("Missing gamma value accepted or misdiagnosed")
    for flag in ("--rmax", "--rmax-frac", "-i", "-c"):
        result = subprocess.run([program, "unused.bp", "--rdf", "0.5", flag, "1"],
                                capture_output=True, text=True)
        if result.returncode == 0 or f"Unknown argument: {flag}" not in result.stderr:
            raise AssertionError(f"Removed option still accepted or misdiagnosed: {flag}")
    for value in ("0", "3", "3.0001", "nan", "inf", "abc", "1oops", "1e999", "1e-999"):
        result = subprocess.run([program, "unused.bp", value], capture_output=True, text=True)
        if result.returncode == 0 or "Sphere multiplier must be" not in result.stderr:
            raise AssertionError(f"Invalid multiplier accepted or misdiagnosed: {value}")
    for args in ([], ["-0.5"]):
        result = subprocess.run([program, "unused.bp", *args], capture_output=True, text=True)
        expected = "Need at least one sphere" if not args else "Unknown argument: -0.5"
        if result.returncode == 0 or expected not in result.stderr:
            raise AssertionError(f"Invalid sphere arguments accepted or misdiagnosed: {args}")
    root = Path(directory) / "varying-box"
    root.mkdir()
    run_case(root, [[(0, 0, 0), (0.5, 0, 0)]] * 3, scaled=True, gamma=2,
             frame_lengths=[(4, 4, 4), (6, 6, 6), (3, 3, 3)])
    root = Path(directory) / "fractional-box"
    root.mkdir()
    run_case(root, [[(-1.3, -1.3, -1.3)], [(-1.3, 0, 0)]], lengths=(2.6, 2.6, 2.6), gamma=1)
print(f"Passed {len(cases) * len(frames) + 6} frame checks across multiple radii and 25 invalid-argument checks.")
