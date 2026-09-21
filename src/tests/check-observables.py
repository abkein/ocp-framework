"""Exercise the executable and ADIOS2 reader against independent small-system sums."""

import csv
import itertools
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


def check_frame(row, rdf, points, lengths, step, cutoff, dr, inscribed, circumscribed):
    close(row["timestep"], step, "timestep")
    close(row["N"], len(points), "N")
    volume = math.prod(lengths)
    close(row["V"], volume, "volume")
    if inscribed:
        radius = min(lengths) / 2
        inside = [p for p in points if math.fsum(v * v for v in p) <= radius**2]
        mean = len(points) * 4 * math.pi * radius**3 / (3 * volume)
        close(row["Ni"], mean, "Ni")
        close(row["Nri"], len(inside), "Nri")
        for key, value in zip(("upexi", "ubexi", "uci"), observables(inside, mean)):
            close(row[key], value, key)
    if circumscribed:
        radius_sq = sum(v * v for v in lengths) / 4
        # Deliberately oversized fixed lattice: independent of production bounds.
        images = [
            tuple(p[a] + shift[a] * lengths[a] for a in range(3))
            for p in points
            for shift in itertools.product(range(-5, 6), repeat=3)
        ]
        inside = [p for p in images if sum(v * v for v in p) <= radius_sq]
        mean = len(points) * 4 * math.pi * radius_sq**1.5 / (3 * volume)
        close(row["Nc"], mean, "Nc")
        close(row["Nrc"], len(inside), "Nrc")
        for key, value in zip(("upexc", "ubexc", "ucc"), observables(inside, mean)):
            close(row[key], value, key)
    if rdf is not None:
        close(row["rmax"], cutoff, "rmax")
        bins = math.ceil(cutoff / dr)
        close(row["nbins"], bins, "nbins")
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


def run_case(root, frames, lengths=(4, 4, 4), center=(0, 0, 0), scaled=False,
             flags=("-i", "-c", "--rdf", "0.5"), cutoff=2, dr=0.5):
    dump = root / "dump.bp"
    header = [int(scaled), *lengths, *center]
    lines = [" ".join(map(str, header))]
    for frame in frames:
        lines.append(str(len(frame)))
        lines.extend(" ".join(map(str, p)) for p in frame)
    subprocess.run([writer, str(dump)], input="\n".join(lines), text=True, check=True)
    # An empty configuration path selects ADIOS2 defaults.
    subprocess.run([program, str(dump), "--config", "", *flags], cwd=root, check=True)
    with (root / "observables.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(frames):
        raise AssertionError("Unexpected frame count")
    rdf = None
    if "--rdf" in flags:
        with (root / "rdf.csv").open() as stream:
            rdf = list(csv.reader(stream))
        if len(rdf) != len(frames):
            raise AssertionError("Unexpected RDF frame count")
    for step, (row, frame) in enumerate(zip(rows, frames)):
        check_frame(row, rdf[step] if rdf else None, frame, lengths, step, cutoff, dr,
                    "-i" in flags, "-c" in flags)


program, writer = map(lambda s: str(Path(s).resolve()), sys.argv[1:])
rng = random.Random(17)
frames = [
    [(-0.5, 0, 0), (0.5, 0, 0), (1.9, 1.9, 1.9), (-1.9, -1.9, -1.9)],
    [(1.9, 1.9, 1.9)],  # Empty inscribed sphere, but nonzero mean population.
    [(-2, -2, -2)],  # Eight boundary images, including the corner neighbor.
    [(-1.9, -1.9, 0)],  # Edge images; pairs of images of the same ion.
    [(-2, 0, 0), (0, 0, 0)],  # Inscribed boundary and RDF's excluded upper edge.
    [(-0.75, -0.75, -0.75), (0.75, 0.75, 0.75)],  # Entirely beyond RDF cutoff.
    [(rng.uniform(-1.9, 1.9), rng.uniform(-1.9, 1.9), rng.uniform(-1.9, 1.9)) for _ in range(9)],
]
cases = [
    {},
    {"scaled": True, "center": (3, -2, 7)},
    {"center": (3, -2, 7)},
    {"flags": ("-c",)},
    {"flags": ("-i",)},
    {"flags": ("--rdf", "0.6", "--rmax", "1.7"), "dr": 0.6, "cutoff": 1.7},
    {"flags": ("--rdf", "0.5", "--rmax-frac", "0.9")},
    {"flags": ("--rdf", "0.5", "--rmax", "10")},
]
with tempfile.TemporaryDirectory(prefix="ocp-regression-") as directory:
    for index, case in enumerate(cases):
        root = Path(directory) / str(index)
        root.mkdir()
        run_case(root, frames, **case)
    root = Path(directory) / "rectangular"
    root.mkdir()
    run_case(root, [[(0, 0, 0), (0.4, 1.8, 2.9)]], lengths=(1, 4, 6), cutoff=0.5)
print("Passed 57 frame checks, including scaled coordinates, periodic images, and RDF cutoffs.")
