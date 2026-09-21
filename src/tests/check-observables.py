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


def check_virial(row, points, mean, density, suffix, gamma):
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
    close(row["virial" + suffix], virial, "force virial " + suffix)
    close(row["virial" + suffix], (float(row["upex" + suffix]) - 2 * float(row["ubex" + suffix])) / 3,
          "virial-energy identity " + suffix)
    close(row["pex" + suffix], total * density / (3 * mean), "excess pressure " + suffix)
    if gamma is None:
        for key in ("p" + suffix, "Z" + suffix):
            if not math.isnan(float(row[key])):
                raise AssertionError(f"{key} must be nan without --gamma")
    else:
        temperature = 1 / gamma
        # Construct pressure from the total force virial and the sphere volume.
        sphere_volume = mean / density
        pressure = (mean * temperature + total / 3) / sphere_volume
        close(row["p" + suffix], pressure, "total pressure " + suffix)
        close(row["Z" + suffix], pressure * sphere_volume / (mean * temperature), "compressibility " + suffix)


def check_frame(row, rdf, points, lengths, step, cutoff, dr, inscribed, circumscribed, gamma):
    if gamma is None:
        if not all(math.isnan(float(row[key])) for key in ("Gamma", "T")):
            raise AssertionError("Missing gamma and temperature must be reported as nan")
    else:
        close(row["Gamma"], gamma, "Gamma")
        close(row["T"], 1 / gamma, "T")
    for suffix, enabled in (("i", inscribed), ("c", circumscribed)):
        if not enabled:
            for key in ("virial" + suffix, "pex" + suffix, "p" + suffix, "Z" + suffix):
                if not math.isnan(float(row[key])):
                    raise AssertionError(f"Disabled quantity {key} must be nan")
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
        check_virial(row, inside, mean, len(points) / volume, "i", gamma)
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
        check_virial(row, inside, mean, len(points) / volume, "c", gamma)
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
             flags=("-i", "-c", "--rdf", "0.5"), cutoff=2, dr=0.5, gamma=None):
    dump = root / "dump.bp"
    header = [int(scaled), *lengths, *center]
    lines = [" ".join(map(str, header))]
    for frame in frames:
        lines.append(str(len(frame)))
        lines.extend(" ".join(map(str, p)) for p in frame)
    subprocess.run([writer, str(dump)], input="\n".join(lines), text=True, check=True)
    # An empty configuration path selects ADIOS2 defaults.
    gamma_args = [] if gamma is None else ["--gamma", str(gamma)]
    subprocess.run([program, str(dump), "--config", "", *flags, *gamma_args], cwd=root, check=True)
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
                    "-i" in flags, "-c" in flags, gamma)


program, writer = map(lambda s: str(Path(s).resolve()), sys.argv[1:])
rng = random.Random(17)
frames = [
    [(0, 0, 0)],  # Zero force virial and zero excess pressure.
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
    {"gamma": 5},
    {"gamma": 0.25, "scaled": True, "center": (3, -2, 7)},
    {"scaled": True, "center": (3, -2, 7)},
    {"center": (3, -2, 7)},
    {"flags": ("-c",), "gamma": 15},
    {"flags": ("-i",), "gamma": 2},
    {"flags": ("--rdf", "0.6", "--rmax", "1.7"), "dr": 0.6, "cutoff": 1.7},
    {"flags": ("--rdf", "0.5", "--rmax-frac", "0.9"), "gamma": 4},
    {"flags": ("--rdf", "0.5", "--rmax", "10")},
]
with tempfile.TemporaryDirectory(prefix="ocp-regression-") as directory:
    for index, case in enumerate(cases):
        root = Path(directory) / str(index)
        root.mkdir()
        run_case(root, frames, **case)
    root = Path(directory) / "rectangular"
    root.mkdir()
    run_case(root, [[(0, 0, 0), (0.4, 1.8, 2.9)]], lengths=(1, 4, 6), cutoff=0.5, gamma=3)
    for value in ("0", "-1", "nan", "inf", "-inf", "1oops", "abc", "1e999", "1e-999"):
        result = subprocess.run([program, "unused.bp", "-i", "--gamma", value],
                                capture_output=True, text=True)
        if result.returncode == 0 or "--gamma must be finite and > 0" not in result.stderr:
            raise AssertionError(f"Invalid gamma accepted or misdiagnosed: {value}")
    result = subprocess.run([program, "unused.bp", "-i", "--gamma"], capture_output=True, text=True)
    if result.returncode == 0 or "Missing value after --gamma" not in result.stderr:
        raise AssertionError("Missing gamma value accepted or misdiagnosed")
print(f"Passed {len(cases) * len(frames) + 1} frame checks and 10 invalid-gamma checks.")
