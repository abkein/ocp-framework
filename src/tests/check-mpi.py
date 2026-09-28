"""Compare MPI output with the numerical reference and one-rank output."""

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile


sys.dont_write_bytecode = True

spec = importlib.util.spec_from_file_location("reference", Path(__file__).with_name("check-observables.py"))
reference = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reference)
reference.program, reference.writer = (str(Path(arg).resolve()) for arg in sys.argv[1:3])
mpiexec, count_flag = sys.argv[3:5]


def launcher(ranks):
    return [mpiexec, count_flag, str(ranks)]


frames = [
    [(-0.5, 0, 0), (0.5, 0, 0)],
    [(0, 0, 0)],
    [(-0.75, 0.25, 0), (0.75, -0.25, 0), (0, 0, 1)],
    [(-1, 0, 0), (1, 0, 0)],
    [(0.1, 0.2, 0.3), (1.1, 0.2, 0.3)],
]
lengths = [(4, 4, 4), (6, 6, 6), (3, 3, 3), (5, 5, 5), (4, 4, 4)]

with tempfile.TemporaryDirectory(prefix="ocp-mpi-") as directory:
    root = Path(directory)
    serial = root / "serial"
    serial.mkdir()
    expected = reference.run_case(serial, frames, frame_lengths=lengths, scaled=True, mass=2,
                                  radii=("0.5", "c", "2.9"))
    for ranks in (2, 3, 4):
        path = root / f"ranks-{ranks}"
        path.mkdir()
        actual = reference.run_case(path, frames, frame_lengths=lengths, scaled=True, mass=2,
                                    radii=("0.5", "c", "2.9"), launcher=launcher(ranks))
        if actual != expected:
            raise AssertionError(f"{ranks}-rank output differs from one-rank output")
        # Repeat with fewer frames in the same files: verifies truncation and
        # idle ranks, including ranks that never own a timestep.
        short = reference.run_case(path, frames[:1], scaled=True, mass=2,
                                   radii=("0.5", "c", "2.9"), launcher=launcher(ranks))
        if short != {key: [rows[:1] for rows in value] if key == "spheres" else value[:1]
                     for key, value in expected.items()}:
            raise AssertionError("A rerun left stale rows or lost the first frame")

    rdf_only = root / "rdf-only"
    rdf_only.mkdir()
    reference.run_case(rdf_only, frames, radii=(), launcher=launcher(3))
    spheres_only = root / "spheres-only"
    spheres_only.mkdir()
    reference.run_case(spheres_only, frames, flags=(), radii=("1",), launcher=launcher(2))

    failure = root / "rank-local-failure"
    failure.mkdir()
    dump = failure / "dump.bp"
    frame = "2\n-0.5 0 0 1 2 3\n0.5 0 0 -1 -2 -3\n"
    subprocess.run([reference.writer, str(dump), "--bad-second-step"],
                   input="0 4 4 4 0 0 0\n" + frame * 2, text=True, check=True)
    result = subprocess.run([*launcher(2), reference.program, str(dump), "1", "--config", ""],
                            cwd=failure, capture_output=True, text=True, timeout=30)
    if result.returncode == 0 or "Rank 1 error:" not in result.stderr or "Expected atoms shape [natoms,7]" not in result.stderr:
        raise AssertionError(f"Rank-local read failure was not reported cleanly: {result.stderr}")

    # Only rank zero writes the JSON. Its failure must also stop peers that
    # may already be waiting inside the next collective file open.
    result = subprocess.run([*launcher(2), reference.program, str(serial / "dump.bp"), "1", "--config", "",
                             "--spheres-out", "missing/spheres.json"],
                            cwd=failure, capture_output=True, text=True, timeout=30)
    if result.returncode == 0 or "Cannot open missing/spheres.json" not in result.stderr:
        raise AssertionError("Rank-zero output failure was not reported cleanly")

print("Passed 2/3/4-rank equivalence, partial batches, idle ranks, truncation, and rank-local failure checks.")
