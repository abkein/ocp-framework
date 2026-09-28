On the input we have an ADIOS2 dump of an MD simulation of one
component plasma performed by LAMMPS. The dump is written with
`dump <ID> all custom/adios <PERIOD> dump.bp id x y z vx vy vz`.
The `atoms` array must have seven columns in the exact order
`ATOM_ID,X,Y,Z,VX,VY,VZ`; old five-column dumps are rejected.
Hereafter I use box=cell, referring to the simulation box.
We have cubic (guaranteed) cell and ion' positions.
We will designate the cell edge length as $L$ and ion number as $N$.
However just in case functions still use all $L_x, L_y, L_z$.
Input ion positions are wrapped into the box (cell).
Coordinate origin is exactly in the cell's center.

A sample LAMMPS input file can be found in [in.lammps](./in.lammps).

# Sphere selection and energies

Supply a positional list of multipliers after the input file, optionally
interspersed with CLI options. Each multiplier $f$ must be finite and in
$(0,3)$ and defines
$$
R_s=f\frac{L_{\min}}{2},\qquad
V_s=\frac{4\pi R_s^3}{3},\qquad
N_s=N\frac{V_s}{L_xL_yL_z}.
$$
The special token `c` means $f=\sqrt{3}$ and
$R_s=L_{\min}\sqrt{3}/2$, the circumsphere radius for a cube. The token is
converted to its numeric multiplier in the JSON output. For a noncubic box
it still uses $L_{\min}$, not the box diagonal.

- Multipliers $(0,1]$ use `inscribed()`: inspect all original ions and select
  those with $r_i^2<R_s^2$. Multiplier 1 gives the original inscribed sphere.
- Multipliers $(1,3)$ use `circumscribed()`: construct periodic images and
  select all points with $r_i^2<R_s^2$. Every selected image is a distinct
  point, including images of the same original ion. The necessary integer
  translations are restricted to $k_x,k_y,k_z\in\{-1,0,1\}$, so at most
  $3\times3\times3$ cells (the original plus 26 images) are processed.

Multiplier 3 is rejected. Every sphere excludes points exactly on its boundary,
including the cube corners for `c`. For `c`, squared-radius selection uses
$3(L_{\min}/2)^2$ directly to avoid rounding the square root in the boundary test.

Spheres are processed sequentially in the supplied order on every frame.
Repeated multipliers are permitted and receive separate output files.
All pair sums count each unordered pair once and use direct distances,
without a minimum-image reduction. The reported energies are
$$
u_c=\frac{3}{5}N_s^{2/3},\qquad
u_{p\mathrm{ex}}=\frac{1}{N_s}\sum_{i<j}\frac{1}{r_{ij}}-u_c,
$$
$$
u_{b\mathrm{ex}}=\frac12\left(\frac{1}{N_s}\sum_i r_i^2-u_c\right).
$$
Normalization uses the unrounded geometric mean $N_s$; the actual selected
count is reported separately. An empty selected sphere with $N>0$ has
$u_{p\mathrm{ex}}=-u_c$ and $u_{b\mathrm{ex}}=-u_c/2$. The formulas assume a
nonempty simulation ($N>0$).

# Virial and pressure contribution

For each sphere, write $N_s$ for its unrounded mean population and
$V_s=4\pi R_s^3/3$ for its volume. The density is
$n=N_s/V_s=N/V$, where $N,V$ refer to the original simulation cell.

Use the same reduced units as the energy formulas and `locp-reference.py`:
the pair potential is $1/r$, and the background force on an ion is
$\mathbf F_{b,i}=-\mathbf r_i$. Positions are measured from the box center.
For the physical OCP interpretation these are ion-sphere units
($a=1$, unit Coulomb coefficient, $k_B=1$); the program does not infer or
convert units from the dump.

For each selected ion or image, include forces only from the other selected
points and the background:
$$
\mathbf F_i=\sum_{\substack{j\in s\\j\ne i}}
\frac{\mathbf r_i-\mathbf r_j}{r_{ij}^3}-\mathbf r_i.
$$
The total force virial and the reported virial per mean ion are
$$
W_s=\sum_{i\in s}\mathbf r_i\cdot\mathbf F_i
=\sum_{\substack{i<j\\i,j\in s}}\frac{1}{r_{ij}}-\sum_{i\in s}r_i^2,
\qquad
w_s=\frac{W_s}{3N_s}
=\frac{u_{p\mathrm{ex}}^{(s)}-2u_{b\mathrm{ex}}^{(s)}}{3}.
$$
Thus the virial contribution to ion pressure is
$$
p_{\mathrm{ex},s}=\frac{W_s}{3V_s}=n w_s.
$$
The factor $1/3$ is already included in $w_s$. The constant energy offsets
$u_c$ cancel and do not contribute to forces. The implementation uses the
raw pair and squared-radius sums before subtracting these offsets.

These are sphere-restricted ion virials: interactions with points outside
the sphere are excluded, and no separate background self-pressure is added.
Total ion pressure additionally includes the kinetic term $nT$, following
the [virial pressure relation](https://docs.lammps.org/stable/compute_pressure.html).
These restricted values need not equal the full periodic Ewald pressure
reported by LAMMPS.

# Temperature

All ions have the same mass, supplied by `--mass m` (default 1). It must be
finite and positive, in the same simulation units as the velocities.
`--gamma` is no longer accepted.

For each frame, compute kinetic temperature from all $N$ original ions,
without counting periodic images:
$$
E_{\mathrm{kin}}=\frac{m}{2}\sum_{i=1}^{N}(v_{x,i}^2+v_{y,i}^2+v_{z,i}^2),
\qquad
T=\frac{2E_{\mathrm{kin}}}{3N-3}.
$$
This matches [LAMMPS's default three-dimensional `compute temp`](https://docs.lammps.org/compute_temp.html)
for the supplied script with no constraint fixes and $k_B=1$. The calculation
uses velocities as stored, without box scaling or center-of-mass subtraction;
the script removes initial center-of-mass motion. The degrees of freedom and
temperature are recalculated from the current atom count every frame.

The center-of-mass velocity components are also calculated from the original
ions in each frame:
$$
(c_{vx},c_{vy},c_{vz})=\frac{1}{N}\sum_{i=1}^{N}(v_{x,i},v_{y,i},v_{z,i}).
$$
The common mass cancels in total momentum divided by total mass. These signed
components are recorded as `cvx,cvy,cvz` in metadata, in simulation velocity
units. For one ion they equal that ion's velocity. They are diagnostics and
are not subtracted from velocities when calculating temperature.

The same current-frame temperature is used for every sphere:
$$
p_s=n(T+w_s)=nT+p_{\mathrm{ex},s},
\qquad
Z_{\mathrm{ion},s}=\frac{p_s V_s}{N_s T}
=\frac{p_s V}{NT}=1+\frac{w_s}{T}.
$$
Here $N_s$ is the same mean population used for energies and virials, not the
instantaneous selected count. Temperature is a whole-cell quantity, not a
separate estimate from the ions selected in each sphere.

For $N\le1$, the temperature convention is undefined: `T`, `p`, and `Z`
are reported as `nan`. For $N>1$ with zero velocities, `T=0`, pressure equals
`pex`, and `Z` is `nan` because its denominator is zero. Energies, virial,
and excess pressure are still calculated in both cases.

# Running and output

For example, calculate four numeric radii plus the exact circumsphere, with
ADIOS2 defaults and ion mass 1:

```sh
ocp-framework dump.bp 0.5 1 c 2 2.9 --config "" --mass 1
```

The flags `-i` and `-c` are replaced by the positional radius list. An RDF-only
run can omit the list; otherwise at least one sphere is required.

Each sphere has its own CSV: `observables-1.csv`, `observables-2.csv`, etc.,
numbered in request order. `--obs-out results.csv` changes the base name to
`results-1.csv`, etc.; a base without an extension gains `.csv`.
Every sphere file has these columns:

```text
timestep,r,Ns,Nrs,upex,ubex,uc,virial,pex,p,Z
```

`r` is the actual radius for that frame, recalculated as the box changes.
`Ns` is the mean population, `Nrs` the actual selected count, `virial` is
$w_s=W_s/(3N_s)$, and `pex` is the virial contribution to pressure.
An empty selected sphere with $N>0$ has zero virial and zero excess pressure.

`spheres.json` (configurable with `--spheres-out`) is a JSON array with one
entry per sphere:

```json
[{"rfrac": 0.5, "filename": "observables-1.csv"}]
```

Here `rfrac` is the numeric multiplier. Filenames are relative to the JSON
file's directory. The JSON contains no radius: each sphere CSV records `r`
for every frame. Currently `r = rfrac * metadata.rmax`; `rmax` remains in
metadata independently. The mapping is written once, including for an empty
trajectory. An RDF-only run produces an empty JSON array.

Each step also writes `timestep,N,V,rmax,nbins,T,cvx,cvy,cvz` to `metadata.csv`; use
`--meta-out path.csv` to change this filename. The files can be joined by
timestep. Neither metadata nor sphere files include `dr` or `Gamma`: `dr` is
the constant CLI bin width, and $\Gamma=1/T$. The metadata always reports
`rmax = 0.5 * min(Lx,Ly,Lz)`; `nbins` is zero when RDF calculation is disabled.

The measured temperature `T` is written to the metadata file. Each sphere
file reports total pressure and compressibility as `p,Z`, using that frame's
`T`. Virial and excess pressure are independent of the velocities and mass.
For positive temperature, zero virial gives $p_s=nT$ and $Z_{\mathrm{ion},s}=1$.

`--rdf dr` produces per-ion pair counts in radial bins of width `dr`, using
minimum-image distances. Distances at or above `rmax` are excluded, including
from a partially filled final bin. The range is fixed at half the shortest box
edge; there are no `--rmax` or `--rmax-frac` options. These counts are not
normalized by density and shell volume to produce $g(r)$.

In each step, `prepare_data()` reads the frame and completes the pending reads.
`run()` then calls `reader.EndStep()` immediately, followed by `modify_data()`
to convert coordinates, extract velocities, and calculate squared distances
from the box center,
then `calculate()` and `write()`.
