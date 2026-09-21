On the input we have an ADIOS2 dump of an MD simulation of one
component plasma performed by LAMMPS. The dump is written with
`dump <ID> all atom/adios <PERIOD> dumps/dump.bp`.
Hereafter I use box=cell, referring to the simulation box.
We have cubic (guaranteed) cell and ion' positions.
We will designate the cell edge length as $L$ and ion number as $N$.
However just in case functions still use all $L_x, L_y, L_z$.
Input ion positions are wrapped into the box (cell).
Coordinate origin is exactly in the cell's center.

A sample LAMMPS input file can be found in [in.lammps](./in.lammps).

# Function `inscribed`:

We imagine a sphere, inscribed into the cube (cubic cell).
The sphere' radius is then $R_{\mathrm{insc}}=L/2$.
The number of ions in the sphere (average) $N_{\mathrm{insc}}=N \pi/6$.
What function should do is for ions inside the inscribed sphere calculate:
$$
u^{(\mathrm{insc})}_{c} = \frac{3}{5}N_{\mathrm{insc}}^{2/3},
$$
$$
u^{(\mathrm{insc})}_{p\mathrm{ex}} = \frac{1}{N_{\mathrm{insc}}}\sum_{i<j}\frac{1}{r_{ij}} - u^{(\mathrm{insc})}_{c},
$$
$$
u^{(\mathrm{insc})}_{b\mathrm{ex}} = \frac{1}{2}\left(\frac{1}{N_{\mathrm{insc}}}\sum_{i}r_i^2 - u^{(\mathrm{insc})}_{c}\right)
$$
Inspect all $N$ ions and select those with $r_i^2 \le R_{\mathrm{insc}}^2$.
Sum over each unordered pair of selected ions exactly once, using direct distances
without a minimum-image reduction. The geometric mean population
$N_{\mathrm{insc}}$ is not rounded; the actual selected count is reported separately.

# Function `circumscribed`:

We imagine a sphere, circumscribed over the cube (cubic cell).
The sphere' radius is then $R_{\mathrm{circ}}=\sqrt{3}L/2$.
The number of ions in the sphere (average) $N_{\mathrm{circ}}=N \frac{\sqrt{3}\pi}{2}$.
We use periodic boundary conditions (PBC), therefore imagine images of ion and treat them as real ions (that's intentional).
What function should do is for ions inside the circumscribed sphere calculate the same properties, as `inscribed`:
$$
u^{(\mathrm{circ})}_c = \frac{3}{5}N_{\mathrm{circ}}^{2/3},
$$
$$
u^{(\mathrm{circ})}_{p\mathrm{ex}} = \frac{1}{N_{\mathrm{circ}}}\sum_{i<j}\frac{1}{r_{ij}} - u^{(\mathrm{circ})}_{c},
$$
$$
u^{(\mathrm{circ})}_{b\mathrm{ex}} = \frac{1}{2}\left(\frac{1}{N_{\mathrm{circ}}}\sum_{i}r_i^2 - u^{(\mathrm{circ})}_{c}\right)
$$
For a cube, consider translations $(k_xL,k_yL,k_zL)$ with each $k$ in
$\{-1,0,1\}$, including the original cell. Select only points satisfying
$r_i^2 \le R_{\mathrm{circ}}^2$. Face and edge neighbors contribute; corner
neighbors can contribute points exactly on the boundary. Treat every selected
image as a distinct point, including different images of the same original ion.
Sum over unordered pairs exactly once using direct distances, without a
minimum-image reduction. Normalize by the unrounded geometric mean population,
not the actual selected count.

For an orthorhombic box, the implementation uses
$R_{\mathrm{insc}}=\min(L_x,L_y,L_z)/2$,
$R_{\mathrm{circ}}=\sqrt{L_x^2+L_y^2+L_z^2}/2$, and mean population
$N(4\pi R^3/3)/(L_xL_yL_z)$. Image bounds expand as necessary along each axis.

An empty selected sphere with $N>0$ still has nonzero $u_c$ and returns
$u_{p\mathrm{ex}}=-u_c$ and $u_{b\mathrm{ex}}=-u_c/2$. The formulas assume a
nonempty simulation ($N>0$).

# Virial and pressure contribution

For either sphere, write $N_s$ for its unrounded mean population and
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

With `--gamma Gamma`, use a constant reduced temperature $T=1/\Gamma$ for
every frame. Total pressure and the ion compressibility factor are then
$$
p_s=n(T+w_s)=\frac{n}{\Gamma}+p_{\mathrm{ex},s},
\qquad
Z_{\mathrm{ion},s}=\frac{p_s V_s}{N_s T}
=\frac{p_s V}{NT}=1+\Gamma w_s.
$$
Here $N_s$ is the same mean population used for energies and virials, not the
instantaneous selected count. This estimates the kinetic contribution from
the supplied constant temperature; it does not measure temperature fluctuations
from the trajectory.

# Running and output

Use `-i` for the inscribed sphere and `-c` for the circumscribed sphere;
the options can be combined. `Ni` and `Nc` in `observables.csv` are mean
populations, while `Nri` and `Nrc` are actual selected counts. The circumscribed
outputs are appended as `Nc,Nrc,upexc,ubexc,ucc`.

The additional columns `viriali,pexi,virialc,pexc` report $w_s$ and
$p_{\mathrm{ex},s}$ for the inscribed and circumscribed spheres respectively.
They are calculated automatically with `-i` and `-c`. A disabled sphere has
`nan` in these columns; an empty selected sphere with $N>0$ has zero virial
and zero excess pressure.

The optional `--gamma Gamma` argument must be finite and positive, with a
representable finite positive reciprocal. It appends values in columns
`Gamma,T,pi,Zi,pc,Zc`: the coupling parameter, inferred temperature, and total
pressure and compressibility factor for each sphere. Without `--gamma`, these
columns contain `nan`; pressure and compressibility also contain `nan` for a
disabled sphere. The existing virial and excess-pressure columns are independent
of `--gamma`. Zero virial gives $p_s=nT$ and $Z_{\mathrm{ion},s}=1$.

For example, using ADIOS2 defaults and $\Gamma=10$:

```sh
ocp-framework dumps/dump.bp --config "" -i -c --gamma 10
```

`--rdf dr` produces per-ion pair counts in radial bins of width `dr`, using
minimum-image distances. Distances at or above `rmax` are excluded, including
from a partially filled final bin. Both `--rmax` and `--rmax-frac` are capped
at half the shortest box edge. These counts are not normalized by density and
shell volume to produce $g(r)$.
