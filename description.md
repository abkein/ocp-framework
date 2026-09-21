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

# Running and output

Use `-i` for the inscribed sphere and `-c` for the circumscribed sphere;
the options can be combined. `Ni` and `Nc` in `observables.csv` are mean
populations, while `Nri` and `Nrc` are actual selected counts. The circumscribed
outputs are appended as `Nc,Nrc,upexc,ubexc,ucc`.

`--rdf dr` produces per-ion pair counts in radial bins of width `dr`, using
minimum-image distances. Distances at or above `rmax` are excluded, including
from a partially filled final bin. Both `--rmax` and `--rmax-frac` are capped
at half the shortest box edge. These counts are not normalized by density and
shell volume to produce $g(r)$.
