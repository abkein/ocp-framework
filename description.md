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
u^{(\mathrm{insc})}_{p\mathrm{ex}} = \frac{1}{N_{\mathrm{insc}}}\sum_{i\ne j}\frac{1}{r_{ij}} - u^{(\mathrm{insc})}_{c},
$$
$$
u^{(\mathrm{insc})}_{b\mathrm{ex}} = \frac{1}{2}\left(\frac{1}{N_{\mathrm{insc}}}\sum_{i}r_i^2 - u^{(\mathrm{insc})}_{c}\right)
$$
So we loop over $i$-th ion ($0\le i<N_{\mathrm{insc}}$), skipping the ion if it is not in the sphere.
Then for each $i$-th ion loop over $j$-th ion ($i+1\le j<N_{\mathrm{insc}}$).
If both ions are inside the sphere we calculate property.

# Function `circumscribed`:

We imagine a sphere, circumscribed over the cube (cubic cell).
The sphere' radius is then $R_{\mathrm{circ}}=L/\sqrt{2}$.
The number of ions in the sphere (average) $N_{\mathrm{circ}}=N \frac{\sqrt{2}\pi}{3}$.
We use periodic boundary conditions (PBC), therefore imagine images of ion and treat them as real ions (that's intentional).
What function should do is for ions inside the circumscribed sphere calculate the same properties, as `inscribed`:
$$
u^{(\mathrm{circ})}_c = \frac{3}{5}N_{\mathrm{circ}}^{2/3},
$$
$$
u^{(\mathrm{circ})}_{p\mathrm{ex}} = \frac{1}{N_{\mathrm{circ}}}\sum_{i\ne j}\frac{1}{r_{ij}} - u^{(\mathrm{circ})}_{c},
$$
$$
u^{(\mathrm{circ})}_{b\mathrm{ex}} = \frac{1}{2}\left(\frac{1}{N_{\mathrm{circ}}}\sum_{i}r_i^2 - u^{(\mathrm{circ})}_{c}\right)
$$
So we first create 6 additional images, therefore we have 7 sets of ions. Loop over each pair and calculate.
