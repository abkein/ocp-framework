#include <adios2/common/ADIOSTypes.h>
#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Engine.h>
#include <adios2/cxx/IO.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct Options {
  std::string bpFile;
  std::string adiosConfig = "adios2_config.xml";
  std::string ioName      = "read_dump";

  int nbins               = 200;
  double rmax             = -1.0;    // if < 0: use rmaxFrac * min(L)
  double rmaxFrac         = 0.5;     // default: min(L)/2
  double cellSize         = -1.0;    // if < 0: use rmax

  bool pbcCom             = true;    // circular mean COM

  std::string obsOut      = "observables.csv";
  std::string rdfOut      = "rdf.csv";
};

struct Box {
  double xlo, xhi, ylo, yhi, zlo, zhi;
  [[nodiscard]] double Lx() const { return xhi - xlo; }
  [[nodiscard]] double Ly() const { return yhi - ylo; }
  [[nodiscard]] double Lz() const { return zhi - zlo; }
  [[nodiscard]] double V() const { return Lx() * Ly() * Lz(); }
};

namespace {
bool starts_with(const std::string& s, const std::string& p)
{
  return s.size() >= p.size() && std::equal(p.begin(), p.end(), s.begin(), s.end());
}

Options parse_args(int argc, char** argv)
{
  Options opt;
  if (argc < 2) {
    throw std::runtime_error("Usage: analyze_lammps_adios <file.bp> [--config adios2_config.xml] [--io read_dump]\n"
                             "                            [--bins 200] [--rmax 0.0] [--rmax-frac 0.5]\n"
                             "                            [--cell-size 0.0] [--no-pbc-com]\n"
                             "                            [--obs-out observables.csv] [--rdf-out rdf.csv]\n");
  }

  opt.bpFile      = argv[1];

  auto need_value = [&](int& i) -> std::string {
    if (i + 1 >= argc) { throw std::runtime_error(std::string("Missing value after ") + argv[i]); }
    return argv[++i];
  };

  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--config") {
      opt.adiosConfig = need_value(i);
    } else if (a == "--io") {
      opt.ioName = need_value(i);
    } else if (a == "--bins") {
      opt.nbins = std::stoi(need_value(i));
    } else if (a == "--rmax") {
      opt.rmax = std::stod(need_value(i));
    } else if (a == "--rmax-frac") {
      opt.rmaxFrac = std::stod(need_value(i));
    } else if (a == "--cell-size") {
      opt.cellSize = std::stod(need_value(i));
    } else if (a == "--no-pbc-com") {
      opt.pbcCom = false;
    } else if (a == "--obs-out") {
      opt.obsOut = need_value(i);
    } else if (a == "--rdf-out") {
      opt.rdfOut = need_value(i);
    } else {
      throw std::runtime_error("Unknown argument: " + a);
    }
  }

  if (opt.nbins <= 0) { throw std::runtime_error("--bins must be > 0"); }
  if (opt.rmaxFrac <= 0.0) { throw std::runtime_error("--rmax-frac must be > 0"); }
  return opt;
}

inline double wrap01(double s)
{
  s -= std::floor(s);
  if (s >= 1.0) { s = 0.0; }
  if (s < 0.0) { s += 1.0; }
  return s;
}

// periodic-aware mean along one dimension using circular mean on [0,1)
double circular_mean_frac(const std::vector<double>& sfrac)
{
  constexpr double two_pi = 6.283185307179586476925286766559;
  double csum = 0.0, ssum = 0.0;
  for (const double s : sfrac) {
    const double ang = two_pi * wrap01(s);
    csum += std::cos(ang);
    ssum += std::sin(ang);
  }
  const double mean_ang = std::atan2(ssum, csum);
  double mean_s         = mean_ang / two_pi;
  mean_s                = wrap01(mean_s);
  return mean_s;
}

inline double min_image(double dx, double L)
{
  // nearest integer multiple
  dx -= L * std::nearbyint(dx / L);
  return dx;
}

void compute_rdf_celllist(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const Box& box, double rmax,
                          int nbins, double cellSize, std::vector<std::uint64_t>& counts_out)
{
  const std::size_t N = x.size();
  counts_out.assign((std::size_t) nbins, 0);

  if (N < 2) { return; }

  const double Lx = box.Lx(), Ly = box.Ly(), Lz = box.Lz();
  if (rmax <= 0.0) { return; }

  // cellSize choice
  if (cellSize <= 0.0) { cellSize = rmax; }

  int nx           = std::max(1, (int) std::floor(Lx / cellSize));
  int ny           = std::max(1, (int) std::floor(Ly / cellSize));
  const int nz     = std::max(1, (int) std::floor(Lz / cellSize));

  const int ncells = nx * ny * nz;
  std::vector<std::vector<int>> cells((std::size_t) ncells);
  cells.shrink_to_fit();

  auto cell_index = [&](int ix, int iy, int iz) -> int {
    return (iz * ny + iy) * nx + ix;
  };

  // build cell lists
  for (std::size_t i = 0; i < N; ++i) {
    double fx = (x[i] - box.xlo) / Lx;
    double fy = (y[i] - box.ylo) / Ly;
    double fz = (z[i] - box.zlo) / Lz;
    fx        = wrap01(fx);
    fy        = wrap01(fy);
    fz        = wrap01(fz);

    int ix    = (int) std::floor(fx * nx);
    int iy    = (int) std::floor(fy * ny);
    int iz    = (int) std::floor(fz * nz);

    if (ix == nx) { ix = nx - 1; }
    if (iy == ny) { iy = ny - 1; }
    if (iz == nz) { iz = nz - 1; }

    cells[(std::size_t) cell_index(ix, iy, iz)].push_back((int) i);
  }

  const double dr = rmax / (double) nbins;

  for (int c = 0; c < ncells; ++c) {
    const int iz  = c / (nx * ny);
    const int rem = c - iz * (nx * ny);
    const int iy  = rem / nx;
    const int ix  = rem - iy * nx;

    const auto& A = cells[(std::size_t) c];
    if (A.empty()) { continue; }

    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (dz < 0) { continue; }
          if (dz == 0 && dy < 0) { continue; }
          if (dz == 0 && dy == 0 && dx < 0) { continue; }

          int jx = ix + dx;
          if (jx < 0) { jx += nx; }
          if (jx >= nx) { jx -= nx; }
          int jy = iy + dy;
          if (jy < 0) { jy += ny; }
          if (jy >= ny) { jy -= ny; }
          int jz = iz + dz;
          if (jz < 0) { jz += nz; }
          if (jz >= nz) { jz -= nz; }

          const int cB  = cell_index(jx, jy, jz);
          const auto& B = cells[(std::size_t) cB];
          if (B.empty()) { continue; }

          if (cB == c) {
            for (std::size_t ia = 0; ia < A.size(); ++ia) {
              const int i = A[ia];
              for (std::size_t ja = ia + 1; ja < A.size(); ++ja) {
                const int j      = A[ja];
                const double dxr = min_image(x[j] - x[i], Lx);
                const double dyr = min_image(y[j] - y[i], Ly);
                const double dzr = min_image(z[j] - z[i], Lz);
                const double r2  = dxr * dxr + dyr * dyr + dzr * dzr;
                if (r2 < rmax * rmax) {
                  const double r = std::sqrt(r2);
                  const int bin  = (int) std::floor(r / dr);
                  if (bin >= 0 && bin < nbins) { counts_out[(std::size_t) bin] += 1; }
                }
              }
            }
          } else {
            for (const int i : A) {
              for (const int j : B) {
                const double dxr = min_image(x[j] - x[i], Lx);
                const double dyr = min_image(y[j] - y[i], Ly);
                const double dzr = min_image(z[j] - z[i], Lz);
                const double r2  = dxr * dxr + dyr * dyr + dzr * dzr;
                if (r2 < rmax * rmax) {
                  const double r = std::sqrt(r2);
                  const int bin  = (int) std::floor(r / dr);
                  if (bin >= 0 && bin < nbins) { counts_out[(std::size_t) bin] += 1; }
                }
              }
            }
          }
        }
      }
    }
  }
}
}    // namespace

int main(int argc, char** argv)
{
  try {
    const Options opt = parse_args(argc, argv);

    adios2::ADIOS adios(opt.adiosConfig);
    adios2::IO io         = adios.DeclareIO(opt.ioName);

    // Read file
    adios2::Engine reader = io.Open(opt.bpFile, adios2::Mode::Read);

    // Attributes (constant across steps)
    int32_t scaled        = 0;
    if (auto a = io.InquireAttribute<int32_t>("scaled"); a) { scaled = a.Data()[0]; }

    int32_t triclinic = 0;
    if (auto a = io.InquireAttribute<int32_t>("triclinic"); a) { triclinic = a.Data()[0]; }

    if (triclinic != 0) { throw std::runtime_error("triclinic=1 is not supported by this program (needs tilt factors)."); }

    // Output files
    std::ofstream obs(opt.obsOut);
    if (!obs) { throw std::runtime_error("Cannot open " + opt.obsOut); }
    obs << "timestep,N,density,com_x,com_y,com_z\n";

    std::ofstream rdf(opt.rdfOut);
    if (!rdf) { throw std::runtime_error("Cannot open " + opt.rdfOut); }
    rdf << "timestep,r,g_r\n";

    // Inquire variables (types based on your bp info)
    auto vAtoms = io.InquireVariable<double>("atoms");
    if (!vAtoms) { throw std::runtime_error("Variable 'atoms' not found (expected double)."); }

    auto vNt = io.InquireVariable<uint64_t>("ntimestep");
    if (!vNt) { throw std::runtime_error("Variable 'ntimestep' not found."); }

    auto vNatoms = io.InquireVariable<uint64_t>("natoms");
    if (!vNatoms) { throw std::runtime_error("Variable 'natoms' not found."); }

    auto vXlo = io.InquireVariable<double>("boxxlo");
    auto vXhi = io.InquireVariable<double>("boxxhi");
    auto vYlo = io.InquireVariable<double>("boxylo");
    auto vYhi = io.InquireVariable<double>("boxyhi");
    auto vZlo = io.InquireVariable<double>("boxzlo");
    auto vZhi = io.InquireVariable<double>("boxzhi");
    if (!vXlo || !vXhi || !vYlo || !vYhi || !vZlo || !vZhi) { throw std::runtime_error("Box bound variables not found (box{x,y,z}{lo,hi})."); }

    std::vector<double> atoms;
    std::vector<double> x, y, z;
    std::vector<double> sfx, sfy, sfz;    // fractional for circular mean COM
    std::vector<std::uint64_t> rdf_counts;

    constexpr double pi = 3.1415926535897932384626433832795;

    while (true) {
      const auto status = reader.BeginStep();
      if (status != adios2::StepStatus::OK) { break; }

      uint64_t timestep = 0, natoms = 0;
      Box box{};

      reader.Get(vNt, timestep);
      reader.Get(vNatoms, natoms);
      reader.Get(vXlo, box.xlo);
      reader.Get(vXhi, box.xhi);
      reader.Get(vYlo, box.ylo);
      reader.Get(vYhi, box.yhi);
      reader.Get(vZlo, box.zlo);
      reader.Get(vZhi, box.zhi);

      // atoms shape is per-step. We rely on natoms (scalar) and ncolumns=5.
      const auto N = static_cast<std::size_t>(natoms);

      // Select [0..N) x [0..5)
      vAtoms.SetSelection({{0, 0}, {N, 5}});
      atoms.resize(N * 5);

      reader.Get(vAtoms, atoms.data());
      reader.PerformGets();
      reader.EndStep();

      const double Lx = box.Lx(), Ly = box.Ly(), Lz = box.Lz();
      const double V = box.V();

      x.resize(N);
      y.resize(N);
      z.resize(N);
      sfx.resize(N);
      sfy.resize(N);
      sfz.resize(N);

      // columns: [id, type, x?, y?, z?] (you said xs,ys,zs in columns 2..4)
      for (std::size_t i = 0; i < N; ++i) {
        const double cx = atoms[i * 5 + 2];
        const double cy = atoms[i * 5 + 3];
        const double cz = atoms[i * 5 + 4];

        if (scaled != 0) {
          // cx,cy,cz are fractional (0..1)
          const double fx = wrap01(cx);
          const double fy = wrap01(cy);
          const double fz = wrap01(cz);

          sfx[i]          = fx;
          sfy[i]          = fy;
          sfz[i]          = fz;
          x[i]            = box.xlo + fx * Lx;
          y[i]            = box.ylo + fy * Ly;
          z[i]            = box.zlo + fz * Lz;
        } else {
          // cx,cy,cz are real coordinates; wrap into box for consistency
          double fx = (cx - box.xlo) / Lx;
          double fy = (cy - box.ylo) / Ly;
          double fz = (cz - box.zlo) / Lz;

          fx        = wrap01(fx);
          fy        = wrap01(fy);
          fz        = wrap01(fz);

          sfx[i]    = fx;
          sfy[i]    = fy;
          sfz[i]    = fz;
          x[i]      = box.xlo + fx * Lx;
          y[i]      = box.ylo + fy * Ly;
          z[i]      = box.zlo + fz * Lz;
        }
      }

      // Density
      const double density = (V > 0.0) ? ((double) N / V) : 0.0;

      // Center of mass (unit masses)
      double comx = 0.0, comy = 0.0, comz = 0.0;
      if (opt.pbcCom) {
        const double mx = circular_mean_frac(sfx);
        const double my = circular_mean_frac(sfy);
        const double mz = circular_mean_frac(sfz);
        comx            = box.xlo + mx * Lx;
        comy            = box.ylo + my * Ly;
        comz            = box.zlo + mz * Lz;
      } else {
        for (std::size_t i = 0; i < N; ++i) {
          comx += x[i];
          comy += y[i];
          comz += z[i];
        }
        comx /= (double) N;
        comy /= (double) N;
        comz /= (double) N;
      }

      obs << timestep << "," << N << "," << density << "," << comx << "," << comy << "," << comz << "\n";

      // RDF parameters for this step
      double rmax = opt.rmax;
      if (rmax <= 0.0) { rmax = opt.rmaxFrac * std::min({Lx, Ly, Lz}); }

      // enforce safe maximum for minimum-image RDF
      const double rmax_limit = 0.5 * std::min({Lx, Ly, Lz});
      rmax                    = std::min(rmax, rmax_limit);

      const double cellSize   = (opt.cellSize > 0.0) ? opt.cellSize : rmax;

      compute_rdf_celllist(x, y, z, box, rmax, opt.nbins, cellSize, rdf_counts);

      const double dr         = rmax / (double) opt.nbins;
      const double norm_denom = (double) N * (double) (N - 1);

      for (int b = 0; b < opt.nbins; ++b) {
        const double r1       = b * dr;
        const double r2       = (b + 1) * dr;
        const double rmid     = 0.5 * (r1 + r2);

        const double shellVol = (4.0 / 3.0) * pi * (r2 * r2 * r2 - r1 * r1 * r1);
        double g              = 0.0;
        if (shellVol > 0.0 && V > 0.0 && N > 1) {
          // g = counts / expected_ideal,
          // expected_ideal = N(N-1)/2 * shellVol / V
          // => g = 2 V counts / (N(N-1) shellVol)
          g = (2.0 * V * (double) rdf_counts[(std::size_t) b]) / (norm_denom * shellVol);
        }

        rdf << timestep << "," << rmid << "," << g << "\n";
      }
    }

    reader.Close();
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}