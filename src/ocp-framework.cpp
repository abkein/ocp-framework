#include "ocp-framework.hpp"

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
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
  [[nodiscard]] inline double min_image(double dx, double L)
  {
    // nearest integer multiple
    dx -= L * std::nearbyint(dx / L);
    return dx;
  }

  [[nodiscard]] inline double N_in_sphere(const size_t N, const Box& box, const double radius)
  {
    const double sphere_volume = (4.0 / 3.0) * std::numbers::pi * radius * radius * radius;
    return static_cast<double>(N) * sphere_volume / box.V;
  }

  // The sums count actual points; normalization uses the mean population.
  inline void normalize_sphere(const double mean_population, double& pair_sum, double& radius_sum, double& uc)
  {
    if (mean_population == 0.0) {
      uc = pair_sum = radius_sum = 0.0;
      return;
    }
    uc         = (3.0 / 5.0) * std::pow(mean_population, 2.0 / 3.0);
    pair_sum   = pair_sum / mean_population - uc;
    radius_sum = 0.5 * (radius_sum / mean_population - uc);
  }

  [[nodiscard]] SphereVirial sphere_virial(const double mean_population, const double pair_sum, const double radius_sum, const double density,
                                           const double gamma)
  {
    SphereVirial result;
    if (mean_population == 0.0) { return result; }
    // Coulomb pairs contribute +1/r; the background force -r contributes -r^2.
    result.virial          = (pair_sum - radius_sum) / (3.0 * mean_population);
    result.excess_pressure = density * result.virial;
    if (std::isfinite(gamma)) {
      result.pressure        = density * (1.0 / gamma + result.virial);
      result.compressibility = 1.0 + gamma * result.virial;
    }
    return result;
  }
}    // namespace

Application::Application(const Options& options) :
    info(options.info), adios(adios2::ADIOS(options.adiosConfig)), calculate_rdf(options.calculate_rdf),
    calculate_inscribed(options.calculate_inscribed), calculate_circumscribed(options.calculate_circumscribed), gamma(options.gamma)
{
  io = adios.DeclareIO(options.ioName);

  if (!(calculate_rdf || calculate_inscribed || calculate_circumscribed)) { throw std::runtime_error("Need at least one calculation specified."); }

  // Read file
  reader  = io.Open(options.bpFile, adios2::Mode::Read);

  // Output files
  out_csv = std::ofstream(options.csv_out);
  if (!out_csv) { throw std::runtime_error("Cannot open " + options.csv_out); }
  out_csv << std::setprecision(std::numeric_limits<double>::max_digits10);
  out_csv << "timestep,N,V,rmax,dr,nbins,Ni,Nri,upexi,ubexi,uci,Nc,Nrc,upexc,ubexc,ucc,viriali,pexi,virialc,pexc,Gamma,T,pi,Zi,pc,Zc\n";

  if (calculate_rdf) {
    out_rdf = std::ofstream(options.rdf_out);
    if (!out_rdf) { throw std::runtime_error("Cannot open " + options.rdf_out); }
    out_rdf << std::setprecision(std::numeric_limits<double>::max_digits10);
  }
}

void Application::prepare_data()
{
  // BP5 exposes metadata only after BeginStep(). Refresh handles per step.
  // Attributes (constant across steps)
  if (auto const a = io.InquireAttribute<int32_t>("scaled"); a) { scaled = a.Data()[0]; }

  if (auto const a = io.InquireAttribute<int32_t>("triclinic"); a) { triclinic = a.Data()[0]; }
  if (triclinic != 0) { throw std::runtime_error("triclinic=1 is not supported by this program (needs tilt factors)."); }

  // Inquire variables
  vAtoms = io.InquireVariable<double>("atoms");
  if (!vAtoms) { throw std::runtime_error("Variable 'atoms' not found (expected double)."); }

  vNt = io.InquireVariable<uint64_t>("ntimestep");
  if (!vNt) { throw std::runtime_error("Variable 'ntimestep' not found."); }

  vNatoms = io.InquireVariable<uint64_t>("natoms");
  if (!vNatoms) { throw std::runtime_error("Variable 'natoms' not found."); }

  vXlo = io.InquireVariable<double>("boxxlo");
  vXhi = io.InquireVariable<double>("boxxhi");
  vYlo = io.InquireVariable<double>("boxylo");
  vYhi = io.InquireVariable<double>("boxyhi");
  vZlo = io.InquireVariable<double>("boxzlo");
  vZhi = io.InquireVariable<double>("boxzhi");
  if (!vXlo || !vXhi || !vYlo || !vYhi || !vZlo || !vZhi) { throw std::runtime_error("Box bound variables not found (box{x,y,z}{lo,hi})."); }

  reader.Get(vNt, timestep);
  uint64_t natoms{};
  reader.Get(vNatoms, natoms);
  N = static_cast<size_t>(natoms);
  double xlo{}, xhi{}, ylo{}, yhi{}, zlo{}, zhi{};
  reader.Get(vXlo, xlo);
  reader.Get(vXhi, xhi);
  reader.Get(vYlo, ylo);
  reader.Get(vYhi, yhi);
  reader.Get(vZlo, zlo);
  reader.Get(vZhi, zhi);

  box = Box(xlo, xhi, ylo, yhi, zlo, zhi);

  // atoms shape is per-step. We rely on natoms (scalar) and ncolumns=5.
  // Select [0..N) x [0..NUM_COLS)
  vAtoms.SetSelection({{0, 0}, {N, COLS::NUM_COLS}});
  atoms.resize(N * COLS::NUM_COLS);

  reader.Get(vAtoms, atoms.data());
  reader.PerformGets();
  reader.EndStep();

  x.reserve(N);
  x.resize(N);
  y.reserve(N);
  y.resize(N);
  z.reserve(N);
  z.resize(N);

  // columns: [id, type, x, y, z]
  if (scaled != 0) {
    // cx,cy,cz are fractional (0..1), we need real
    for (size_t i = 0; i < N; ++i) {
      x[i] = box.xlo + box.Lx * atoms[COLS::NUM_COLS * i + COLS::X];
      y[i] = box.ylo + box.Ly * atoms[COLS::NUM_COLS * i + COLS::Y];
      z[i] = box.zlo + box.Lz * atoms[COLS::NUM_COLS * i + COLS::Z];
    }
  } else {
    // Coordinates are already wrapped into the box by the input contract.
    for (size_t i = 0; i < N; ++i) {
      x[i] = atoms[COLS::NUM_COLS * i + COLS::X];
      y[i] = atoms[COLS::NUM_COLS * i + COLS::Y];
      z[i] = atoms[COLS::NUM_COLS * i + COLS::Z];
    }
  }
}

void Application::inscribed()
{
  const double radius = 0.5 * box.L_min;
  const double Rsq    = radius * radius;
  Ni                  = N_in_sphere(N, box, radius);

  upexi               = 0;
  ubexi               = 0;

  idlist.reserve(N);
  idlist.clear();
  for (size_t i = 0; i < N; ++i) {
    if (sr[i] <= Rsq) {
      ubexi += sr[i];
      idlist.emplace_back(i);
    }
  }
  Nri = idlist.size();

  for (size_t i = 0; i < Nri; ++i) {
    const auto ii   = idlist[i];
    const double xi = x[ii];
    const double yi = y[ii];
    const double zi = z[ii];
    for (size_t j = i + 1; j < Nri; ++j) {
      const auto jj    = idlist[j];
      const double dx  = xi - x[jj];
      const double dy  = yi - y[jj];
      const double dz  = zi - z[jj];
      const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
      upexi += 1.0 / rij;
    }
  }

  viriali = sphere_virial(Ni, upexi, ubexi, static_cast<double>(N) / box.V, gamma);
  normalize_sphere(Ni, upexi, ubexi, uci);
}

void Application::circumscribed()
{
  const double Rsq    = 0.25 * (box.Lx * box.Lx + box.Ly * box.Ly + box.Lz * box.Lz);
  const double radius = std::sqrt(Rsq);
  Nc                  = N_in_sphere(N, box, radius);
  upexc               = 0.0;
  ubexc               = 0.0;

  struct Point {
    double x, y, z;
  };
  std::vector<Point> points;
  points.reserve(N);

  // An image at shift k can intersect the sphere only if
  // (|k| - 1/2) L <= radius. For cubes this gives shifts -1, 0, 1.
  const int nx = static_cast<int>(std::floor(radius / box.Lx + 0.5));
  const int ny = static_cast<int>(std::floor(radius / box.Ly + 0.5));
  const int nz = static_cast<int>(std::floor(radius / box.Lz + 0.5));
  for (int ix = -nx; ix <= nx; ++ix) {
    for (int iy = -ny; iy <= ny; ++iy) {
      for (int iz = -nz; iz <= nz; ++iz) {
        for (size_t i = 0; i < N; ++i) {
          const Point p{.x = x[i] + ix * box.Lx, .y = y[i] + iy * box.Ly, .z = z[i] + iz * box.Lz};
          const double rsq = p.x * p.x + p.y * p.y + p.z * p.z;
          if (rsq <= Rsq) {
            points.push_back(p);
            ubexc += rsq;
          }
        }
      }
    }
  }
  Nrc = points.size();

  // Every selected image is a distinct point, including images of the same
  // original ion. Use direct distances, without a minimum-image reduction.
  for (size_t i = 0; i < Nrc; ++i) {
    for (size_t j = i + 1; j < Nrc; ++j) {
      const double dx = points[i].x - points[j].x;
      const double dy = points[i].y - points[j].y;
      const double dz = points[i].z - points[j].z;
      upexc += 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
  virialc = sphere_virial(Nc, upexc, ubexc, static_cast<double>(N) / box.V, gamma);
  normalize_sphere(Nc, upexc, ubexc, ucc);
}

void Application::calculate()
{
  if (calculate_rdf) { compute_rdf_celllist(); }

  sr.resize(N);
  for (size_t i = 0; i < N; ++i) {
    x[i] -= box.cx;
    y[i] -= box.cy;
    z[i] -= box.cz;
    sr[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
  }

  if (calculate_inscribed) { inscribed(); }
  if (calculate_circumscribed) { circumscribed(); }
}

void Application::write()
{
  out_csv << timestep << "," << N << "," << box.V << "," << rmax << "," << dr << "," << nbins << "," << Ni << "," << Nri << "," << upexi << ","
          << ubexi << "," << uci << "," << Nc << "," << Nrc << "," << upexc << "," << ubexc << "," << ucc << "," << viriali.virial << ","
          << viriali.excess_pressure << "," << virialc.virial << "," << virialc.excess_pressure << "," << gamma << "," << 1.0 / gamma << ","
          << viriali.pressure << "," << viriali.compressibility << "," << virialc.pressure << "," << virialc.compressibility << "\n";
  out_csv.flush();
  if (calculate_rdf) {
    out_rdf << timestep;
    for (size_t i = 0; i < nbins; ++i) { out_rdf << "," << counts[i]; }
    out_rdf << "\n";
    out_rdf.flush();
  }
}

void Application::run()
{
  while (true) {
    const auto status = reader.BeginStep();
    if (status != adios2::StepStatus::OK) { break; }

    prepare_data();

    calculate();
    write();
  }
}

void Application::compute_rdf_celllist()
{
  // enforce safe maximum for minimum-image RDF
  const double requested_rmax = info.rmax <= 0.0 ? info.rmaxFrac * box.L_min : info.rmax;
  rmax                        = std::min(requested_rmax, 0.5 * box.L_min);
  dr                          = info.dr;
  nbins                       = std::max<size_t>(1, static_cast<size_t>(std::ceil(rmax / dr)));
  counts.resize(nbins);
  counts.assign(nbins, 0);
  if (N < 2) { return; }

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = i + 1; j < N; ++j) {
      const double dx       = min_image(x[i] - x[j], box.Lx);
      const double dy       = min_image(y[i] - y[j], box.Ly);
      const double dz       = min_image(z[i] - z[j], box.Lz);
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (distance >= rmax) { continue; }
      const auto ibin = static_cast<size_t>(distance / dr);    // floor
      counts[ibin < nbins ? ibin : nbins - 1] += 2;
    }
  }
  for (size_t i = 0; i < nbins; ++i) { counts[i] /= static_cast<double>(N); }
}

Application::~Application()
{
  out_csv.close();
  if (calculate_rdf) { out_rdf.close(); }
  reader.Close();
}

int main(int argc, char** argv)
{
  try {
    Application app(Options(argc, argv));
    app.run();

  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 1;
  }
}
