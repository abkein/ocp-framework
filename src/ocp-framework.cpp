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
#include <iostream>
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

  [[nodiscard]] inline size_t N_in_sphere(const size_t N, const Box& box, const double radius)
  {
    const double sphere_volume = (4.0 / 3.0) * PI * radius * radius * radius;
    const double N_in          = static_cast<double>(N) * sphere_volume / box.V;
    return static_cast<size_t>(std::llround(N_in));
  }
}    // namespace

Application::Application(const Options& opt)
{
  calculate_rdf           = opt.calculate_rdf;
  calculate_inscribed     = opt.calculate_inscribed;
  if (!(calculate_rdf || calculate_inscribed)) { throw std::runtime_error("Need at least one calculation specified."); }

  info   = opt.info;
  adios  = adios2::ADIOS(opt.adiosConfig);
  io     = adios.DeclareIO(opt.ioName);

  // Read file
  reader = io.Open(opt.bpFile, adios2::Mode::Read);

  // Attributes (constant across steps)
  if (auto a = io.InquireAttribute<int32_t>("scaled"); a) { scaled = a.Data()[0]; }

  if (auto a = io.InquireAttribute<int32_t>("triclinic"); a) { triclinic = a.Data()[0]; }
  if (triclinic != 0) { throw std::runtime_error("triclinic=1 is not supported by this program (needs tilt factors)."); }

  // Output files
  obs = std::ofstream(opt.obsOut);
  if (!obs) { throw std::runtime_error("Cannot open " + opt.obsOut); }
  obs << "timestep,N,V,rmax,dr,nbins,Ni,Nri,upexi,ubexi,uci\n";

  if (calculate_rdf) {
    rdf = std::ofstream(opt.rdfOut);
    if (!rdf) { throw std::runtime_error("Cannot open " + opt.rdfOut); }
  }

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
}

void Application::prepare_data()
{
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
    // cx,cy,cz are real coordinates; wrap into box for consistency
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
      // if (rij == 0.0) { continue; }
      upexi += 1.0 / rij;
    }
  }

  if (Nri == 0) {
    uci   = 0.0;
    upexi = 0.0;
    ubexi = 0.0;
    return;
  }

  uci = (3.0 / 5.0) * std::pow(Nri, 2.0 / 3.0);
  upexi /= static_cast<double>(Nri);
  upexi -= uci;

  ubexi /= static_cast<double>(Nri);
  ubexi -= uci;
  ubexi *= 0.5;
}

void Application::calculate()
{
  if (calculate_rdf) { compute_rdf_celllist(); }

  // for (size_t i = 0; i < N; ++i) { x[i] -= box.cx; }
  // for (size_t i = 0; i < N; ++i) { y[i] -= box.cy; }
  // for (size_t i = 0; i < N; ++i) { z[i] -= box.cz; }
  // for (size_t i = 0; i < N; ++i) { sr[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i]; }

  sr.resize(N);
  for (size_t i = 0; i < N; ++i) {
    x[i] -= box.cx;
    y[i] -= box.cy;
    z[i] -= box.cz;
    sr[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
  }

  if (calculate_inscribed) { inscribed(); }
}

void Application::write()
{
  obs << timestep << "," << N << "," << box.V << "," << rmax << "," << dr << "," << nbins << "," << Ni << "," << Nri << "," << upexi << "," << ubexi
      << "," << uci << "\n";
  obs.flush();
  if (calculate_rdf) {
    rdf << timestep;
    for (size_t i = 0; i < nbins; ++i) { rdf << "," << counts[i]; }
    rdf << "\n";
    rdf.flush();
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
  rmax  = info.rmax <= 0.0 ? info.rmaxFrac * box.L_min : std::min(info.rmax, 0.5 * box.L_min);
  dr    = info.dr;
  nbins = std::max<size_t>(1, static_cast<size_t>(std::ceil(rmax / dr)));
  counts.resize(nbins);
  counts.assign(nbins, 0);
  if (N < 2) { return; }

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = i + 1; j < N; ++j) {
      const double dx       = min_image(x[i] - x[j], box.Lx);
      const double dy       = min_image(y[i] - y[j], box.Ly);
      const double dz       = min_image(z[i] - z[j], box.Lz);
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      const auto ibin       = static_cast<size_t>(distance / dr);    // floor
      counts[ibin < nbins ? ibin : nbins - 1] += 2;
    }
  }
  for (size_t i = 0; i < nbins; ++i) { counts[i] /= static_cast<double>(N); }
}

Application::~Application()
{
  obs.close();
  if (calculate_rdf) { rdf.close(); }
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
