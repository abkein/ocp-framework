#include "analyze_lammps_adios.hpp"

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
}    // namespace

Application::Application(const Options& opt)
{
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
  obs << "timestep,N,V,rmax,dr,nbins\n";

  rdf = std::ofstream(opt.rdfOut);
  if (!rdf) { throw std::runtime_error("Cannot open " + opt.rdfOut); }

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
  reader.Get(vXlo, box.xlo);
  reader.Get(vXhi, box.xhi);
  reader.Get(vYlo, box.ylo);
  reader.Get(vYhi, box.yhi);
  reader.Get(vZlo, box.zlo);
  reader.Get(vZhi, box.zhi);

  // atoms shape is per-step. We rely on natoms (scalar) and ncolumns=5.
  // Select [0..N) x [0..NUM_COLS)
  vAtoms.SetSelection({{0, 0}, {N, COLS::NUM_COLS}});
  atoms.resize(N * COLS::NUM_COLS);

  reader.Get(vAtoms, atoms.data());
  reader.PerformGets();
  reader.EndStep();

  x.resize(N);
  y.resize(N);
  z.resize(N);

  const double Lx = box.Lx();
  const double Ly = box.Ly();
  const double Lz = box.Lz();

  // columns: [id, type, x, y, z]
  if (scaled != 0) {
    // cx,cy,cz are fractional (0..1)
    for (size_t i = 0; i < N; ++i) {
      const double cx = atoms[i * COLS::NUM_COLS + COLS::X];
      const double cy = atoms[i * COLS::NUM_COLS + COLS::Y];
      const double cz = atoms[i * COLS::NUM_COLS + COLS::Z];

      x[i]            = box.xlo + cx * Lx;
      y[i]            = box.ylo + cy * Ly;
      z[i]            = box.zlo + cz * Lz;
    }
  } else {
    // cx,cy,cz are real coordinates; wrap into box for consistency
    for (size_t i = 0; i < N; ++i) {
      x[i] = atoms[i * COLS::NUM_COLS + COLS::X];
      y[i] = atoms[i * COLS::NUM_COLS + COLS::Y];
      z[i] = atoms[i * COLS::NUM_COLS + COLS::Z];
    }
  }
}

void Application::calculate()
{
  const double Lx = box.Lx();
  const double Ly = box.Ly();
  const double Lz = box.Lz();

  double rmax     = info.rmax;
  if (rmax <= 0.0) { rmax = info.rmaxFrac * std::sqrt(Lx * Lx + Ly * Ly + Lz * Lz); }

  // enforce safe maximum for minimum-image RDF
  rmax             = std::min(rmax, 0.5 * std::min({Lx, Ly, Lz}));

  const auto nbins = static_cast<size_t>(rmax / info.dr);
  const double dr  = rmax / static_cast<double>(nbins);

  obs << timestep << "," << N << "," << box.V() << "," << rmax << "," << dr << "," << nbins << "\n";

  compute_rdf_celllist(dr, nbins);

  rdf << timestep << ",";
  for (size_t i = 0; i < nbins; ++i) { rdf << "," << counts[i] << "\n"; }
}

void Application::run()
{
  while (true) {
    const auto status = reader.BeginStep();
    if (status != adios2::StepStatus::OK) { break; }

    prepare_data();

    calculate();
  }
}

void Application::compute_rdf_celllist(const double dr, const size_t nbins)
{
  counts.resize(nbins);
  counts.assign(nbins, 0);
  if (N < 2) { return; }
  if (dr <= 0.0) { return; }
  const double Lx = box.Lx();
  const double Ly = box.Ly();
  const double Lz = box.Lz();

  for (size_t i = 0; i < N - 1; ++i) {
    for (size_t j = i; j < N; ++j) {
      const double dx       = min_image(x[i] - x[j], Lx);
      const double dy       = min_image(y[i] - y[j], Ly);
      const double dz       = min_image(z[i] - z[j], Lz);
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      const auto ibin       = static_cast<size_t>(distance / dr);    // floor
      if (ibin >= nbins) { continue; }
      counts[ibin] += 2;
    }
  }
  for (size_t i = 0; i < nbins; ++i) { counts[i] /= static_cast<double>(N); }
}

Application::~Application()
{
  obs.close();
  rdf.close();
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
