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
  obs << "timestep,N,V,rmax,dr,nbins,Ni,Nri,upexi,ubexi,uci,Nc,Nrc,upexc,ubexc,ucc\n";

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

  x.resize(N);
  y.resize(N);
  z.resize(N);

  // columns: [id, type, x, y, z]
  if (scaled != 0) {
    // cx,cy,cz are fractional (0..1)
    for (size_t i = 0; i < N; ++i) {
      const double cx = atoms[i * COLS::NUM_COLS + COLS::X];
      const double cy = atoms[i * COLS::NUM_COLS + COLS::Y];
      const double cz = atoms[i * COLS::NUM_COLS + COLS::Z];

      x[i]            = box.xlo + cx * box.Lx;
      y[i]            = box.ylo + cy * box.Ly;
      z[i]            = box.zlo + cz * box.Lz;
    }
  } else {
    // cx,cy,cz are real coordinates; wrap into box for consistency
    for (size_t i = 0; i < N; ++i) {
      x[i] = atoms[i * COLS::NUM_COLS + COLS::X];
      y[i] = atoms[i * COLS::NUM_COLS + COLS::Y];
      z[i] = atoms[i * COLS::NUM_COLS + COLS::Z];
    }
  }
  ninside_original.assign(N, true);
}

void Application::inscribed()
{
  Ni               = static_cast<size_t>(PI * static_cast<double>(N) * (1.0 / 6.0));
  const double Rsq = 0.5 * 0.5 * (box.Lx * box.Lx + box.Ly * box.Ly + box.Lz * box.Lz);

  Nri              = 0;
  upexi            = 0;
  ubexi            = 0;

  for (size_t i = 0; i < N; ++i) {
    const double risq = fcx[i] * fcx[i] + fcy[i] * fcy[i] + fcz[i] * fcz[i];
    ubexi += risq;
    if (risq > Rsq) { continue; }
    ++Nri;
    for (size_t j = i + 1; j < N; ++j) {
      const double rjsq = fcx[j] * fcx[j] + fcy[j] * fcy[j] + fcz[j] * fcz[j];
      if (rjsq > Rsq) { continue; }
      const double dx  = x[i] - x[j];
      const double dy  = y[i] - y[j];
      const double dz  = z[i] - z[j];
      const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (rij == 0.0) { continue; }
      upexi += 1.0 / rij;
    }
  }

  if (Nri == 0) {
    uci   = 0.0;
    upexi = 0.0;
    ubexi = 0.0;
    return;
  }

  uci = std::pow(Nri * Nri, (1.0 / 3.0)) * (3.0 / 5.0);
  upexi /= static_cast<double>(Nri);
  ubexi /= static_cast<double>(Nri);
  ubexi *= 0.5;
}

void Application::circumscribed()
{
  Nc               = static_cast<size_t>(PI * static_cast<double>(N) * SQRTofTWO_over_THREE);
  const double Rsq = 0.5 * (box.Lx * box.Lx + box.Ly * box.Ly + box.Lz * box.Lz);

  Nrc              = 0;
  ninside1.assign(N, true);
  cx1.resize(N);
  cy1.resize(N);
  cz1.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cxi      = fcx[i] - box.Lx;
    const double distance = cxi * cxi + fcy[i] * fcy[i] + fcz[i] * fcz[i];
    ninside1[i]           = distance > Rsq;
    if (!ninside1[i]) {
      cx1[i] = cxi;
      cy1[i] = fcy[i];
      cz1[i] = fcz[i];
    }
  }

  ninside2.assign(N, true);
  cx2.resize(N);
  cy2.resize(N);
  cz2.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cxi      = fcx[i] + box.Lx;
    const double distance = cxi * cxi + fcy[i] * fcy[i] + fcz[i] * fcz[i];
    ninside2[i]           = distance > Rsq;
    if (!ninside2[i]) {
      cx2[i] = cxi;
      cy2[i] = fcy[i];
      cz2[i] = fcz[i];
    }
  }

  ninside3.assign(N, true);
  cx3.resize(N);
  cy3.resize(N);
  cz3.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cyi      = fcy[i] - box.Ly;
    const double distance = fcx[i] * fcx[i] + cyi * cyi + fcz[i] * fcz[i];
    ninside3[i]           = distance > Rsq;
    if (!ninside3[i]) {
      cx3[i] = fcx[i];
      cy3[i] = cyi;
      cz3[i] = fcz[i];
    }
  }

  ninside4.assign(N, true);
  cx4.resize(N);
  cy4.resize(N);
  cz4.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cyi      = fcy[i] + box.Ly;
    const double distance = fcx[i] * fcx[i] + cyi * cyi + fcz[i] * fcz[i];
    ninside4[i]           = distance > Rsq;
    if (!ninside4[i]) {
      cx4[i] = fcx[i];
      cy4[i] = cyi;
      cz4[i] = fcz[i];
    }
  }

  ninside5.assign(N, true);
  cx5.resize(N);
  cy5.resize(N);
  cz5.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double czi      = fcz[i] - box.Lz;
    const double distance = fcx[i] * fcx[i] + fcy[i] * fcy[i] + czi * czi;
    ninside5[i]           = distance > Rsq;
    if (!ninside5[i]) {
      cx5[i] = fcx[i];
      cy5[i] = fcy[i];
      cz5[i] = czi;
    }
  }

  ninside6.assign(N, true);
  cx6.resize(N);
  cy6.resize(N);
  cz6.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double czi      = fcz[i] + box.Lz;
    const double distance = fcx[i] * fcx[i] + fcy[i] * fcy[i] + czi * czi;
    ninside6[i]           = distance > Rsq;
    if (!ninside6[i]) {
      cx6[i] = fcx[i];
      cy6[i] = fcy[i];
      cz6[i] = czi;
    }
  }

  upexc = 0;
  ubexc = 0;

  for (int imgi = 0; imgi < 7; ++imgi) {
    const std::vector<bool>& i_ninside = *(ninside[imgi]);
    const std::vector<double>& ximgi   = *(vcs[imgi][0]);
    const std::vector<double>& yimgi   = *(vcs[imgi][1]);
    const std::vector<double>& zimgi   = *(vcs[imgi][2]);
    for (size_t i = 0; i < N; ++i) {
      if (i_ninside[i]) { continue; }
      ++Nrc;
      const double xi   = ximgi[i];
      const double yi   = yimgi[i];
      const double zi   = zimgi[i];
      const double risq = xi * xi + yi * yi + zi * zi;
      ubexc += risq;
      for (int imgj = 0; imgj < 7; ++imgj) {
        const std::vector<bool>& j_ninside = *(ninside[imgj]);
        const std::vector<double>& ximgj   = *(vcs[imgj][0]);
        const std::vector<double>& yimgj   = *(vcs[imgj][1]);
        const std::vector<double>& zimgj   = *(vcs[imgj][2]);
        for (size_t j = i; j < N; ++j) {
          if (j_ninside[j]) { continue; }
          const double dx  = xi - ximgj[j];
          const double dy  = yi - yimgj[j];
          const double dz  = zi - zimgj[j];
          const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
          if (rij == 0.0) { continue; }
          upexc += 1.0 / rij;
        }
      }
    }
  }

  if (Nrc == 0) {
    ucc   = 0.0;
    upexc = 0.0;
    ubexc = 0.0;
    return;
  }

  ucc = std::pow(Nrc * Nrc, (1.0 / 3.0)) * (3.0 / 5.0);
  upexc /= static_cast<double>(Nrc);
  ubexc /= static_cast<double>(Nrc);
  ubexc *= 0.5;
}

void Application::calculate()
{
  fcx.resize(N);
  for (size_t i = 0; i < N; ++i) { fcx[i] = x[i] - box.cx; }
  fcy.resize(N);
  for (size_t i = 0; i < N; ++i) { fcy[i] = y[i] - box.cy; }
  fcz.resize(N);
  for (size_t i = 0; i < N; ++i) { fcz[i] = z[i] - box.cz; }

  rmax = info.rmax;
  if (rmax <= 0.0) {
    rmax = info.rmaxFrac * std::sqrt(box.Lx * box.Lx + box.Ly * box.Ly + box.Lz * box.Lz) * 0.5;
  }    // 0.5 b/c diameter to radius conversion

  // enforce safe maximum for minimum-image RDF
  rmax  = std::min(rmax, 0.5 * std::min({box.Lx, box.Ly, box.Lz}));

  if (rmax <= 0.0) { throw std::runtime_error("Computed rmax must be > 0."); }
  dr = info.dr;
  if (dr <= 0.0) { dr = rmax; }
  nbins = std::max<size_t>(1, static_cast<size_t>(std::ceil(rmax / dr)));
  dr    = rmax / static_cast<double>(nbins);

  compute_rdf_celllist();
  inscribed();
  circumscribed();
}

void Application::write()
{
  obs << timestep << "," << N << "," << box.V << "," << rmax << "," << dr << "," << nbins << "," << Ni << "," << Nri << "," << upexi << "," << ubexi
      << "," << uci << "," << Nc << "," << Nrc << "," << upexc << "," << ubexc << "," << ucc << "\n";
  rdf << timestep;
  for (size_t i = 0; i < nbins; ++i) { rdf << "," << counts[i]; }
  rdf << "\n";
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
  counts.resize(nbins);
  counts.assign(nbins, 0);
  if (N < 2) { return; }
  if (dr <= 0.0) { return; }

  for (size_t i = 0; i < N - 1; ++i) {
    for (size_t j = i + 1; j < N; ++j) {
      const double dx       = min_image(x[i] - x[j], box.Lx);
      const double dy       = min_image(y[i] - y[j], box.Ly);
      const double dz       = min_image(z[i] - z[j], box.Lz);
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
