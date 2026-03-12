#include "analyze_lammps_adios.hpp"

#include <adios2/common/ADIOSTypes.h>
#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Engine.h>
#include <adios2/cxx/IO.h>

#include <algorithm>
#include <array>
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

  [[nodiscard]] inline double squared_norm(const double x, const double y, const double z)
  {
    return x * x + y * y + z * z;
  }

  [[nodiscard]] inline size_t N_in_sphere(const size_t N, const Box& box, const double radius)
  {
    const double sphere_volume = (4.0 / 3.0) * PI * radius * radius * radius;
    const double N_in          = static_cast<double>(N) * sphere_volume / box.V;
    return static_cast<size_t>(std::llround(N_in));
  }

  struct SpherePoint {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };
}    // namespace

Application::Application(const Options& opt)
{
  calculate_rdf           = opt.calculate_rdf;
  calculate_inscribed     = opt.calculate_inscribed;
  calculate_circumscribed = opt.calculate_circumscribed;
  if (!(calculate_rdf || calculate_inscribed || calculate_circumscribed)) { throw std::runtime_error("Need at least one calculation specified."); }

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
  // N3 = 3 * N;
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

  const size_t Nreserve = (calculate_circumscribed ? 19 : 1) * N;
  x.reserve(Nreserve);
  x.resize(N);
  y.reserve(Nreserve);
  y.resize(N);
  z.reserve(Nreserve);
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

  for (size_t i = 0; i < N; ++i) {
    if (sr[i] <= Rsq) {
      ubexi += sr[i];
      idlist.emplace_back(i);
    }
  }
  Nri = idlist.size();

  for (size_t i = 0; i < Nri; ++i) {
    for (size_t j = i + 1; j < Nri; ++j) {
      const double dx  = x[idlist[i]] - x[idlist[j]];
      const double dy  = y[idlist[i]] - y[idlist[j]];
      const double dz  = z[idlist[i]] - z[idlist[j]];
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

namespace {
  inline void empl1(const double Rsq, double& acc, const size_t i, std::vector<double>& av, const double a, std::vector<double>& bv,
                    std::vector<double>& sbv, std::vector<double>& cv, std::vector<double>& scv)
  {
    const double rsq = a * a + sbv[i] + scv[i];
    if (rsq <= Rsq) {
      acc += rsq;
      av.emplace_back(a);
      bv.emplace_back(bv[i]);
      cv.emplace_back(cv[i]);
    }
  }

  inline void empl2(const double Rsq, double& acc, const size_t i, std::vector<double>& av, const double a, std::vector<double>& bv, const double b,
                    std::vector<double>& cv, std::vector<double>& scv)
  {
    const double rsq = a * a + b * b + scv[i];
    if (rsq <= Rsq) {
      acc += rsq;
      av.emplace_back(a);
      bv.emplace_back(b);
      cv.emplace_back(cv[i]);
    }
  }
}    // namespace

void Application::circumscribed()
{
  const double Rsq = (1.0 / 4.0) * squared_norm(box.Lx, box.Ly, box.Lz);
  Nc               = N_in_sphere(N, box, std::sqrt(Rsq));

  upexc            = 0;
  ubexc            = 0;

  for (size_t i = 0; i < N; ++i) {
    ubexc += sx[i] + sy[i] + sz[i];

    empl1(Rsq, ubexc, i, x, x[i] - box.Lx, y, sy, z, sz);
    empl1(Rsq, ubexc, i, x, x[i] + box.Lx, y, sy, z, sz);
    empl1(Rsq, ubexc, i, y, y[i] - box.Ly, x, sx, z, sz);
    empl1(Rsq, ubexc, i, y, y[i] + box.Ly, x, sx, z, sz);
    empl1(Rsq, ubexc, i, z, z[i] - box.Lz, x, sx, y, sy);
    empl1(Rsq, ubexc, i, z, z[i] + box.Lz, x, sx, y, sy);

    empl2(Rsq, ubexc, i, x, x[i] - box.Lx, y, y[i] - box.Ly, z, sz);
    empl2(Rsq, ubexc, i, x, x[i] - box.Lx, y, y[i] + box.Ly, z, sz);
    empl2(Rsq, ubexc, i, x, x[i] + box.Lx, y, y[i] - box.Ly, z, sz);
    empl2(Rsq, ubexc, i, x, x[i] + box.Lx, y, y[i] + box.Ly, z, sz);

    empl2(Rsq, ubexc, i, x, x[i] - box.Lx, z, z[i] - box.Lz, y, sy);
    empl2(Rsq, ubexc, i, x, x[i] - box.Lx, z, z[i] + box.Lz, y, sy);
    empl2(Rsq, ubexc, i, x, x[i] + box.Lx, z, z[i] - box.Lz, y, sy);
    empl2(Rsq, ubexc, i, x, x[i] + box.Lx, z, z[i] + box.Lz, y, sy);

    empl2(Rsq, ubexc, i, y, y[i] - box.Ly, z, z[i] - box.Lz, x, sx);
    empl2(Rsq, ubexc, i, y, y[i] - box.Ly, z, z[i] + box.Lz, x, sx);
    empl2(Rsq, ubexc, i, y, y[i] + box.Ly, z, z[i] - box.Lz, x, sx);
    empl2(Rsq, ubexc, i, y, y[i] + box.Ly, z, z[i] + box.Lz, x, sx);
  }

  Nrc = x.size();

  for (size_t i = 0; i < N; ++i) {
    const double xc = x[i];
    const double xm = x[i] - box.Lx;
    const double xp = x[i] + box.Lx;
    const size_t id = N + 18 * i;
    x[id + 0]       = xm;
    x[id + 1]       = xm;
    x[id + 2]       = xm;
    x[id + 3]       = xm;
    x[id + 4]       = xm;
    x[id + 5]       = xp;
    x[id + 6]       = xp;
    x[id + 7]       = xp;
    x[id + 8]       = xp;
    x[id + 9]       = xp;
    x[id + 10]      = xc;
    x[id + 11]      = xc;
    x[id + 12]      = xc;
    x[id + 13]      = xc;
    x[id + 14]      = xc;
    x[id + 15]      = xc;
    x[id + 16]      = xc;
    x[id + 17]      = xc;
  }

  for (size_t i = 0; i < N; ++i) {
    const double yc = y[i];
    const double ym = y[i] - box.Ly;
    const double yp = y[i] + box.Ly;
    const size_t id = N + 18 * i;
    y[id + 0]       = ym;
    y[id + 1]       = yp;
    y[id + 2]       = yc;
    y[id + 3]       = yc;
    y[id + 4]       = yc;
    y[id + 5]       = yc;
    y[id + 6]       = ym;
    y[id + 7]       = yp;
    y[id + 8]       = yc;
    y[id + 9]       = yc;
    y[id + 10]      = yc;
    y[id + 11]      = yc;
    y[id + 12]      = ym;
    y[id + 13]      = ym;
    y[id + 14]      = ym;
    y[id + 15]      = yp;
    y[id + 16]      = yp;
    y[id + 17]      = yp;
  }

  for (size_t i = 0; i < N; ++i) {
    const double zc = z[i];
    const double zm = z[i] - box.Lz;
    const double zp = z[i] + box.Lz;
    const size_t id = N + 18 * i;
    z[id + 0]       = zc;
    z[id + 1]       = zc;
    z[id + 2]       = zc;
    z[id + 3]       = zm;
    z[id + 4]       = zp;
    z[id + 5]       = zc;
    z[id + 6]       = zc;
    z[id + 7]       = zc;
    z[id + 8]       = zm;
    z[id + 9]       = zp;
    z[id + 10]      = zm;
    z[id + 11]      = zp;
    z[id + 12]      = zm;
    z[id + 13]      = zp;
    z[id + 14]      = zc;
    z[id + 15]      = zc;
    z[id + 16]      = zm;
    z[id + 17]      = zp;
  }

  std::vector<SpherePoint> points;
  points.reserve(7 * N);

  ninside_original.assign(N, true);
  for (size_t i = 0; i < N; ++i) {
    const double distance = squared_norm(fcx[i], fcy[i], fcz[i]);
    ninside_original[i]   = distance > Rsq;
    if (!ninside_original[i]) { points.push_back({fcx[i], fcy[i], fcz[i]}); }
  }

  ninside1.assign(N, true);
  cx1.resize(N);
  cy1.resize(N);
  cz1.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cxi      = fcx[i] - box.Lx;
    const double distance = squared_norm(cxi, fcy[i], fcz[i]);
    ninside1[i]           = distance > Rsq;
    if (!ninside1[i]) {
      cx1[i] = cxi;
      cy1[i] = fcy[i];
      cz1[i] = fcz[i];
      points.push_back({cx1[i], cy1[i], cz1[i]});
    }
  }

  ninside2.assign(N, true);
  cx2.resize(N);
  cy2.resize(N);
  cz2.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cxi      = fcx[i] + box.Lx;
    const double distance = squared_norm(cxi, fcy[i], fcz[i]);
    ninside2[i]           = distance > Rsq;
    if (!ninside2[i]) {
      cx2[i] = cxi;
      cy2[i] = fcy[i];
      cz2[i] = fcz[i];
      points.push_back({cx2[i], cy2[i], cz2[i]});
    }
  }

  ninside3.assign(N, true);
  cx3.resize(N);
  cy3.resize(N);
  cz3.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cyi      = fcy[i] - box.Ly;
    const double distance = squared_norm(fcx[i], cyi, fcz[i]);
    ninside3[i]           = distance > Rsq;
    if (!ninside3[i]) {
      cx3[i] = fcx[i];
      cy3[i] = cyi;
      cz3[i] = fcz[i];
      points.push_back({cx3[i], cy3[i], cz3[i]});
    }
  }

  ninside4.assign(N, true);
  cx4.resize(N);
  cy4.resize(N);
  cz4.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double cyi      = fcy[i] + box.Ly;
    const double distance = squared_norm(fcx[i], cyi, fcz[i]);
    ninside4[i]           = distance > Rsq;
    if (!ninside4[i]) {
      cx4[i] = fcx[i];
      cy4[i] = cyi;
      cz4[i] = fcz[i];
      points.push_back({cx4[i], cy4[i], cz4[i]});
    }
  }

  ninside5.assign(N, true);
  cx5.resize(N);
  cy5.resize(N);
  cz5.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double czi      = fcz[i] - box.Lz;
    const double distance = squared_norm(fcx[i], fcy[i], czi);
    ninside5[i]           = distance > Rsq;
    if (!ninside5[i]) {
      cx5[i] = fcx[i];
      cy5[i] = fcy[i];
      cz5[i] = czi;
      points.push_back({cx5[i], cy5[i], cz5[i]});
    }
  }

  ninside6.assign(N, true);
  cx6.resize(N);
  cy6.resize(N);
  cz6.resize(N);
  for (size_t i = 0; i < N; ++i) {
    const double czi      = fcz[i] + box.Lz;
    const double distance = squared_norm(fcx[i], fcy[i], czi);
    ninside6[i]           = distance > Rsq;
    if (!ninside6[i]) {
      cx6[i] = fcx[i];
      cy6[i] = fcy[i];
      cz6[i] = czi;
      points.push_back({cx6[i], cy6[i], cz6[i]});
    }
  }

  Nrc   = points.size();
  upexc = 0;
  ubexc = 0;
  for (const SpherePoint& point : points) { ubexc += squared_norm(point.x, point.y, point.z); }

  if (Nrc == 0) {
    ucc   = 0.0;
    upexc = 0.0;
    ubexc = 0.0;
    return;
  }

  double pair_sum_images = 0.0;
  for (size_t i = 0; i + 1 < points.size(); ++i) {
    for (size_t j = i + 1; j < points.size(); ++j) {
      const double dx  = points[i].x - points[j].x;
      const double dy  = points[i].y - points[j].y;
      const double dz  = points[i].z - points[j].z;
      const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (rij == 0.0) { continue; }
      pair_sum_images += 1.0 / rij;
    }
  }

  ucc   = (3.0 / 5.0) * std::pow(Nrc, 2.0 / 3.0);
  upexc = pair_sum_images / static_cast<double>(Nrc);
  upexc -= ucc;

  ubexc /= static_cast<double>(Nrc);
  ubexc -= ucc;
  ubexc *= 0.5;
}

void Application::calculate()
{
  if (calculate_rdf) { compute_rdf_celllist(); }

  // for (size_t i = 0; i < N; ++i) { x[i] -= box.cx; }
  // for (size_t i = 0; i < N; ++i) { y[i] -= box.cy; }
  // for (size_t i = 0; i < N; ++i) { z[i] -= box.cz; }
  // for (size_t i = 0; i < N; ++i) { sr[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i]; }

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
  obs << timestep << "," << N << "," << box.V << "," << rmax << "," << dr << "," << nbins << "," << Ni << "," << Nri << "," << upexi << "," << ubexi
      << "," << uci << "," << Nc << "," << Nrc << "," << upexc << "," << ubexc << "," << ucc << "\n";
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
