#include "ocp-framework.hpp"

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
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
  // JSON filenames may contain quotes, backslashes, and control characters.
  [[nodiscard]] std::string json_string(const std::string& value)
  {
    constexpr std::array<char, 16> hex{{'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'}};
    std::string result = "\"";
    for (const unsigned char c : value) {
      if (c == '"' || c == '\\') {
        result += '\\';
        result += static_cast<char>(c);
      } else if (c < 0x20) {
        result += "\\u00";
        result += hex[c >> 4];
        result += hex[c & 15];
      } else {
        result += static_cast<char>(c);
      }
    }
    result += '"';
    return result;
  }

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
    adios(adios2::ADIOS(options.adiosConfig)), calculate_rdf(options.calculate_rdf), dr(options.dr), spheres_out(options.spheres_out),
    gamma(options.gamma)
{
  io     = adios.DeclareIO(options.ioName);
  reader = io.Open(options.bpFile, adios2::Mode::Read);

  // Numbered filenames preserve request order, including repeated multipliers.
  const std::filesystem::path base(options.csv_out);
  const std::string extension = base.has_extension() ? base.extension().string() : ".csv";
  for (const auto& spec : options.spheres) {
    Sphere sphere{.spec = spec};
    sphere.filename = (base.parent_path() / (base.stem().string() + "-" + std::to_string(spheres.size() + 1) + extension)).string();
    sphere.out_csv.open(sphere.filename);
    if (!sphere.out_csv) { throw std::runtime_error("Cannot open " + sphere.filename); }
    sphere.out_csv << std::setprecision(std::numeric_limits<double>::max_digits10);
    sphere.out_csv << "timestep,r,Ns,Nrs,upex,ubex,uc,virial,pex,p,Z\n";
    spheres.push_back(std::move(sphere));
  }
  write_sphere_index();

  out_meta = std::ofstream(options.metadata_out);
  if (!out_meta) { throw std::runtime_error("Cannot open " + options.metadata_out); }
  out_meta << std::setprecision(std::numeric_limits<double>::max_digits10);
  out_meta << "timestep,N,V,rmax,nbins,T\n";

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
}

void Application::modify_data()
{
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

  sr.resize(N);
  for (size_t i = 0; i < N; ++i) {
    x[i] -= box.cx;
    y[i] -= box.cy;
    z[i] -= box.cz;
    sr[i] = x[i] * x[i] + y[i] * y[i] + z[i] * z[i];
  }
}

void Application::inscribed(Sphere& sphere)
{
  const double Rsq = sphere.radius_squared;

  idlist.reserve(N);
  idlist.clear();
  for (size_t i = 0; i < N; ++i) {
    if (sr[i] < Rsq) {
      sphere.ubex += sr[i];
      idlist.emplace_back(i);
    }
  }
  sphere.count = idlist.size();

  for (size_t i = 0; i < sphere.count; ++i) {
    const auto ii   = idlist[i];
    const double xi = x[ii];
    const double yi = y[ii];
    const double zi = z[ii];
    for (size_t j = i + 1; j < sphere.count; ++j) {
      const auto jj    = idlist[j];
      const double dx  = xi - x[jj];
      const double dy  = yi - y[jj];
      const double dz  = zi - z[jj];
      const double rij = std::sqrt(dx * dx + dy * dy + dz * dz);
      sphere.upex += 1.0 / rij;
    }
  }
}

void Application::circumscribed(Sphere& sphere)
{
  const double Rsq = sphere.radius_squared;

  struct Point {
    double x, y, z;
  };
  std::vector<Point> points;
  points.reserve(N);

  // rfrac < 3 and a strict sphere boundary require at most the 26 neighboring cells.
  for (int ix = -1; ix <= 1; ++ix) {
    for (int iy = -1; iy <= 1; ++iy) {
      for (int iz = -1; iz <= 1; ++iz) {
        for (size_t i = 0; i < N; ++i) {
          const Point p{.x = x[i] + ix * box.Lx, .y = y[i] + iy * box.Ly, .z = z[i] + iz * box.Lz};
          const double rsq = p.x * p.x + p.y * p.y + p.z * p.z;
          if (rsq < Rsq) {
            points.push_back(p);
            sphere.ubex += rsq;
          }
        }
      }
    }
  }
  sphere.count = points.size();

  // Every selected image is a distinct point, including images of the same
  // original ion. Use direct distances, without a minimum-image reduction.
  for (size_t i = 0; i < sphere.count; ++i) {
    for (size_t j = i + 1; j < sphere.count; ++j) {
      const double dx = points[i].x - points[j].x;
      const double dy = points[i].y - points[j].y;
      const double dz = points[i].z - points[j].z;
      sphere.upex += 1.0 / std::sqrt(dx * dx + dy * dy + dz * dz);
    }
  }
}

void Application::calculate()
{
  rmax = 0.5 * box.L_min;
  if (calculate_rdf) { compute_rdf_celllist(); }

  for (auto& sphere : spheres) {
    sphere.radius          = sphere.spec.rfrac * rmax;
    sphere.radius_squared  = sphere.spec.exact_circumscribed ? 3.0 * (rmax * rmax) : sphere.radius * sphere.radius;
    sphere.mean_population = N_in_sphere(N, box, sphere.radius);
    sphere.upex = sphere.ubex = 0.0;
    if (sphere.spec.rfrac <= 1.0) {
      inscribed(sphere);
    } else {
      circumscribed(sphere);
    }
    sphere.virial = sphere_virial(sphere.mean_population, sphere.upex, sphere.ubex, static_cast<double>(N) / box.V, gamma);
    normalize_sphere(sphere.mean_population, sphere.upex, sphere.ubex, sphere.uc);
  }
}

void Application::write()
{
  out_meta << timestep << "," << N << "," << box.V << "," << rmax << "," << nbins << "," << 1.0 / gamma << "\n";
  out_meta.flush();
  for (auto& sphere : spheres) {
    sphere.out_csv << timestep << "," << sphere.radius << "," << sphere.mean_population << "," << sphere.count << "," << sphere.upex << ","
                   << sphere.ubex << "," << sphere.uc << "," << sphere.virial.virial << "," << sphere.virial.excess_pressure << ","
                   << sphere.virial.pressure << "," << sphere.virial.compressibility << "\n";
    sphere.out_csv.flush();
  }
  if (calculate_rdf) {
    out_rdf << timestep;
    for (size_t i = 0; i < nbins; ++i) { out_rdf << "," << counts[i]; }
    out_rdf << "\n";
    out_rdf.flush();
  }
}

void Application::write_sphere_index()
{
  std::ofstream index(spheres_out);
  if (!index) { throw std::runtime_error("Cannot open " + spheres_out); }
  index << std::setprecision(std::numeric_limits<double>::max_digits10) << "[\n";
  const auto directory = std::filesystem::absolute(spheres_out).parent_path().lexically_normal();
  for (size_t i = 0; i < spheres.size(); ++i) {
    const auto& sphere  = spheres[i];
    const auto filename = std::filesystem::absolute(sphere.filename).lexically_normal().lexically_relative(directory).generic_string();
    index << "  {\"rfrac\": " << sphere.spec.rfrac << ", \"filename\": " << json_string(filename) << "}";
    index << (i + 1 == spheres.size() ? "\n" : ",\n");
  }
  index << "]\n";
}

void Application::run()
{
  while (true) {
    const auto status = reader.BeginStep();
    if (status != adios2::StepStatus::OK) { break; }

    prepare_data();
    reader.EndStep();

    modify_data();
    calculate();
    write();
  }
}

void Application::compute_rdf_celllist()
{
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
      if (distance >= rmax) { continue; }
      const auto ibin = static_cast<size_t>(distance / dr);    // floor
      counts[ibin < nbins ? ibin : nbins - 1] += 2;
    }
  }
  for (size_t i = 0; i < nbins; ++i) { counts[i] /= static_cast<double>(N); }
}

Application::~Application()
{
  for (auto& sphere : spheres) { sphere.out_csv.close(); }
  out_meta.close();
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
