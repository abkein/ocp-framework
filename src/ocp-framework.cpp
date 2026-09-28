#include "ocp-framework.hpp"
#include "mpi-csv.hpp"

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
#include <mpi.h>
#include <numbers>
#include <sstream>
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
                                           const double temperature)
  {
    SphereVirial result;
    if (mean_population == 0.0) { return result; }
    // Coulomb pairs contribute +1/r; the background force -r contributes -r^2.
    result.virial          = (pair_sum - radius_sum) / (3.0 * mean_population);
    result.excess_pressure = density * result.virial;
    if (std::isfinite(temperature)) {
      result.pressure = density * (temperature + result.virial);
      if (temperature > 0.0) { result.compressibility = 1.0 + result.virial / temperature; }
    }
    return result;
  }
}    // namespace

Application::Application(const Options& options) :
    adios(adios2::ADIOS(options.adiosConfig, MPI_COMM_SELF)), calculate_rdf(options.calculate_rdf), dr(options.dr), spheres_out(options.spheres_out),
    mass(options.mass)
{
  mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "MPI_Comm_rank");
  mpi_check(MPI_Comm_size(MPI_COMM_WORLD, &world_size), "MPI_Comm_size");
  io     = adios.DeclareIO(options.ioName);
  reader = io.Open(options.bpFile, adios2::Mode::ReadRandomAccess);

  // Numbered filenames preserve request order, including repeated multipliers.
  const std::filesystem::path base(options.csv_out);
  const std::string extension = base.has_extension() ? base.extension().string() : ".csv";
  for (const auto& spec : options.spheres) {
    Sphere sphere{.spec     = spec,
                  .filename = (base.parent_path() / (base.stem().string() + "-" + std::to_string(spheres.size() + 1) + extension)).string()};
    sphere.out_csv.open(sphere.filename, "timestep,r,Ns,Nrs,upex,ubex,uc,virial,pex,p,Z\n");
    spheres.push_back(std::move(sphere));
  }
  if (rank == 0) { write_sphere_index(); }

  out_meta.open(options.metadata_out, "timestep,N,V,rmax,nbins,T,cvx,cvy,cvz\n");
  if (calculate_rdf) { out_rdf.open(options.rdf_out, ""); }
}

void Application::prepare_data(const size_t step)
{
  // Random-access mode exposes metadata for all steps after Open().
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

  vNt.SetStepSelection({step, 1});
  vNatoms.SetStepSelection({step, 1});
  vAtoms.SetStepSelection({step, 1});
  for (auto* variable : {&vXlo, &vXhi, &vYlo, &vYhi, &vZlo, &vZhi}) { variable->SetStepSelection({step, 1}); }

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

  box              = Box(xlo, xhi, ylo, yhi, zlo, zhi);

  // Reject old five-column dumps instead of interpreting type/coordinates as velocities.
  const auto shape = vAtoms.Shape(step);
  if (shape.size() != 2 || shape[0] != N || shape[1] != COLS::NUM_COLS) {
    throw std::runtime_error("Expected atoms shape [natoms,7] with columns id,x,y,z,vx,vy,vz");
  }
  // The atom count and shape may change each step.
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

  // columns: [id, x, y, z, vx, vy, vz]
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
  vx.resize(N);
  vy.resize(N);
  vz.resize(N);
  for (size_t i = 0; i < N; ++i) {
    vx[i] = atoms[COLS::NUM_COLS * i + COLS::VX];
    vy[i] = atoms[COLS::NUM_COLS * i + COLS::VY];
    vz[i] = atoms[COLS::NUM_COLS * i + COLS::VZ];
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

void Application::compute_temperature()
{
  temperature = std::numeric_limits<double>::quiet_NaN();
  if (N <= 1) { return; }

  // LAMMPS thermo_temp in 3D, k_B=1, no constraint fixes: 3N-3 DOFs.
  // Use the original velocities without box scaling or COM subtraction.
  double squared_speed_sum = 0.0;
  for (size_t i = 0; i < N; ++i) { squared_speed_sum += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]; }
  temperature = mass * squared_speed_sum / (3.0 * static_cast<double>(N - 1));
}

void Application::compute_center_of_mass_velocity()
{
  cvx = cvy = cvz = std::numeric_limits<double>::quiet_NaN();
  if (N == 0) { return; }

  // Equal ion masses cancel in total momentum / total mass.
  cvx = cvy = cvz = 0.0;
  for (size_t i = 0; i < N; ++i) {
    cvx += vx[i];
    cvy += vy[i];
    cvz += vz[i];
  }
  cvx /= static_cast<double>(N);
  cvy /= static_cast<double>(N);
  cvz /= static_cast<double>(N);
}

void Application::calculate()
{
  compute_temperature();
  compute_center_of_mass_velocity();
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
    sphere.virial = sphere_virial(sphere.mean_population, sphere.upex, sphere.ubex, static_cast<double>(N) / box.V, temperature);
    normalize_sphere(sphere.mean_population, sphere.upex, sphere.ubex, sphere.uc);
  }
}

void Application::write(const bool active)
{
  std::ostringstream row;
  row << std::setprecision(std::numeric_limits<double>::max_digits10);
  if (active) {
    row << timestep << "," << N << "," << box.V << "," << rmax << "," << nbins << "," << temperature << "," << cvx << "," << cvy << "," << cvz
        << "\n";
  }
  out_meta.append(row.str());
  for (auto& sphere : spheres) {
    row.str("");
    if (active) {
      row << timestep << "," << sphere.radius << "," << sphere.mean_population << "," << sphere.count << "," << sphere.upex << "," << sphere.ubex
          << "," << sphere.uc << "," << sphere.virial.virial << "," << sphere.virial.excess_pressure << "," << sphere.virial.pressure << ","
          << sphere.virial.compressibility << "\n";
    }
    sphere.out_csv.append(row.str());
  }
  if (calculate_rdf) {
    row.str("");
    if (active) {
      row << timestep;
      for (size_t i = 0; i < nbins; ++i) { row << "," << counts[i]; }
      row << "\n";
    }
    out_rdf.append(row.str());
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
  const size_t steps  = reader.Steps();
  const size_t stride = static_cast<size_t>(world_size);
  for (size_t batch = 0; batch < steps; batch += stride) {
    const size_t step = batch + static_cast<size_t>(rank);
    const bool active = step < steps;
    if (active) {
      prepare_data(step);
      modify_data();
      calculate();
    }
    // All ranks participate, including those without a frame in the last batch.
    write(active);
  }
  reader.Close();
  for (auto& sphere : spheres) { sphere.out_csv.close(); }
  out_meta.close();
  if (calculate_rdf) { out_rdf.close(); }
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

int main(int argc, char** argv)
{
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) { return 1; }
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  try {
    mpi_check(MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN), "MPI_Comm_set_errhandler");
    Application app(Options(argc, argv));
    app.run();
  } catch (const std::exception& e) {
    std::cerr << "Rank " << rank << " error: " << e.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  return MPI_Finalize() == MPI_SUCCESS ? 0 : 1;
}
