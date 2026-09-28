#ifndef OCP_FRAMEWORK_HPP
#define OCP_FRAMEWORK_HPP

#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Engine.h>
#include <adios2/cxx/IO.h>
#include <adios2/cxx/Variable.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using size_t = std::size_t;

struct SphereSpec {
  double rfrac;
  bool exact_circumscribed = false;    // c: use 3*(L_min/2)^2 for the boundary test
};

struct Options {
  Options(int argc, char** argv)
  {
    const std::string usage = "Usage: ocp-framework <file.bp> [rfrac ...] [--config adios2_config.xml] [--io read_dump]\n"
                              "                              [--obs-out observables.csv] [--rdf-out rdf.csv]\n"
                              "                              [--meta-out metadata.csv] [--rdf dr]\n"
                              "                              [--spheres-out spheres.json] [--gamma Gamma]\n"
                              "  rfrac: multiplier in (0,3), or c for sqrt(3)\n";
    if (argc < 2) { throw std::runtime_error(usage); }
    const std::string firstarg = argv[1];
    if ((firstarg == "--help") || (firstarg == "-h")) { throw std::runtime_error(usage); }

    bpFile = firstarg;

    for (int i = 2; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "--config") {
        adiosConfig = need_value(argc, argv, i);
      } else if (a == "--io") {
        ioName = need_value(argc, argv, i);
      } else if (a == "--rdf") {
        calculate_rdf = true;
        dr            = std::stod(need_value(argc, argv, i));
        if (dr <= 0.0) { throw std::runtime_error("dr must be > 0"); }
      } else if (a == "--obs-out") {
        csv_out = need_value(argc, argv, i);
      } else if (a == "--rdf-out") {
        rdf_out = need_value(argc, argv, i);
      } else if (a == "--meta-out") {
        metadata_out = need_value(argc, argv, i);
      } else if (a == "--spheres-out") {
        spheres_out = need_value(argc, argv, i);
      } else if (a == "--gamma") {
        const auto value  = need_value(argc, argv, i);
        const char* error = "--gamma must be finite and > 0, with a finite positive 1/Gamma";
        size_t consumed   = 0;
        try {
          gamma = std::stod(value, &consumed);
        } catch (const std::logic_error&) {
          throw std::runtime_error(error);
        }
        if (consumed != value.size() || !std::isfinite(gamma) || gamma <= 0.0 || !std::isfinite(1.0 / gamma) || 1.0 / gamma <= 0.0) {
          throw std::runtime_error(error);
        }
      } else if (a == "c") {
        spheres.push_back({std::numbers::sqrt3, true});
      } else if (a.starts_with("-")) {
        throw std::runtime_error("Unknown argument: " + a);
      } else {
        const char* error = "Sphere multiplier must be finite and in (0,3), or c";
        size_t consumed   = 0;
        double rfrac{};
        try {
          rfrac = std::stod(a, &consumed);
        } catch (const std::logic_error&) {
          throw std::runtime_error(error);
        }
        if (consumed != a.size() || !std::isfinite(rfrac) || rfrac <= 0.0 || rfrac >= 3.0) { throw std::runtime_error(error); }
        spheres.push_back({rfrac});
      }
    }
    if (!calculate_rdf && spheres.empty()) { throw std::runtime_error("Need at least one sphere or --rdf specified."); }
  }

  double dr    = -1.0;
  double gamma = std::numeric_limits<double>::quiet_NaN();

  std::string bpFile;
  std::string adiosConfig  = "adios2_config.xml";
  std::string ioName       = "read_dump";

  std::string csv_out      = "observables.csv";
  std::string rdf_out      = "rdf.csv";
  std::string metadata_out = "metadata.csv";
  std::string spheres_out  = "spheres.json";
  std::vector<SphereSpec> spheres;

  bool calculate_rdf = false;

  [[nodiscard]] static std::string need_value(int argc, const char* const* const argv, int& i)
  {
    if (i + 1 >= argc) { throw std::runtime_error(std::string("Missing value after ") + argv[i]); }
    return argv[++i];
  }
};

struct Box {
  Box() = default;
  constexpr Box(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi) :
      xlo(xlo), xhi(xhi), ylo(ylo), yhi(yhi), zlo(zlo), zhi(zhi), Lx(xhi - xlo), Ly(yhi - ylo), Lz(zhi - zlo), V(Lx * Ly * Lz), cx(xlo + 0.5 * Lx),
      cy(ylo + 0.5 * Ly), cz(zlo + 0.5 * Lz), L_average((Lx + Ly + Lz) * (1.0 / 3.0)), L_min(std::min({Lx, Ly, Lz}))
  {
  }
  double xlo{}, xhi{}, ylo{}, yhi{}, zlo{}, zhi{};
  double Lx{}, Ly{}, Lz{};
  double V{};
  double cx{}, cy{}, cz{};    // box center, therefore center of inscribed/circumscribed spheres
  double L_average{};         // average edge length
  double L_min{};             // min edge length
};

struct SphereVirial {
  // Virial per mean ion, including the 1/3 used in the pressure equation.
  double virial          = std::numeric_limits<double>::quiet_NaN();
  double excess_pressure = std::numeric_limits<double>::quiet_NaN();
  double pressure        = std::numeric_limits<double>::quiet_NaN();
  double compressibility = std::numeric_limits<double>::quiet_NaN();
};

struct Sphere {
  SphereSpec spec;
  std::string filename;
  std::ofstream out_csv;
  double radius{}, radius_squared{}, mean_population{};
  size_t count{};
  double upex{}, ubex{}, uc{};
  SphereVirial virial;
};

struct Application {
  explicit Application(const Options& opt);
  Application(const Application&) = delete;
  Application(Application&&)      = delete;
  ~Application();

  void operator=(const Application&) = delete;
  void operator=(Application&&)      = delete;

  void run();
  void prepare_data();
  void modify_data();
  void calculate();
  void write();
  void compute_rdf_celllist();
  void inscribed(Sphere& sphere);
  void circumscribed(Sphere& sphere);
  void write_sphere_index();

  adios2::ADIOS adios;
  adios2::IO io;
  adios2::Engine reader;
  int32_t scaled    = 0;    // 0 if coordinates are dimensional (real); 1 if dimensionless — x->x/Lx, y->y/Ly, z->z/Lz
  int32_t triclinic = 0;
  std::ofstream out_rdf;
  std::ofstream out_meta;

  bool calculate_rdf = false;

  adios2::Variable<double> vAtoms;
  adios2::Variable<uint64_t> vNt;
  adios2::Variable<uint64_t> vNatoms;
  adios2::Variable<double> vXlo;
  adios2::Variable<double> vXhi;
  adios2::Variable<double> vYlo;
  adios2::Variable<double> vYhi;
  adios2::Variable<double> vZlo;
  adios2::Variable<double> vZhi;

  std::vector<double> atoms;
  std::vector<double> x, y, z;    // raw dimensional coordinates
  // std::vector<double> sx, sy, sz;
  std::vector<double> sr;    // squared distances from cell center
  std::vector<size_t> idlist;
  std::vector<double> counts;    // radial distribution function (rdf), but only raw counts, unnormalized

  uint64_t timestep = 0;
  size_t N          = 0;    // atom number in cell
  Box box;                  // cell
  double rmax  = 0;         // max distance to calculate rdf
  size_t nbins = 0;         // number of bins for rdf
  double dr    = 0;         // bin width

  std::vector<Sphere> spheres;
  std::string spheres_out;
  double gamma = std::numeric_limits<double>::quiet_NaN();
};

namespace COLS {
  constexpr int NUM_COLS  = 5;
  constexpr int ATOM_ID   = 0;
  constexpr int ATOM_TYPE = 1;
  constexpr int X         = 2;
  constexpr int Y         = 3;
  constexpr int Z         = 4;
}    // namespace COLS

struct CenterOfMass {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  CenterOfMass(const Box& box, const size_t N, const std::vector<double>& sfx, const std::vector<double>& sfy, const std::vector<double>& sfz)
  {
    if (N == 0) { return; }

    for (size_t i = 0; i < N; ++i) {
      x += sfx[i];
      y += sfy[i];
      z += sfz[i];
    }
    x *= box.Lx / static_cast<double>(N);
    y *= box.Ly / static_cast<double>(N);
    z *= box.Lz / static_cast<double>(N);
  }
};

#endif    // !OCP_FRAMEWORK_HPP
