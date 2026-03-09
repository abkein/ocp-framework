#ifndef ANALYZE_LAMMPS_ADIOS_H
#define ANALYZE_LAMMPS_ADIOS_H

// #include <adios2/common/ADIOSTypes.h>
#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Engine.h>
#include <adios2/cxx/IO.h>
#include <adios2/cxx/Variable.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using size_t                             = std::size_t;

static constexpr double PI               = 3.141592653589793;
static const double SQRTofTWO_over_THREE = std::sqrt(2.0) / 3.0;
static const double ONEoverSQRTofTWO     = 1.0 / std::sqrt(2.0);

struct SDInfo {
  double rmax     = -1.0;    // if < 0: use rmaxFrac * min(L)
  double rmaxFrac = 0.5;     // default: min(L)/2
  double dr       = -1.0;    // if < 0: use rmax
};

struct Options {
  Options(int argc, char** argv)
  {
    if (argc < 2) {
      throw std::runtime_error("Usage: analyze_lammps_adios <file.bp> [--config adios2_config.xml] [--io read_dump]\n"
                               "                                      [--rmax 0.0] [--rmax-frac 0.5] [--dr 0.0]\n"
                               "                                      [--obs-out observables.csv] [--rdf-out rdf.csv]\n");
    }

    bpFile = argv[1];

    for (int i = 2; i < argc; ++i) {
      const std::string a = argv[i];
      if (a == "--config") {
        adiosConfig = need_value(argc, argv, i);
      } else if (a == "--io") {
        ioName = need_value(argc, argv, i);
      } else if (a == "--rmax") {
        info.rmax = std::stod(need_value(argc, argv, i));
      } else if (a == "--rmax-frac") {
        info.rmaxFrac = std::stod(need_value(argc, argv, i));
      } else if (a == "--dr") {
        info.dr = std::stod(need_value(argc, argv, i));
      } else if (a == "--obs-out") {
        obsOut = need_value(argc, argv, i);
      } else if (a == "--rdf-out") {
        rdfOut = need_value(argc, argv, i);
      } else {
        throw std::runtime_error("Unknown argument: " + a);
      }
    }

    if (info.rmaxFrac <= 0.0) { throw std::runtime_error("--rmax-frac must be > 0"); }
  }

  std::string bpFile;
  std::string adiosConfig = "adios2_config.xml";
  std::string ioName      = "read_dump";

  SDInfo info{};

  std::string obsOut = "observables.csv";
  std::string rdfOut = "rdf.csv";

  [[nodiscard]] static std::string need_value(int argc, const char* const* const argv, int& i)
  {
    if (i + 1 >= argc) { throw std::runtime_error(std::string("Missing value after ") + argv[i]); }
    return argv[++i];
  }
};

struct Box {
  Box() = default;
  Box(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi) :
      xlo(xlo), xhi(xhi), ylo(ylo), yhi(yhi), zlo(zlo), zhi(zhi), Lx(xhi - xlo), Ly(yhi - ylo), Lz(zhi - zlo)
  {
    V  = Lx * Ly * Lz;
    cx = xlo + 0.5 * Lx;
    cy = ylo + 0.5 * Ly;
    cz = zlo + 0.5 * Lz;
  }
  double xlo{}, xhi{}, ylo{}, yhi{}, zlo{}, zhi{};
  double Lx{}, Ly{}, Lz{};
  double V{};
  double cx{}, cy{}, cz{};
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
  void calculate();
  void write();
  void compute_rdf_celllist();
  void inscribed();
  void circumscribed();

  SDInfo info{};

  adios2::ADIOS adios;
  adios2::IO io;
  adios2::Engine reader;
  int32_t scaled    = 0;
  int32_t triclinic = 0;
  std::ofstream obs;
  std::ofstream rdf;

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
  std::vector<double> x, y, z;
  std::vector<double> fcx, fcy, fcz;
  std::vector<double> cx1, cy1, cz1;
  std::vector<double> cx2, cy2, cz2;
  std::vector<double> cx3, cy3, cz3;
  std::vector<double> cx4, cy4, cz4;
  std::vector<double> cx5, cy5, cz5;
  std::vector<double> cx6, cy6, cz6;
  std::array<std::array<std::vector<double>*, 3>, 7> vcs = {
      std::array<std::vector<double>*, 3>({&fcx, &fcy, &fcz}), std::array<std::vector<double>*, 3>({&cx1, &cy1, &cz1}),
      std::array<std::vector<double>*, 3>({&cx2, &cy2, &cz2}), std::array<std::vector<double>*, 3>({&cx3, &cy3, &cz3}),
      std::array<std::vector<double>*, 3>({&cx4, &cy4, &cz4}), std::array<std::vector<double>*, 3>({&cx5, &cy5, &cz5}),
      std::array<std::vector<double>*, 3>({&cx6, &cy6, &cz6})};
  std::vector<bool> ninside_original;
  std::vector<bool> ninside1;
  std::vector<bool> ninside2;
  std::vector<bool> ninside3;
  std::vector<bool> ninside4;
  std::vector<bool> ninside5;
  std::vector<bool> ninside6;
  std::array<std::vector<bool>*, 7> ninside = {&ninside_original, &ninside1, &ninside2, &ninside3, &ninside4, &ninside5, &ninside6};
  std::vector<double> counts;

  uint64_t timestep = 0;
  size_t N          = 0;
  Box box;
  double rmax  = 0;
  size_t nbins = 0;
  double dr    = 0;

  size_t Ni    = 0;
  size_t Nc    = 0;
  size_t Nri    = 0;
  size_t Nrc    = 0;
  double upexi = 0;
  double ubexi = 0;
  double uci   = 0;
  double upexc = 0;
  double ubexc = 0;
  double ucc   = 0;
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

#endif    // !ANALYZE_LAMMPS_ADIOS_H
