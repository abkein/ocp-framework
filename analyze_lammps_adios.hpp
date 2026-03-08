#ifndef ANALYZE_LAMMPS_ADIOS_H
#define ANALYZE_LAMMPS_ADIOS_H

// #include <adios2/common/ADIOSTypes.h>
#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Engine.h>
#include <adios2/cxx/IO.h>
#include <adios2/cxx/Variable.h>

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using size_t = std::size_t;

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
                               "                                      [--rmax 0.0] [--rmax-frac 0.5] [--cell-size 0.0]\n"
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
  double xlo, xhi, ylo, yhi, zlo, zhi;
  [[nodiscard]] double Lx() const { return xhi - xlo; }
  [[nodiscard]] double Ly() const { return yhi - ylo; }
  [[nodiscard]] double Lz() const { return zhi - zlo; }
  [[nodiscard]] double V() const { return Lx() * Ly() * Lz(); }
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
  void compute_rdf_celllist(const double dr, const size_t nbins);

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
  std::vector<double> counts;

  uint64_t timestep = 0;
  size_t N          = 0;
  Box box{};
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
    x *= box.Lx() / static_cast<double>(N);
    y *= box.Ly() / static_cast<double>(N);
    z *= box.Lz() / static_cast<double>(N);
  }
};

#endif    // !ANALYZE_LAMMPS_ADIOS_H
