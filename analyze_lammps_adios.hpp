#ifndef ANALYZE_LAMMPS_ADIOS_H
#define ANALYZE_LAMMPS_ADIOS_H

#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Engine.h>
#include <adios2/cxx/IO.h>
#include <adios2/cxx/Variable.h>

#include <algorithm>
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
static const double SQRTofTHREE_over_TWO = 0.5 * std ::sqrt(2.0);
// static const double ONEoverSQRTofTWO     = 1.0 / std::sqrt(2.0);

struct SDInfo {
  double rmax     = -1.0;    // if < 0: use rmaxFrac * min(L)
  double rmaxFrac = 0.5;     // default: min(L)/2
  double dr       = -1.0;    //
};

struct Options {
  Options(int argc, char** argv)
  {
    const std::string usage = "Usage: analyze_lammps_adios <file.bp> [--config adios2_config.xml] [--io read_dump]\n"
                              "                                      [--obs-out observables.csv] [--rdf-out rdf.csv]\n"
                              "                                      [--rdf dr] [--rmax 0.0] [--rmax-frac 0.5]\n"
                              "                                      [-i] [-c]\n";
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
      } else if (a == "--rmax") {
        info.rmax = std::stod(need_value(argc, argv, i));
        if (info.rmax <= 0.0) { throw std::runtime_error("--rmax must be > 0"); }
      } else if (a == "--rmax-frac") {
        info.rmaxFrac = std::stod(need_value(argc, argv, i));
        if (info.rmaxFrac <= 0.0) { throw std::runtime_error("--rmax-frac must be > 0"); }
      } else if (a == "--rdf") {
        calculate_rdf = true;
        info.dr       = std::stod(need_value(argc, argv, i));
        if (info.dr <= 0.0) { throw std::runtime_error("dr must be > 0"); }
      } else if (a == "--obs-out") {
        obsOut = need_value(argc, argv, i);
      } else if (a == "--rdf-out") {
        rdfOut = need_value(argc, argv, i);
      } else if (a == "-i") {
        calculate_inscribed = true;
      } else if (a == "-c") {
        calculate_circumscribed = true;
      } else {
        throw std::runtime_error("Unknown argument: " + a);
      }
    }
  }

  SDInfo info{};

  std::string bpFile;
  std::string adiosConfig      = "adios2_config.xml";
  std::string ioName           = "read_dump";

  std::string obsOut           = "observables.csv";
  std::string rdfOut           = "rdf.csv";

  bool calculate_rdf           = false;
  bool calculate_inscribed     = false;
  bool calculate_circumscribed = false;

  [[nodiscard]] static std::string need_value(int argc, const char* const* const argv, int& i)
  {
    if (i + 1 >= argc) { throw std::runtime_error(std::string("Missing value after ") + argv[i]); }
    return argv[++i];
  }
};

struct Box {
  Box() = default;
  constexpr Box(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi) :
      xlo(xlo), xhi(xhi), ylo(ylo), yhi(yhi), zlo(zlo), zhi(zhi), Lx(xhi - xlo), Ly(yhi - ylo), Lz(zhi - zlo)
  {
    V         = Lx * Ly * Lz;
    cx        = xlo + 0.5 * Lx;
    cy        = ylo + 0.5 * Ly;
    cz        = zlo + 0.5 * Lz;
    L_average = (Lx + Ly + Lz) * (1.0 / 3.0);
    L_min     = std::min({Lx, Ly, Lz});
  }
  double xlo{}, xhi{}, ylo{}, yhi{}, zlo{}, zhi{};
  double Lx{}, Ly{}, Lz{};
  double V{};
  double cx{}, cy{}, cz{};    // box center, therefore center of inscribed/circumscribed spheres
  double L_average{};         // average edge length
  double L_min{};             // min edge length
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
  int32_t scaled    = 0;    // 0 if coordinates are dimensional (real); 1 if dimensionless — x->x/Lx, y->y/Ly, z->z/Lz
  int32_t triclinic = 0;
  std::ofstream obs;
  std::ofstream rdf;

  bool calculate_rdf           = false;
  bool calculate_inscribed     = false;
  bool calculate_circumscribed = false;

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
  std::vector<double> sx, sy, sz;
  std::vector<double> sr;
  std::vector<size_t> idlist;
  std::vector<double> coords;
  std::vector<double> temp_coords;
  std::vector<bool> inside;
  std::vector<double> fcx, fcy, fcz;    // atom dimensional coordinates relative to the cell center
  std::vector<double> cx1, cy1, cz1;    // same as fc, but for image x->x-Lx
  std::vector<double> cx2, cy2, cz2;    // same as fc, but for image x->x+Lx
  std::vector<double> cx3, cy3, cz3;    // same as fc, but for image y->y-Ly
  std::vector<double> cx4, cy4, cz4;    // same as fc, but for image y->y+Ly
  std::vector<double> cx5, cy5, cz5;    // same as fc, but for image z->z-Lz
  std::vector<double> cx6, cy6, cz6;    // same as fc, but for image z->z+Lz
  std::array<std::array<std::vector<double>*, 3>, 7> vcs = {
      std::array<std::vector<double>*, 3>({&fcx, &fcy, &fcz}), std::array<std::vector<double>*, 3>({&cx1, &cy1, &cz1}),
      std::array<std::vector<double>*, 3>({&cx2, &cy2, &cz2}), std::array<std::vector<double>*, 3>({&cx3, &cy3, &cz3}),
      std::array<std::vector<double>*, 3>({&cx4, &cy4, &cz4}), std::array<std::vector<double>*, 3>({&cx5, &cy5, &cz5}),
      std::array<std::vector<double>*, 3>({&cx6, &cy6, &cz6})};    // for convinience of iteration
  std::vector<bool> ninside_original;                              // ninside[i] == false -> inside; whether atom is inside circumscribed sphere
  std::vector<bool> ninside1;    // ninside[i] == false -> inside; whether image x->x-Lx is inside circumscribed sphere
  std::vector<bool> ninside2;    // ninside[i] == false -> inside; whether image x->x+Lx is inside circumscribed sphere
  std::vector<bool> ninside3;    // ninside[i] == false -> inside; whether image y->y-Ly is inside circumscribed sphere
  std::vector<bool> ninside5;    // ninside[i] == false -> inside; whether image y->y+Ly is inside circumscribed sphere
  std::vector<bool> ninside4;    // ninside[i] == false -> inside; whether image z->z-Lz is inside circumscribed sphere
  std::vector<bool> ninside6;    // ninside[i] == false -> inside; whether image z->z+Lz is inside circumscribed sphere
  std::array<std::vector<bool>*, 7> ninside = {&ninside_original, &ninside1, &ninside2, &ninside3,
                                               &ninside4,         &ninside5, &ninside6};    // for convinience of iteration
  std::vector<double> counts;    // radial distribution function (rdf), but only raw counts, unnormalized

  uint64_t timestep = 0;
  size_t N          = 0;    // atom number in cell
  size_t N3         = 0;
  Box box;             // cell
  double rmax  = 0;    // max distance to calculate rdf
  size_t nbins = 0;    // number of bins for rdf
  double dr    = 0;    // bin width

  size_t Ni    = 0;    // number of atoms inside inscribed sphere, calculated from geometry
  size_t Nri   = 0;    // number of atoms inside inscribed sphere, counted from dataa
  size_t Nc    = 0;    // number of atoms inside circumscribed sphere, calculated from geometry
  size_t Nrc   = 0;    // number of atoms inside circumscribed sphere, counted from dataa
  double upexi = 0;    // energies for inscribed sphere
  double ubexi = 0;    // energies for inscribed sphere
  double uci   = 0;    // energies for inscribed sphere
  double upexc = 0;    // energies for circumscribed sphere
  double ubexc = 0;    // energies for circumscribed sphere
  double ucc   = 0;    // energies for circumscribed sphere
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
