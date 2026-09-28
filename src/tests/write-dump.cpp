#include <adios2/common/ADIOSTypes.h>
#include <adios2/cxx/ADIOS.h>
#include <adios2/cxx/Variable.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Input: scaled, lengths[3], center[3], then frames of N followed by N (x,y,z,vx,vy,vz) rows.
int main(int argc, char** argv)
{
  if (argc < 2) { return 1; }
  bool vary_box  = false;
  size_t columns = 7;
  for (int i = 2; i < argc; ++i) {
    const std::string option(argv[i]);
    if (option == "--vary-box") {
      vary_box = true;
    } else if (option == "--legacy-columns") {
      columns = 5;
    } else {
      return 1;
    }
  }
  int32_t scaled = 0;
  std::array<double, 3> lengths{}, center{};
  std::cin >> scaled;
  for (auto& value : lengths) { std::cin >> value; }
  for (auto& value : center) { std::cin >> value; }

  adios2::ADIOS adios;
  auto io = adios.DeclareIO("fixture");
  io.SetEngine("BP5");
  io.DefineAttribute<int32_t>("scaled", scaled);
  io.DefineAttribute<int32_t>("triclinic", 0);
  auto atoms    = io.DefineVariable<double>("atoms", {1, columns}, {0, 0}, {1, columns});
  auto natoms   = io.DefineVariable<uint64_t>("natoms");
  auto timestep = io.DefineVariable<uint64_t>("ntimestep");
  const std::array<const char*, 6> names{"boxxlo", "boxxhi", "boxylo", "boxyhi", "boxzlo", "boxzhi"};
  std::vector<adios2::Variable<double>> bounds;
  bounds.reserve(names.size());
  for (const auto* const name : names) { bounds.push_back(io.DefineVariable<double>(name)); }
  auto writer = io.Open(argv[1], adios2::Mode::Write);
  uint64_t n{}, step = 0;
  while (std::cin >> n) {
    if (vary_box) {
      for (auto& length : lengths) { std::cin >> length; }
    }
    std::vector<double> data(n * columns);
    for (size_t i = 0; i < n; ++i) {
      data[columns * i] = i + 1;
      if (columns == 5) { data[columns * i + 1] = 1; }
      const size_t position_offset = columns == 7 ? 1 : 2;
      for (size_t axis = 0; axis < 3; ++axis) {
        double coordinate{};
        if (!(std::cin >> coordinate)) { throw std::runtime_error("Incomplete frame"); }
        data[columns * i + position_offset + axis] = (scaled != 0) ? coordinate / lengths[axis] + 0.5 : coordinate + center[axis];
      }
      for (size_t axis = 0; axis < 3; ++axis) {
        double velocity{};
        if (!(std::cin >> velocity)) { throw std::runtime_error("Incomplete velocities"); }
        if (columns == 7) { data[columns * i + 4 + axis] = velocity; }
      }
    }
    writer.BeginStep();
    writer.Put(natoms, n, adios2::Mode::Sync);
    writer.Put(timestep, step, adios2::Mode::Sync);
    for (size_t axis = 0; axis < 3; ++axis) {
      writer.Put(bounds[2 * axis], center[axis] - lengths[axis] / 2, adios2::Mode::Sync);
      writer.Put(bounds[2 * axis + 1], center[axis] + lengths[axis] / 2, adios2::Mode::Sync);
    }
    atoms.SetShape({n, columns});
    atoms.SetSelection({{0, 0}, {n, columns}});
    writer.Put(atoms, data.data(), adios2::Mode::Sync);
    writer.EndStep();
    ++step;
  }
  writer.Close();
}
