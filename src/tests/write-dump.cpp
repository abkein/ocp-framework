#include <adios2.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

// Input: scaled, lengths[3], center[3], then frames of N followed by N xyz triples.
int main(int argc, char** argv)
{
  if (argc != 2) { return 1; }
  int32_t scaled;
  std::array<double, 3> lengths, center;
  std::cin >> scaled;
  for (auto& value : lengths) { std::cin >> value; }
  for (auto& value : center) { std::cin >> value; }

  adios2::ADIOS adios;
  auto io = adios.DeclareIO("fixture");
  io.SetEngine("BP5");
  io.DefineAttribute<int32_t>("scaled", scaled);
  io.DefineAttribute<int32_t>("triclinic", 0);
  auto atoms    = io.DefineVariable<double>("atoms", {1, 5}, {0, 0}, {1, 5});
  auto natoms   = io.DefineVariable<uint64_t>("natoms");
  auto timestep = io.DefineVariable<uint64_t>("ntimestep");
  const std::array<const char*, 6> names{"boxxlo", "boxxhi", "boxylo", "boxyhi", "boxzlo", "boxzhi"};
  std::vector<adios2::Variable<double>> bounds;
  for (const auto name : names) { bounds.push_back(io.DefineVariable<double>(name)); }
  auto writer = io.Open(argv[1], adios2::Mode::Write);
  uint64_t n, step = 0;
  while (std::cin >> n) {
    std::vector<double> data(n * 5);
    for (size_t i = 0; i < n; ++i) {
      data[5 * i]     = i + 1;
      data[5 * i + 1] = 1;
      for (size_t axis = 0; axis < 3; ++axis) {
        double coordinate;
        if (!(std::cin >> coordinate)) { throw std::runtime_error("Incomplete frame"); }
        data[5 * i + 2 + axis] = scaled ? coordinate / lengths[axis] + 0.5 : coordinate + center[axis];
      }
    }
    writer.BeginStep();
    writer.Put(natoms, n, adios2::Mode::Sync);
    writer.Put(timestep, step, adios2::Mode::Sync);
    for (size_t axis = 0; axis < 3; ++axis) {
      writer.Put(bounds[2 * axis], center[axis] - lengths[axis] / 2, adios2::Mode::Sync);
      writer.Put(bounds[2 * axis + 1], center[axis] + lengths[axis] / 2, adios2::Mode::Sync);
    }
    atoms.SetShape({n, 5});
    atoms.SetSelection({{0, 0}, {n, 5}});
    writer.Put(atoms, data.data(), adios2::Mode::Sync);
    writer.EndStep();
    ++step;
  }
  writer.Close();
}
