#ifndef OCP_MPI_CSV_HPP
#define OCP_MPI_CSV_HPP

#include <mpi.h>

#include <string>

void mpi_check(int result, const char* operation);

// All ranks open, append, and close files in the same order. Empty rows are
// valid, allowing ranks without a frame to participate in a partial batch.
class MpiCsv {
 public:
  MpiCsv()                         = default;
  MpiCsv(const MpiCsv&)            = delete;
  MpiCsv& operator=(const MpiCsv&) = delete;
  MpiCsv(MpiCsv&& other) noexcept;
  ~MpiCsv();

  void open(const std::string& filename, const std::string& header);
  void append(const std::string& row);
  void close();

 private:
  MPI_File file     = MPI_FILE_NULL;
  MPI_Offset offset = 0;
  int rank          = 0;
};

#endif
