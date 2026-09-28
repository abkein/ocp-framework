#include "mpi-csv.hpp"

#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>

void mpi_check(const int result, const char* operation)
{
  if (result == MPI_SUCCESS) { return; }
  char message[MPI_MAX_ERROR_STRING];
  int length = 0;
  MPI_Error_string(result, message, &length);
  throw std::runtime_error(std::string(operation) + ": " + std::string(message, length));
}

MpiCsv::MpiCsv(MpiCsv&& other) noexcept : file(std::exchange(other.file, MPI_FILE_NULL)), offset(other.offset), rank(other.rank) {}

MpiCsv::~MpiCsv()
{
  // Never enter a collective while unwinding a rank-local exception: main()
  // aborts the communicator so peers cannot hang in their next collective.
  if (file != MPI_FILE_NULL && std::uncaught_exceptions() == 0) { MPI_File_close(&file); }
}

void MpiCsv::open(const std::string& filename, const std::string& header)
{
  mpi_check(MPI_Comm_rank(MPI_COMM_WORLD, &rank), "MPI_Comm_rank");
  mpi_check(MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file), "MPI_File_open");
  mpi_check(MPI_File_set_errhandler(file, MPI_ERRORS_RETURN), "MPI_File_set_errhandler");
  mpi_check(MPI_File_set_size(file, 0), "MPI_File_set_size");
  offset = 0;
  // Exactly one header, with the same collective path as data rows.
  append(rank == 0 ? header : std::string{});
}

void MpiCsv::append(const std::string& row)
{
  if (row.size() > static_cast<size_t>(std::numeric_limits<int>::max())) { throw std::runtime_error("CSV row exceeds the MPI count limit"); }
  const MPI_Offset bytes = static_cast<MPI_Offset>(row.size());
  MPI_Offset prefix = 0, total = 0;
  mpi_check(MPI_Exscan(&bytes, &prefix, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD), "MPI_Exscan");
  if (rank == 0) { prefix = 0; }
  mpi_check(MPI_Allreduce(&bytes, &total, 1, MPI_OFFSET, MPI_SUM, MPI_COMM_WORLD), "MPI_Allreduce");
  if (total < 0 || offset > std::numeric_limits<MPI_Offset>::max() - total) { throw std::runtime_error("CSV file exceeds the MPI offset limit"); }
  mpi_check(MPI_File_write_at_all(file, offset + prefix, row.data(), static_cast<int>(row.size()), MPI_BYTE, MPI_STATUS_IGNORE),
            "MPI_File_write_at_all");
  offset += total;
}

void MpiCsv::close()
{
  if (file != MPI_FILE_NULL) { mpi_check(MPI_File_close(&file), "MPI_File_close"); }
}
