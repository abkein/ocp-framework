{
  lib,
  stdenv,
  cmake,
  adios2,
  python3,
  mpi,
  mpiCheckPhaseHook,
}:

stdenv.mkDerivation {
  pname = "ocp-framework";
  version = "unstable";

  src = lib.cleanSource ./src;

  nativeBuildInputs = [
    cmake
    python3
  ];
  buildInputs = [
    adios2
    mpi
  ];
  nativeCheckInputs = [ mpiCheckPhaseHook ];
  doCheck = true;
  checkPhase = ''
    runHook preCheck
    ctest --output-on-failure
    runHook postCheck
  '';

  # The CMake project does not define an install target.
  installPhase = ''
    runHook preInstall
    install -Dm755 ocp-framework "$out/bin/ocp-framework"
    runHook postInstall
  '';

  meta = {
    description = "Analyze molecular dynamics simulation dumps in ADIOS2 format";
    mainProgram = "ocp-framework";
    platforms = lib.platforms.linux;
  };
}
