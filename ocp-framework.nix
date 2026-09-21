{
  lib,
  stdenv,
  cmake,
  adios2,
}:

stdenv.mkDerivation {
  pname = "ocp-framework";
  version = "unstable";

  src = lib.cleanSource ./src;

  nativeBuildInputs = [ cmake ];
  buildInputs = [ adios2 ];

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
