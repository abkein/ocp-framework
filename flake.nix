{
  inputs = {
    system-flake.url = "path:/home/kein/nixos-configuration";
    nixpkgs.follows = "system-flake/nixpkgs";
    flake-utils.follows = "system-flake/flake-utils";
  };

  outputs =
    { nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };
        ocp-framework = pkgs.callPackage ./ocp-framework.nix { };
      in
      {
        packages = {
          ocp-framework = ocp-framework;
          default = ocp-framework;
        };
      }
    );
}
