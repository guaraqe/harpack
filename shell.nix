let
  nixpkgs-src = builtins.fetchTarball {
    name = "nixos-19.03";
    url = "https://github.com/nixos/nixpkgs/archive/d011e4749457af484adf2e90062c83a44ad072a4.tar.gz";
    sha256 = "1nvhya0v98agidgxv83sb7xyb559qagx9f6iqknpxhxsv09j2qsp";
  };

  pkgs = import nixpkgs-src {};

in

(pkgs.haskellPackages.callPackage ./. {
  arpack = pkgs.arpack.override {
    openblas = pkgs.openblasCompat;
  };
}).env
