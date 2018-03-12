with (import <nixpkgs> {});

(haskellPackages.callPackage ./. { arpack = arpack; }).env
