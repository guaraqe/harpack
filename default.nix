{ mkDerivation, arpack, base, monad-loops, openblasCompat, stdenv
, storable-complex, tasty, tasty-hunit, tasty-quickcheck
, tasty-smallcheck, vector
}:
mkDerivation {
  pname = "harpack";
  version = "0.1.0";
  src = ./.;
  libraryHaskellDepends = [
    base monad-loops storable-complex tasty tasty-hunit vector
  ];
  librarySystemDepends = [ arpack openblasCompat ];
  testHaskellDepends = [
    base tasty tasty-hunit tasty-quickcheck tasty-smallcheck vector
  ];
  description = "An interface to ARPACK for sparse eigenvalue problems";
  license = stdenv.lib.licenses.bsd3;
}
