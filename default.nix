with import <nixpkgs> {};

stdenv.mkDerivation rec {
  name = "astrochem-env";
  env = buildEnv { name = name; paths = buildInputs; };
  python2Env = python2.withPackages (ps: with ps; [ numpy h5py ]);
  buildInputs = [
    python2Env
    gfortran
    hdf5
    autoconf
    automake
    libtool
    sundials
  ];
}
