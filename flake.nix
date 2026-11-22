{
  description = "Bayesian local FDR for sparse count data: code for Gauran, Park, Rattsev, Peterson, Kann, Park (2022), Ann. Appl. Stat. 16(3):1459-1475";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/f3aa7d3d72c2bc66e0f30f4262d7ec2e33c9eab2";

  outputs =
    { self, nixpkgs }:
    let
      lib = nixpkgs.lib;
      systems = [
        "aarch64-darwin"
        "x86_64-darwin"
        "aarch64-linux"
        "x86_64-linux"
      ];
      forEachSystem = f: lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});

      # R with the packages the scripts load: ggplot2 draws Figures S1 and S3,
      # openxlsx reads the output files of the original runs
      # (archive/reduce_archive.R), testthat runs the unit tests, Rcpp
      # compiles the C++ files of rcpp/src/.
      rEnv =
        pkgs:
        pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            ggplot2
            openxlsx
            testthat
            Rcpp
          ];
        };

      # Python with the packages of python/localfdr; numba compiles the random
      # number generator and the loops that follow R's order of operations.
      pyEnv =
        pkgs:
        pkgs.python3.withPackages (
          ps: with ps; [
            numpy
            pandas
            matplotlib
            numba
            pytest
          ]
        );

      # The BLAS and LAPACK library R is linked to; python/localfdr/rbase.py
      # calls it for its matrix products and decompositions. R is built on
      # pkgs.blas, whose provider is OpenBLAS with 32-bit integers on every
      # system; pkgs.openblas itself has 64-bit integers on x86_64, another
      # library with another calling convention.
      blasLibrary =
        pkgs: "${pkgs.blas.provider}/lib/libopenblas${pkgs.stdenv.hostPlatform.extensions.sharedLibrary}";

      # One command per table and figure of the paper.
      scripts = {
        table1 = "table1_simulation";
        table2 = "table2_prior_sensitivity";
        table3 = "table3_protein_domains";
        fig1 = "fig1_domain_barplots";
        figS1 = "figS1_fdr_tpr_well_separated";
        figS2 = "figS2_histograms";
        figS3 = "figS3_fdr_table";
      };

      # Run the scripts from the flake source; results/ and figures/ are
      # written to the directory the command is started in. The scripts run
      # their chains and replications in parallel R processes (--cores), so
      # each process gets one BLAS thread.
      runner =
        pkgs: name: names:
        pkgs.writeShellScript "localfdr-${name}" ''
          set -euo pipefail
          export LOCALFDR_OUTPUT="''${LOCALFDR_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${rEnv pkgs}/bin/Rscript "scripts/$script.R" "$@"
          done
        '';
      app = pkgs: name: names: {
        type = "app";
        program = "${runner pkgs name names}";
        meta.description = "Build ${name} of the paper";
      };

      # The same commands with the Python scripts, as py-<name>.
      pyRunner =
        pkgs: name: names:
        pkgs.writeShellScript "localfdr-py-${name}" ''
          set -euo pipefail
          export LOCALFDR_OUTPUT="''${LOCALFDR_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
          export LOCALFDR_BLAS=${blasLibrary pkgs}
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${pyEnv pkgs}/bin/python "python/scripts/$script.py" "$@"
          done
        '';
      pyApp = pkgs: name: names: {
        type = "app";
        program = "${pyRunner pkgs name names}";
        meta.description = "Build ${name} of the paper with the Python scripts";
      };

      # The same commands with the Rcpp scripts, as rcpp-<name>. The C++
      # files are compiled on the first run (R calls the compiler it was
      # built with, through make) and kept in the temporary folder of the
      # machine, or in LOCALFDR_RCPP_CACHE.
      compilerPath = pkgs: "${pkgs.gnumake}/bin:${pkgs.stdenv.cc}/bin";
      rcppRunner =
        pkgs: name: names:
        pkgs.writeShellScript "localfdr-rcpp-${name}" ''
          set -euo pipefail
          export LOCALFDR_OUTPUT="''${LOCALFDR_OUTPUT:-$PWD}"
          export OPENBLAS_NUM_THREADS=1
          export PATH=${compilerPath pkgs}:$PATH
          cd ${self}
          for script in ${lib.concatStringsSep " " names}; do
            ${rEnv pkgs}/bin/Rscript "rcpp/scripts/$script.R" "$@"
          done
        '';
      rcppApp = pkgs: name: names: {
        type = "app";
        program = "${rcppRunner pkgs name names}";
        meta.description = "Build ${name} of the paper with the Rcpp scripts";
      };

      # The unit tests: testthat over tests/testthat/ and rcpp/tests/testthat/
      # and pytest over python/tests/, started in the current directory. The
      # Python tests marked against_r run R and compare with it, so R is on
      # the PATH too.
      unitTests = pkgs: ''
        export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
        export LOCALFDR_BLAS=${blasLibrary pkgs}
        export PYTHONDONTWRITEBYTECODE=1
        export PATH=${rEnv pkgs}/bin:${pyEnv pkgs}/bin:${compilerPath pkgs}:$PATH
        unit_r() { Rscript tests/testthat.R; }
        unit_rcpp() { Rscript rcpp/tests/testthat.R; }
        unit_python() { python -m pytest -p no:cacheprovider python/tests "$@"; }
      '';

      # The text `nix run .#help` prints: every command of the flake.
      helpText =
        let
          row =
            command: text:
            "  nix run .#${command}${lib.fixedWidthString (12 - builtins.stringLength command) " " ""}${text}";
          rows = lib.mapAttrsToList (name: script: row name "scripts/${script}.R") scripts;
          pyRows = lib.mapAttrsToList (name: script: row "py-${name}" "python/scripts/${script}.py") scripts;
          rcppRows = lib.mapAttrsToList (
            name: script: row "rcpp-${name}" "rcpp/scripts/${script}.R"
          ) scripts;
        in
        lib.concatStringsSep "\n" (
          [
            "BLISSFULHotspots: code of Gauran, Park, Rattsev, Peterson, Kann, Park (2022), Ann. Appl. Stat. 16(3):1459-1475"
            ""
            "One command per table and figure of the paper (R):"
          ]
          ++ rows
          ++ [
            (row "all" "every table and figure")
            ""
            "The same outputs with the Python scripts:"
          ]
          ++ pyRows
          ++ [
            (row "py-all" "every table and figure")
            ""
            "The same outputs with the Rcpp scripts (the sampler in C++):"
          ]
          ++ rcppRows
          ++ [
            (row "rcpp-all" "every table and figure")
            ""
            "Other commands:"
            (row "test" "unit tests of R/, rcpp/src/ and python/localfdr/ in the current directory")
            (row "versions" "print the software versions")
            (row "help" "print this list (also: nix run . -- --help)")
            "  nix flake check       build every output, compare it with the paper and R with Python and Rcpp, run the unit tests"
            "  nix develop           a shell with the same R and Python"
            ""
            "Results are written to results/ and figures/ of the current directory (or of LOCALFDR_OUTPUT)."
            "Every table and figure command takes the options listed in the header of its script, after --."
            "The simulation studies (table1, table2, figS1, figS3) read the archived runs unless --rerun is"
            "given, for example"
            "  nix run .#table1 -- --rerun --runs=1:10 --cores=4"
            "The test command passes its arguments to pytest, for example: nix run .#test -- -k rng"
          ]
        );
      help = pkgs: {
        type = "app";
        program = "${pkgs.writeShellScript "localfdr-help" ''
          printf '%b\n' "${helpText}"
        ''}";
        meta.description = "List every command";
      };
    in
    {
      devShells = forEachSystem (pkgs: {
        default = pkgs.mkShell {
          packages = [
            (rEnv pkgs)
            (pyEnv pkgs)
          ];
          OPENBLAS_NUM_THREADS = 1;
          VECLIB_MAXIMUM_THREADS = 1;
          LOCALFDR_BLAS = blasLibrary pkgs;
        };
      });

      apps = forEachSystem (
        pkgs:
        lib.mapAttrs (name: script: app pkgs name [ script ]) scripts
        // lib.mapAttrs' (name: script: lib.nameValuePair "py-${name}" (pyApp pkgs name [ script ])) scripts
        // lib.mapAttrs' (
          name: script: lib.nameValuePair "rcpp-${name}" (rcppApp pkgs name [ script ])
        ) scripts
        // {
          all = app pkgs "all" (lib.attrValues scripts);
          py-all = pyApp pkgs "all" (lib.attrValues scripts);
          rcpp-all = rcppApp pkgs "all" (lib.attrValues scripts);
          versions = {
            type = "app";
            program = "${pkgs.writeShellScript "localfdr-versions" ''
              export PATH=${rEnv pkgs}/bin:${pyEnv pkgs}/bin:$PATH
              export LOCALFDR_BLAS=${blasLibrary pkgs}
              cd ${self}
              exec bash checks/versions.sh "$@"
            ''}";
            meta.description = "Print the software versions";
          };
          test = {
            type = "app";
            program = "${pkgs.writeShellScript "localfdr-test" ''
              set -euo pipefail
              ${unitTests pkgs}
              unit_r
              unit_rcpp
              unit_python "$@"
            ''}";
            meta.description = "Run the unit tests of the R, the Rcpp and the Python code in the current directory";
          };
          # `nix run .#help`, `nix run . -- --help` and `nix run .` print the
          # list of commands.
          help = help pkgs;
          default = help pkgs;
        }
      );

      # Build every table and figure with the default settings and compare
      # each with the paper; the result is the comparison (checks/RESULTS.md).
      checks = forEachSystem (pkgs: {
        outputs =
          pkgs.runCommand "localfdr-outputs"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                pkgs.python3
              ];
            }
            ''
              export HOME=$TMPDIR LOCALFDR_OUTPUT=$TMPDIR OPENBLAS_NUM_THREADS=1
              cd ${self}
              bash checks/run_checks.sh $TMPDIR/RESULTS.md || { cat $TMPDIR/RESULTS.md; exit 1; }
              cp $TMPDIR/RESULTS.md $out
            '';

        # The unit tests of every function in R/ (testthat, tests/testthat/).
        unit-r = pkgs.runCommand "localfdr-unit-r" { } ''
          export HOME=$TMPDIR
          ${unitTests pkgs}
          cd ${self}
          unit_r | tee $TMPDIR/unit-r.txt
          cp $TMPDIR/unit-r.txt $out
        '';

        # The unit tests of the Rcpp version (testthat, rcpp/tests/testthat/):
        # every C++ function of rcpp/src/ against its R function.
        unit-rcpp = pkgs.runCommand "localfdr-unit-rcpp" { } ''
          export HOME=$TMPDIR LOCALFDR_RCPP_CACHE=$TMPDIR/compiled
          ${unitTests pkgs}
          cd ${self}
          unit_rcpp | tee $TMPDIR/unit-rcpp.txt
          cp $TMPDIR/unit-rcpp.txt $out
        '';

        # The unit tests of every function in python/localfdr/ (pytest,
        # python/tests/), with the comparisons with R (marker against_r).
        unit-python = pkgs.runCommand "localfdr-unit-python" { } ''
          export HOME=$TMPDIR NUMBA_CACHE_DIR=$TMPDIR
          ${unitTests pkgs}
          cd ${self}
          unit_python | tee $TMPDIR/unit-python.txt
          cp $TMPDIR/unit-python.txt $out
        '';

        # Build every table and figure with the R and with the Python
        # scripts and compare the two; the result is the comparison
        # (checks/PYTHON.md, with the stored comparison of the reruns). Also
        # compare the software versions with those listed in README.md.
        python =
          pkgs.runCommand "localfdr-python"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                (pyEnv pkgs)
              ];
              LOCALFDR_BLAS = blasLibrary pkgs;
            }
            ''
              export HOME=$TMPDIR NUMBA_CACHE_DIR=$TMPDIR PYTHONDONTWRITEBYTECODE=1
              export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 CORES=4
              cd ${self}
              bash checks/check_python.sh archive $TMPDIR/work $TMPDIR/PYTHON.md $TMPDIR \
                || { cat $TMPDIR/PYTHON.md; exit 1; }
              bash checks/versions.sh --check
              cp $TMPDIR/PYTHON.md $out
            '';

        # Build every table and figure with the R and with the Rcpp scripts
        # and compare the two; the result is the comparison (checks/RCPP.md,
        # with the stored comparison of the reruns).
        rcpp =
          pkgs.runCommand "localfdr-rcpp"
            {
              nativeBuildInputs = [
                (rEnv pkgs)
                pkgs.python3
              ];
            }
            ''
              export HOME=$TMPDIR OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 CORES=4
              cd ${self}
              bash checks/check_rcpp.sh archive $TMPDIR/work $TMPDIR/RCPP.md $TMPDIR \
                || { cat $TMPDIR/RCPP.md; exit 1; }
              cp $TMPDIR/RCPP.md $out
            '';
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
