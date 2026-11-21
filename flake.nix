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
      # (archive/reduce_archive.R), testthat runs the unit tests.
      rEnv =
        pkgs:
        pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            ggplot2
            openxlsx
            testthat
          ];
        };

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

      # The unit tests: testthat over tests/testthat/, started in the current
      # directory.
      # PYTHON: the tests of the Python version are added here (a function
      # unit_python next to unit_r, and its interpreter on PATH).
      unitTests = pkgs: ''
        export OPENBLAS_NUM_THREADS=1
        export PATH=${rEnv pkgs}/bin:$PATH
        unit_r() { Rscript tests/testthat.R; }
      '';

      # The text `nix run .#help` prints: every command of the flake.
      helpText =
        let
          row =
            command: text:
            "  nix run .#${command}${lib.fixedWidthString (12 - builtins.stringLength command) " " ""}${text}";
          rows = lib.mapAttrsToList (name: script: row name "scripts/${script}.R") scripts;
        in
        lib.concatStringsSep "\n" (
          [
            "BayesLocalFDR: code of Gauran, Park, Rattsev, Peterson, Kann, Park (2022), Ann. Appl. Stat. 16(3):1459-1475"
            ""
            "One command per table and figure of the paper (R):"
          ]
          ++ rows
          ++ [
            (row "all" "every table and figure")
            ""
            # PYTHON: the commands of the Python version are listed here.
            "Other commands:"
            (row "test" "unit tests of R/ in the current directory")
            (row "versions" "print the software versions")
            (row "help" "print this list (also: nix run . -- --help)")
            "  nix flake check       build every output, compare it with the paper, run the unit tests"
            "  nix develop           a shell with the same R and Python"
            ""
            "Results are written to results/ and figures/ of the current directory (or of LOCALFDR_OUTPUT)."
            "Every table and figure command takes the options listed in the header of its script, after --."
            "The simulation studies (table1, table2, figS1, figS3) read the archived runs unless --rerun is"
            "given, for example"
            "  nix run .#table1 -- --rerun --runs=1:10 --cores=4"
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
            pkgs.python3
          ];
          OPENBLAS_NUM_THREADS = 1;
        };
      });

      apps = forEachSystem (
        pkgs:
        lib.mapAttrs (name: script: app pkgs name [ script ]) scripts
        // {
          all = app pkgs "all" (lib.attrValues scripts);
          versions = {
            type = "app";
            program = "${pkgs.writeShellScript "localfdr-versions" ''
              export PATH=${rEnv pkgs}/bin:${pkgs.python3}/bin:$PATH
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
              # PYTHON: unit_python "$@" is added here.
            ''}";
            meta.description = "Run the unit tests of the R code in the current directory";
          };
          # `nix run .#help`, `nix run . -- --help` and `nix run .` print the
          # list of commands.
          help = help pkgs;
          default = help pkgs;
        }
      );

      # Build every table and figure with the default settings and compare
      # each with the paper; the result is the comparison (checks/RESULTS.md).
      # Also compare the software versions with those listed in README.md.
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
              bash checks/versions.sh --check
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

        # PYTHON: the checks of the Python version (its unit tests and the
        # comparison of its outputs with those of R) are added here.
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
