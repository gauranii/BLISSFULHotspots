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
      # (archive/reduce_archive.R).
      rEnv =
        pkgs:
        pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            ggplot2
            openxlsx
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
      });

      formatter = forEachSystem (pkgs: pkgs.nixfmt-rfc-style);
    };
}
