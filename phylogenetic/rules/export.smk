"""
This part of the workflow collects the phylogenetic tree and annotations to
export a Nextstrain dataset.
REQUIRED INPUTS:
    metadata        = results/{serotype}/metadata.tsv
    tree            = results/{serotype}/{gene}/tree.nwk
    branch_lengths  = results/{serotype}/{gene}/branch-lengths.json
    node_data       = results/{serotype}/{gene}/*.json
OUTPUTS:
    auspice_json      = auspice/dengue_{serotype}_{gene}.json
    frequencies_json  = auspice/dengue_{serotype}_{gene}_tip-frequencies.json
    There are optional sidecar JSON files that can be exported as part of the dataset.
    See Nextstrain's data format docs for more details on sidecar files:
    https://docs.nextstrain.org/page/reference/data-formats.html
This part of the workflow usually includes the following steps:
    - augur export v2
    - augur frequencies
See Augur's usage docs for these commands for more details.
"""

import json
import shutil


rule colors:
    input:
        color_schemes = "defaults/color_schemes.tsv",
        color_orderings = "defaults/color_orderings.tsv",
        metadata = "results/{serotype}/metadata.tsv",
        manual_colors = "defaults/colors.tsv"
    output:
        colors = "results/{serotype}/colors.tsv"
    benchmark:
        "benchmarks/{serotype}/colors.txt"
    shell:
        """
        python3 scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}
        cat {input.manual_colors} >> {output.colors}
        """


rule prepare_auspice_config:
    """Prepare the auspice config file for each serotypes"""
    output:
        auspice_config="results/defaults/{serotype}/{gene}/auspice_config.json",
    benchmark:
        "benchmarks/{serotype}/{gene}/prepare_auspice_config.txt"
    params:
        replace_clade_key=lambda wildcard: r"clade_membership" if wildcard.gene in ['genome'] else r"major_lineage",
        replace_clade_title=lambda wildcard: r"Serotype" if wildcard.serotype in ['all'] else r"Genotype (Nextclade)",
    run:
        data = {
            "title": "Real-time tracking of dengue virus evolution",
            "maintainers": [
              {"name": "the Nextstrain team", "url": "https://nextstrain.org/team"}
            ],
            "data_provenance": [
              {
                "name": "GenBank",
                "url": "https://www.ncbi.nlm.nih.gov/genbank/"
              }
            ],
            "build_url": "https://github.com/nextstrain/dengue",
            "colorings": [
              {
                "key": "gt",
                "title": "Genotype",
                "type": "categorical"
              },
              {
                "key": "num_date",
                "title": "Date",
                "type": "continuous"
              },
              {
                "key": "country",
                "title": "Country",
                "type": "categorical"
              },
              {
                "key": "region",
                "title": "Region",
                "type": "categorical"
              },
              {
                "key": "deme",
                "title": "Deme (country in Asia, region elsewhere)",
                "type": "categorical"
              },
              {
                "key": "serotype_genbank",
                "title": "Serotype (Genbank metadata)",
                "type": "categorical"
              },
              {
                "key": "genotype",
                "title": "Genotype (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "major_lineage",
                "title": "Major lineage (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "minor_lineage",
                "title": "Minor lineage (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "host",
                "title": "Host Species",
                "type": "categorical"
              },
              {
                "key": "host_genus",
                "title": "Host Genus",
                "type": "categorical"
              },
              {
                "key": "host_type",
                "title": "Host Type",
                "type": "categorical"
              }
            ],
            "geo_resolutions": [
              "deme",
              "country",
              "region"
            ],
            "display_defaults": {
              "map_triplicate": True,
              "color_by": "deme",
              "geo_resolution": "deme",
              "transmission_lines": True,
              "tip_label": "strain"
            },
            "filters": [
              "deme",
              "country",
              "region",
              "author"
            ],
            "panels": [
              "tree",
              "map",
              "entropy",
              "frequencies"
            ],
            "metadata_columns": [
              "accession",
              "strain",
              "url"
            ]
          }

        # During genome/dengue_all workflows, clade membership represents Serotype
        # While genome/dengue_denvX workflows, clade_membership represents the more detailed Genotype
        if params.replace_clade_key == 'clade_membership':
            if wildcards.gene in ['genome'] and wildcards.serotype in ['all']:
                clade_membership_title="Serotype (Nextstrain)"
            else:
                clade_membership_title="Genotype (Nextstrain)"

            data["colorings"].append({
                "key": "clade_membership",
                "title": clade_membership_title,
                "type": "categorical"
            })
        else:
            # During E/dengue_all workflows, default color by Serotype
            if wildcards.serotype in ['all']:
                data["display_defaults"]["color_by"]="serotype_genbank"

        with open(output.auspice_config, 'w') as fh:
            json.dump(data, fh, indent=2)


rule export:
    """Exporting data files for auspice"""
    input:
        tree = "results/{serotype}/{gene}/tree.nwk",
        metadata = "results/{serotype}/metadata.tsv",
        branch_lengths = "results/{serotype}/{gene}/branch-lengths.json",
        traits = "results/{serotype}/{gene}/traits.json",
        clades = lambda wildcard: "results/{serotype}/{gene}/clades.json" if wildcard.gene in ['genome'] else [],
        nt_muts = "results/{serotype}/{gene}/nt-muts.json",
        aa_muts = "results/{serotype}/{gene}/aa-muts.json",
        description = config["export"]["description"],
        auspice_config = "results/defaults/{serotype}/{gene}/auspice_config.json",
        colors = "results/{serotype}/colors.tsv",
        lat_longs = "results/{serotype}/lat_longs_deme.tsv",
    output:
        auspice_json = "auspice/dengue_{serotype}_{gene}.json"
    benchmark:
        "benchmarks/{serotype}/{gene}/export.txt"
    params:
        strain_id = config.get("strain_id_field", "strain"),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --node-data {input.branch_lengths} {input.traits} {input.clades} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence-inline \
            --output {output.auspice_json}
        """

rule export_country_view:
    """
    Writing a second Auspice dataset defaulting to the standard
    country-everywhere phylogeography

    Both datasets come from the same tree and the same `augur traits` run, so
    this only swaps which geo-resolution Auspice opens on: the `deme` dataset
    opens on country-in-Asia/region-elsewhere demes, this one on countries
    worldwide. Either dataset can still be switched to the other resolution
    from the Auspice sidebar.
    """
    input:
        auspice_json = "auspice/dengue_{serotype}_{gene}.json",
        tip_freq = "auspice/dengue_{serotype}_{gene}_tip-frequencies.json",
    output:
        auspice_json = "auspice/dengue_{serotype}_{gene}_country.json",
        tip_freq = "auspice/dengue_{serotype}_{gene}_country_tip-frequencies.json",
    benchmark:
        "benchmarks/{serotype}/{gene}/export_country_view.txt"
    run:
        with open(input.auspice_json) as fh:
            data = json.load(fh)

        display_defaults = data.setdefault("meta", {}).setdefault("display_defaults", {})
        display_defaults["geo_resolution"] = "country"
        display_defaults["color_by"] = "country"
        display_defaults["transmission_lines"] = True

        with open(output.auspice_json, "w") as fh:
            json.dump(data, fh)

        # Auspice finds the sidecar by dataset name, so it needs a copy.
        shutil.copyfile(input.tip_freq, output.tip_freq)


rule tip_frequencies:
    """
    Estimating KDE frequencies for tips
    """
    input:
        tree = "results/{serotype}/{gene}/tree.nwk",
        metadata = "results/{serotype}/metadata.tsv",
    output:
        tip_freq = "auspice/dengue_{serotype}_{gene}_tip-frequencies.json"
    benchmark:
        "benchmarks/{serotype}/{gene}/tip_frequencies.txt"
    params:
        strain_id = config["strain_id_field"],
        min_date = config["tip_frequencies"]["min_date"],
        max_date = config["tip_frequencies"]["max_date"],
        narrow_bandwidth = config["tip_frequencies"]["narrow_bandwidth"],
        wide_bandwidth = config["tip_frequencies"]["wide_bandwidth"]
    shell:
        r"""
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --output {output.tip_freq}
        """
