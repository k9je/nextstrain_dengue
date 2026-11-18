"""
This part of the workflow defines the phylogeographic *demes* and their map
coordinates.

A deme is the geographic unit whose ancestral states `augur traits`
reconstructs, and hence the unit between which Auspice draws transmission
lines on the map. Here a deme is:

    * each COUNTRY inside Asia
    * each REGION outside Asia, with Nextstrain's coarse "North America"
      optionally split via defaults/deme_overrides.tsv

REQUIRED INPUTS:

    metadata = results/{serotype}/metadata_raw.tsv   (from merge_inputs.smk)

OUTPUTS:

    metadata  = results/{serotype}/metadata.tsv            (+ `deme` column)
    lat_longs = results/{serotype}/lat_longs_deme.tsv      (deme coordinates)

`add_deme` deliberately produces the canonical results/{serotype}/metadata.tsv
that every downstream rule already consumes, so filter, refine, traits, export
and tip_frequencies all see the `deme` column without further changes.
"""


rule add_deme:
    """
    Adding the `deme` column: country within {params.country_level_region}, region elsewhere
    """
    input:
        metadata = "results/{serotype}/metadata_raw.tsv",
        overrides = config["deme"]["overrides"],
    output:
        metadata = "results/{serotype}/metadata.tsv",
    params:
        country_level_region = config["deme"]["country_level_region"],
        deme_column = config["deme"]["column"],
    log:
        "logs/{serotype}/add_deme.txt",
    benchmark:
        "benchmarks/{serotype}/add_deme.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/add-deme.py \
            --metadata {input.metadata:q} \
            --overrides {input.overrides:q} \
            --country-level-region {params.country_level_region:q} \
            --deme-column {params.deme_column:q} \
            --output {output.metadata:q}
        """


rule deme_lat_longs:
    """
    Building map coordinates for every deme present in the metadata
    """
    input:
        metadata = "results/{serotype}/metadata.tsv",
        extra_lat_longs = config["deme"]["extra_lat_longs"],
    output:
        lat_longs = "results/{serotype}/lat_longs_deme.tsv",
    params:
        deme_column = config["deme"]["column"],
    log:
        "logs/{serotype}/deme_lat_longs.txt",
    benchmark:
        "benchmarks/{serotype}/deme_lat_longs.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/make-deme-lat-longs.py \
            --metadata {input.metadata:q} \
            --extra-lat-longs {input.extra_lat_longs:q} \
            --deme-column {params.deme_column:q} \
            --output {output.lat_longs:q}
        """
