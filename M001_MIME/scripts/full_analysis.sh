snakemake -j24 --configfile /home/moritz/projects/mosaic/M001_MIME/config_nd_metadata/M001_config.json --use-conda --until all_libs

python sourmash_nd_context.py

snakemake -j24 --configfile /home/moritz/projects/mosaic/M001_MIME/config_nd_metadata/M001_config.json --use-conda --until all_libs
