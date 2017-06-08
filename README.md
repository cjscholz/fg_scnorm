# fg_scnorm
This docker container is a wrapper around the SCnorm R-package (https://github.com/rhondabacher/SCnorm) to use the algorithm in the FAST Genomics pipeline

## Analysis Configuration
Analysis parameters for SCnorm can be adusted by making changes to the './config/SCnorm_config.yml' file.

## Data Input/Output
The input file directory needs to be mounted as '/fastgenomics/input', the output file directory needs to be mounted as '/fastgenomics/output'.
