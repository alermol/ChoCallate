# Feature Roadmap

This document outlines the future direction of our project. It's a living document and will be updated as our priorities evolve.

## Portability

- [ ] Docker image
- [ ] Singularity image

## Ployploid calling

- [ ] Integration of [Octopus](https://github.com/luntergroup/octopus) into the pipeline

## Mapping and calling using long reads

- [ ] Integration of [Clair3](https://github.com/HKU-BAL/Clair3) into the pipeline
    - [ ] Train models for plants
- [ ] Integration of [longshot](https://github.com/pjedge/longshot) into the pipeline
- [ ] Integration of [DeepVariant](https://github.com/google/deepvariant) into the pipeline

## Improvement of consensus file

- [ ] Expandind set of tags in consensus file based on tags from individual callers

## General improvements of calling

- [ ] Integration of [longdust](https://github.com/lh3/longdust) in order to exclude "bad" regions
- [ ] Integration of "exclude_bed" and "include_bed" region
- [ ] Integration of [adaptive filtering of variant calls](https://github.com/yyren/FVC) into the pipeline