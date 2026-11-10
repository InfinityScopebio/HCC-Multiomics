# Third-party notices

The MIT License in `LICENSE` applies to repository-original code and documents.
It does not replace the licenses of upstream software or fixed upstream sources.

## CoVarNet

`code/pipelines/20_scrna_cm02/covarnet_v030_vendor.R` contains vendored
CoVarNet v0.3.0 functions, with the upstream commit recorded as
`d4c52a4b2ea821dc4b75861c3772e21232d1ee51`. The vendored file preserves the
historical frequency, correlation, filtering, and network behavior used by the
analysis. CoVarNet remains subject to its own copyright and license terms.

## Fixed analysis sources

`code/vendor/` contains fixed analysis sources retained for provenance and
numerical compatibility. Inclusion does not imply that those sources are
covered by the repository-original MIT grant. Their original copyright and
license notices remain applicable. Review the corresponding upstream terms
before redistribution.

## Scientific software and knowledge resources

The workflows call external R/Python packages and read separately supplied
MSigDB, PROGENy, CollecTRI, PTMsigDB, OmniPath, CellChat, DepMap/PRISM, TCGA,
GTEx, GEO, and other data resources. Those packages and datasets are not
relicensed by this repository. Their respective licenses, citation
requirements, and data-use conditions apply.
