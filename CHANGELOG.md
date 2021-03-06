# Changelog

All notable changes to this project will be documented in this file.

## [2.0.10] - 2019-07-29
### Added
- Add genome database for intergenic region (hg19, hg38)
- Add functions for modeling intergenic region
- Add warning message if there is no variant past the feature selection
- Add warning message if there is no variant remained in previous step
### Fixed
- Update the UCSC mysql host
- Fix error when modeling empty feature set (step5)
- Fix error when modeling empty feature set (step6)
- Refine PRS plot

## [2.0.9] - 2019-12-31
### Added
- Add standalone AppGenEpi.app
### Fixed
- Fix the bug of joblib multithreading

## [2.0.8] - 2019-12-26
### Fixed
- Fix unsupported image format for ploting PRS

## [2.0.7] - 2019-12-26
### Fixed
- Fix executable python for AppGenEpi 

## [2.0.6] - 2019-12-25
### Fixed
- Fix file path bug for ploting PRS 
### Removed
- Remove unused PIL import

## [2.0.5] - 2019-12-24
### Added
- Add AppGenEpi (GUI)
- Add polygenic risk score calculation for classifier

## [2.0.4] - 2019-12-12
### Added
- Add sliding windows scanning to deal with mega genes

## [2.0.3] - 2019-10-13
### Added
- Add GenEpi's documentation to Read the Docs
- Add model persistance for classifier and regressor

## [2.0.2] - 2019-09-18
### Added
- Add multiprocessing for gene batch running

## [2.0.1] - 2019-07-12
### Added
- Support self-defiend genome regions (for any species)
- Add argument parser (GenEpi -h)
- Add "GenEpi" to bin when installing
- Add output log
- Add quick test
- Add change log in Changelog.md

### Changed
- Change feature encoder for reducing memory usage

### Removed
- Remove random forest ensemble
