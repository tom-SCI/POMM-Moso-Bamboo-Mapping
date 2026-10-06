
# POMM: Phenology-driven On-/Off-year Moso Bamboo Forest Mapping

Google Earth Engine source code for automatic training-sample generation and cross-regional mapping of on-/off-year Moso bamboo forests.

## Scripts

- `01_POMM_Automatic_Sample_Generation.js`  
  Generates automatic On-/Off-year samples using Mapbase, NABAI, regional P20/P80 thresholds, spectral purification, and spatial erosion.

- `02_POMM_Workflow_Reapplication.js`  
  Implements workflow re-application in a target region using locally generated automatic samples and Random Forest classification.

- `03_POMM_Model_Transfer.js`  
  Implements direct transfer of the Deqing-trained Random Forest model to a target region.

## Data

The scripts use:

- Sentinel-2 SR Harmonized
- AlphaEarth Foundations (AEF)
- ESA WorldCover
- SRTM DEM
- User-provided ROI and sample assets

Class labels:

```text
0 = On-year
1 = Off-year
2 = Other / non-MBF
```

## Usage

Run the scripts in the Google Earth Engine JavaScript Code Editor and replace the user-specific asset paths with your own GEE assets.

For workflow re-application, run `01_POMM_Automatic_Sample_Generation.js` first to generate local automatic samples, then use them in `02_POMM_Workflow_Reapplication.js`.

## Citation

Li, X., Li, L., Li, N., et al.  
*Mapping on- and off-year Moso bamboo forests across regions with automatically generated training samples.*  
Remote Sensing of Environment.

Citation information will be updated after publication.
