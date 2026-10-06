var roi = ee.FeatureCollection("projects/ideal-8895/assets/longyou/LY"),
    dn = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lydnfield200"),
    xn = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lyxnfield_200");


// -------------------- PARAMETERS --------------------
var REGION_NAME = 'LY';
var YEAR_A = 2024;
var YEAR_B = 2025;
var REF_MONTH = 5;

var CLASS_ON = 0;
var CLASS_OFF = 1;

var NDVI_WIN_MIN = 0.50;
var MIN_SLOPE_DEG = 5;
var BI_SPR_LO = 0.50;
var BI_SPR_HI = 0.85;

var LOW_PERCENTILE = 20;
var HIGH_PERCENTILE = 80;
var NABAI_EPS = 0.01;

var PURITY_BANDS = ['ndvi_spr','lswi_spr','re3_spr','nir_spr','re4_spr'];
var PURITY_KERNEL_RADIUS_PX = 1;
var IMPURITY_PERCENTILE = 30;
var IMPURITY_FALLBACK = 0.50;
var CORE_ERODE_RADIUS_PX = 1;

var N_ON = 200;
var N_OFF = 200;
var SAMPLE_SEED_ON = 101;
var SAMPLE_SEED_OFF = 202;
var SAMPLE_SPACING_M = 100;
var SAMPLE_OVERSAMPLE_FACTOR = 10;
var SAMPLE_TILE_SCALE = 4;

var EXPORT_FOLDER = 'POMM_FINAL_2026';
var region = roi.geometry();
Map.centerObject(region, 9);

// ============================================================
// 1. SENTINEL-2
// ============================================================

function maskS2(img){
  var scl = img.select('SCL');
  var clear = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
    .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(clear).copyProperties(img, ['system:time_start']);
}

function prepS2(img){
  var sr10 = img.select(['B2','B3','B4','B8']).divide(10000)
    .rename(['B2_sr','B3_sr','B4_sr','B8_sr']);
  var sr20 = img.select(['B5','B6','B7','B8A','B11','B12']).divide(10000)
    .resample('bilinear')
    .rename(['B5_sr','B6_sr','B7_sr','B8A_sr','B11_sr','B12_sr']);

  img = img.addBands(sr10).addBands(sr20);

  var b3 = img.select('B3_sr');
  var b4 = img.select('B4_sr');
  var b8 = img.select('B8_sr');
  var b11 = img.select('B11_sr');

  var ndvi = b8.subtract(b4).divide(b8.add(b4).add(1e-6)).rename('ndvi');
  var ndwi = b3.subtract(b8).divide(b3.add(b8).add(1e-6)).rename('ndwi');
  var lswi = b8.subtract(b11).divide(b8.add(b11).add(1e-6)).rename('lswi');
  var bi = ndvi.subtract(lswi).divide(ndvi.add(lswi).add(1e-6)).rename('BI');

  return img.addBands([ndvi, ndwi, lswi, bi]);
}

var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(region)
  .filterDate(ee.Date.fromYMD(YEAR_A,4,1), ee.Date.fromYMD(YEAR_B,REF_MONTH,1).advance(1,'month'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL'])
  .map(maskS2)
  .map(prepS2);

var p10 = ee.Image(s2.first()).select('B2_sr').projection();

function composite(start,end){
  return s2.filterDate(start,end).median().clip(region).setDefaultProjection(p10);
}

function monthComposite(year,month){
  var start = ee.Date.fromYMD(year,month,1);
  return composite(start,start.advance(1,'month'));
}

var spring = composite(ee.Date.fromYMD(YEAR_A,4,1), ee.Date.fromYMD(YEAR_A,6,1));
var winter = composite(ee.Date.fromYMD(YEAR_A,12,1), ee.Date.fromYMD(YEAR_A+1,3,1));

var biSpr = s2.filterDate(ee.Date.fromYMD(YEAR_A,4,1), ee.Date.fromYMD(YEAR_A,6,1))
  .select('BI').median().clip(region).rename('bi_spr').setDefaultProjection(p10);

var monthA = monthComposite(YEAR_A, REF_MONTH);
var monthB = monthComposite(YEAR_B, REF_MONTH);

// ============================================================
// 2. MAPBASE
// ============================================================

var dem = ee.Image('USGS/SRTMGL1_003');
var slope = ee.Terrain.slope(dem);

var treeMask = ee.Image('ESA/WorldCover/v200/2021').select('Map').eq(10).selfMask();

var auxMask = ee.Image.constant(1).clip(region)
  .updateMask(treeMask)
  .updateMask(winter.select('ndvi').gte(NDVI_WIN_MIN))
  .updateMask(winter.select('ndwi').lt(0).and(winter.select('lswi').gt(0)))
  .updateMask(slope.gte(MIN_SLOPE_DEG))
  .selfMask().setDefaultProjection(p10);

var biMask = biSpr.gte(BI_SPR_LO).and(biSpr.lte(BI_SPR_HI))
  .selfMask().setDefaultProjection(p10);

var mapbase = ee.Image.constant(1).clip(region).setDefaultProjection(p10)
  .updateMask(auxMask).updateMask(biMask)
  .rename('Mapbase').selfMask();

// ============================================================
// 3. NABAI — P20 / P80
// ============================================================

function ycbiBase(img){
  return img.select('B8A_sr').add(img.select('B8_sr')).add(img.select('B7_sr'));
}

var ycbiA = ycbiBase(monthA);
var ycbiB = ycbiBase(monthB);

var nabai = ycbiB.subtract(ycbiA)
  .divide(ycbiB.add(ycbiA).add(NABAI_EPS))
  .rename('NABAI').updateMask(mapbase).setDefaultProjection(p10);

var nabaiStats = nabai.reduceRegion({
  reducer: ee.Reducer.percentile([LOW_PERCENTILE,HIGH_PERCENTILE]),
  geometry: region,
  scale: 10,
  crs: p10,
  tileScale: 8,
  maxPixels: 1e13
});

var lowThreshold = ee.Number(nabaiStats.get('NABAI_p' + LOW_PERCENTILE));
var highThreshold = ee.Number(nabaiStats.get('NABAI_p' + HIGH_PERCENTILE));

print('NABAI P20 / P80:', lowThreshold, highThreshold);

var rawOn = nabai.gte(highThreshold).selfMask();
var rawOff = nabai.lte(lowThreshold).selfMask();

// ============================================================
// 4. P30 SPECTRAL PURIFICATION
// ============================================================

var purityImage = spring.select(
  ['ndvi','lswi','B7_sr','B8_sr','B8A_sr'],
  PURITY_BANDS
).updateMask(mapbase).setDefaultProjection(p10);

function zScore(img,bands,mask){
  var selected = img.select(bands).updateMask(mask);

  var stats = selected.reduceRegion({
    reducer: ee.Reducer.mean().combine({
      reducer2: ee.Reducer.stdDev(),
      sharedInputs: true
    }),
    geometry: region,
    scale: 30,
    bestEffort: true,
    tileScale: 8,
    maxPixels: 1e13
  });

  var images = ee.List(bands).map(function(b){
    b = ee.String(b);
    var meanKey = b.cat('_mean');
    var sdKey = b.cat('_stdDev');
    var mean = ee.Number(ee.Algorithms.If(stats.contains(meanKey), stats.get(meanKey), 0));
    var sd = ee.Number(ee.Algorithms.If(stats.contains(sdKey), stats.get(sdKey), 1));
    sd = ee.Number(ee.Algorithms.If(sd.gt(1e-6), sd, 1));
    return selected.select(b).subtract(mean).divide(sd).rename(b);
  });

  return ee.ImageCollection.fromImages(images).toBands()
    .rename(bands).setDefaultProjection(p10);
}

function localImpurity(img){
  var kernel = ee.Kernel.square({
    radius: PURITY_KERNEL_RADIUS_PX,
    units: 'pixels',
    normalize: false
  });

  var localMean = img.reduceNeighborhood({
    reducer: ee.Reducer.mean(),
    kernel: kernel,
    inputWeight: 'kernel',
    skipMasked: true
  }).rename(img.bandNames());

  return img.subtract(localMean).pow(2)
    .reduce(ee.Reducer.sum()).sqrt().rename('impurity');
}

var impurity = localImpurity(zScore(purityImage,PURITY_BANDS,mapbase))
  .updateMask(mapbase);

function impurityThreshold(mask,name){
  var stats = impurity.updateMask(mask).rename(name).reduceRegion({
    reducer: ee.Reducer.percentile([IMPURITY_PERCENTILE]),
    geometry: region,
    scale: 30,
    bestEffort: true,
    tileScale: 8,
    maxPixels: 1e13
  });

  var key = name + '_p' + IMPURITY_PERCENTILE;
  return ee.Number(ee.Algorithms.If(
    stats.contains(key),
    stats.get(key),
    IMPURITY_FALLBACK
  ));
}

var impurityOnThreshold = impurityThreshold(rawOn,'on');
var impurityOffThreshold = impurityThreshold(rawOff,'off');

print('P30 impurity thresholds On / Off:', impurityOnThreshold, impurityOffThreshold);

var pureOn = rawOn.updateMask(impurity.lte(impurityOnThreshold)).selfMask();
var pureOff = rawOff.updateMask(impurity.lte(impurityOffThreshold)).selfMask();

// ============================================================
// 5. 3×3 SPATIAL EROSION
// No connected-component / area threshold
// ============================================================

var coreOn = pureOn.unmask(0).focalMin({
  radius: CORE_ERODE_RADIUS_PX,
  units: 'pixels',
  kernelType: 'square',
  iterations: 1
}).updateMask(mapbase).rename('core_on').selfMask();

var coreOff = pureOff.unmask(0).focalMin({
  radius: CORE_ERODE_RADIUS_PX,
  units: 'pixels',
  kernelType: 'square',
  iterations: 1
}).updateMask(mapbase).rename('core_off').selfMask();

// ============================================================
// 6. CORE AREA
// ============================================================

function areaHa(mask){
  var areaImage = ee.Image.pixelArea().updateMask(mask.selfMask()).rename('area');

  var stats = areaImage.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: region,
    scale: 10,
    crs: p10,
    tileScale: 8,
    maxPixels: 1e13
  });

  return ee.Number(ee.Algorithms.If(
    stats.contains('area'),
    stats.get('area'),
    0
  )).divide(10000);
}

var coreOnAreaHa = areaHa(coreOn);
var coreOffAreaHa = areaHa(coreOff);

print('Core On area (ha):', coreOnAreaHa);
print('Core Off area (ha):', coreOffAreaHa);

// ============================================================
// 7. AUTOMATIC SAMPLING
// ============================================================

function autoSample(mask,classValue,n,seed){
  var cls = ee.Image.constant(classValue).rename('class').toInt16();
  var xy = ee.Image.pixelCoordinates(p10).rename(['px','py']).toInt64();

  return cls.addBands(xy).updateMask(mask).stratifiedSample({
    numPoints: n,
    classBand: 'class',
    region: region,
    scale: 10,
    projection: p10,
    seed: seed,
    dropNulls: true,
    tileScale: 4,
    geometries: true
  }).map(function(f){
    return ee.Feature(f.geometry()).set({
      class: classValue,
      source: 'auto',
      region: REGION_NAME,
      label_year: YEAR_A,
      year_a: YEAR_A,
      year_b: YEAR_B
    });
  });
}

var generatedOn = autoSample(coreOn,CLASS_ON,N_ON,SAMPLE_SEED_ON);
var generatedOff = autoSample(coreOff,CLASS_OFF,N_OFF,SAMPLE_SEED_OFF);
var generatedSamples = generatedOn.merge(generatedOff);

print('Generated On:', generatedOn.size());
print('Generated Off:', generatedOff.size());
print('Generated total:', generatedSamples.size());

// ============================================================
// 8. FIELD SAMPLES — ONLY FOR SAD / ED
// ============================================================

function centroidClass(fc,classValue){
  return ee.FeatureCollection(fc).map(function(f){
    return ee.Feature(f.geometry().centroid(1)).set({
      class: classValue,
      source: 'field'
    });
  });
}

var fieldOn = centroidClass(dn,CLASS_ON);
var fieldOff = centroidClass(xn,CLASS_OFF);

print('Field On / Off:', fieldOn.size(), fieldOff.size());

// ============================================================
// 9. SAD / ED — CLASS-WISE MEAN SPECTRAL VECTORS
// ============================================================

// ============================================================
// 9. SAD / ED — CLASS-WISE MEAN SPECTRAL VECTORS
// ============================================================

var featureImage = spring.select(
  ['B4_sr','ndvi','lswi','B7_sr','B8_sr','B8A_sr','B11_sr'],
  ['red','ndvi_spr','lswi_spr','re3_spr','nir_spr','re4_spr','swir1_spr']
).setDefaultProjection(p10);

var featureBands = featureImage.bandNames();

function meanVector(fc,img,bands){
  var sampled = img.select(bands).sampleRegions({
    collection: fc,
    scale: 10,
    projection: p10,
    geometries: false,
    tileScale: 4
  });

  print('Valid feature samples:', sampled.size());

  return ee.Array(bands.map(function(b){
    b = ee.String(b);
    var stats = sampled.reduceColumns(ee.Reducer.mean(), [b]);
    var value = stats.get('mean');
    return ee.Number(ee.Algorithms.If(value, value, 0));
  }));
}

function cosineSimilarity(a,b){
  var dot = a.multiply(b).reduce('sum',[0]);
  var normA = a.pow(2).reduce('sum',[0]).sqrt();
  var normB = b.pow(2).reduce('sum',[0]).sqrt();
  return dot.divide(normA.multiply(normB));
}

function spectralAngle(a,b){
  return cosineSimilarity(a,b).max(-1).min(1).acos();
}

function euclideanDistance(a,b){
  return a.subtract(b).pow(2).reduce('sum',[0]).sqrt();
}

var vectorFieldOn = meanVector(fieldOn,featureImage,featureBands);
var vectorGeneratedOn = meanVector(generatedOn,featureImage,featureBands);
var vectorFieldOff = meanVector(fieldOff,featureImage,featureBands);
var vectorGeneratedOff = meanVector(generatedOff,featureImage,featureBands);

var onCosine = cosineSimilarity(vectorFieldOn,vectorGeneratedOn);
var onSAD = spectralAngle(vectorFieldOn,vectorGeneratedOn);
var onED = euclideanDistance(vectorFieldOn,vectorGeneratedOn);

var offCosine = cosineSimilarity(vectorFieldOff,vectorGeneratedOff);
var offSAD = spectralAngle(vectorFieldOff,vectorGeneratedOff);
var offED = euclideanDistance(vectorFieldOff,vectorGeneratedOff);

print('ON field vs generated - cosine / SAD(rad) / ED:',
  onCosine,onSAD,onED);

print('OFF field vs generated - cosine / SAD(rad) / ED:',
  offCosine,offSAD,offED);

// ============================================================
// 10. SUMMARY
// ============================================================

var summary = ee.Feature(null,{
  region: REGION_NAME,
  year_a: YEAR_A,
  year_b: YEAR_B,
  ref_month: REF_MONTH,

  ndvi_win_min: NDVI_WIN_MIN,
  slope_min_deg: MIN_SLOPE_DEG,
  bi_low: BI_SPR_LO,
  bi_high: BI_SPR_HI,

  nabai_low_percentile: LOW_PERCENTILE,
  nabai_high_percentile: HIGH_PERCENTILE,
  nabai_low_threshold: lowThreshold,
  nabai_high_threshold: highThreshold,

  impurity_percentile: IMPURITY_PERCENTILE,
  impurity_on_threshold: impurityOnThreshold,
  impurity_off_threshold: impurityOffThreshold,

  erosion_kernel_px: CORE_ERODE_RADIUS_PX * 2 + 1,

  core_on_area_ha: coreOnAreaHa,
  core_off_area_ha: coreOffAreaHa,

  sample_spacing_m: SAMPLE_SPACING_M,
  requested_on_n: N_ON,
  requested_off_n: N_OFF,
  auto_on_n: generatedOn.size(),
  auto_off_n: generatedOff.size(),

field_on_n: fieldOn.size(),
auto_on_n_sad: generatedOn.size(),
field_off_n: fieldOff.size(),
auto_off_n_sad: generatedOff.size(),

on_cosine: onCosine,
on_SAD_rad: onSAD,
on_ED: onED,

off_cosine: offCosine,
off_SAD_rad: offSAD,
off_ED: offED
});

print('FINAL SAMPLE SUMMARY:', summary);

// ============================================================
// 11. VISUALIZATION
// ============================================================

Map.addLayer(mapbase, {palette:['ffff00']}, 'Mapbase', false);
Map.addLayer(nabai, {min:-0.15,max:0.15,palette:['0000ff','ffffff','ff0000']}, 'NABAI', false);
Map.addLayer(rawOn, {palette:['00ff88']}, 'Raw On P80', false);
Map.addLayer(rawOff, {palette:['ff5555']}, 'Raw Off P20', false);
Map.addLayer(impurity, {min:0,max:2}, 'Local impurity', false);
Map.addLayer(coreOn, {palette:['00ff88']}, 'Final On core', false);
Map.addLayer(coreOff, {palette:['ff5555']}, 'Final Off core', false);
Map.addLayer(generatedOn, {color:'00ff88'}, 'Generated On', true);
Map.addLayer(generatedOff, {color:'ff5555'}, 'Generated Off', true);
Map.addLayer(fieldOn, {color:'006400'}, 'Field On', false);
Map.addLayer(fieldOff, {color:'8b0000'}, 'Field Off', false);

// ============================================================
// 12. EXPORT
// ============================================================

Export.table.toDrive({
  collection: generatedOn,
  description: REGION_NAME + '_POMM_On_' + YEAR_A + '_N' + N_ON,
  folder: REGION_NAME +'sample_generation',
  fileNamePrefix: REGION_NAME + '_POMM_On_' + YEAR_A + '_N' + N_ON,
  fileFormat: 'SHP'
});

Export.table.toDrive({
  collection: generatedOff,
  description: REGION_NAME + '_POMM_Off_' + YEAR_A + '_N' + N_OFF,
  folder: REGION_NAME +'sample_generation',
  fileNamePrefix: REGION_NAME + '_POMM_Off_' + YEAR_A + '_N' + N_OFF,
  fileFormat: 'SHP'
});

Export.table.toDrive({
  collection: generatedSamples,
  description: REGION_NAME + '_POMM_AutoSamples_' + YEAR_A,
  folder: REGION_NAME +'sample_generation',
  fileNamePrefix: REGION_NAME + '_POMM_AutoSamples_' + YEAR_A,
  fileFormat: 'SHP'
});

Export.table.toDrive({
  collection: ee.FeatureCollection([summary]),
  description: REGION_NAME + '_POMM_SampleSummary',
  folder: REGION_NAME +'_POMM_SampleSummary',
  fileNamePrefix: REGION_NAME + '_POMM_SampleSummary',
  fileFormat: 'CSV'
});