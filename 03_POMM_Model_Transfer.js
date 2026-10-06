var roiDQ = ee.FeatureCollection("projects/ideal-8895/assets/DQ/deqing"),
    dqGenOnRaw = ee.FeatureCollection("projects/ideal-8895/assets/DQ/dqdngen_200"),
    dqGenOffRaw = ee.FeatureCollection("projects/ideal-8895/assets/DQ/dqxngen_200"),
    dqOtherRaw = ee.FeatureCollection("projects/ideal-8895/assets/DQ/qt"),
    targetOnRaw = ee.FeatureCollection("projects/ideal-8895/assets/anji/ajdnfield200"),
    targetOffRaw = ee.FeatureCollection("projects/ideal-8895/assets/anji/ajxnfield200"),
    targetOtherRaw = ee.FeatureCollection("projects/ideal-8895/assets/anji/qtfield200"),
    roiTarget = ee.FeatureCollection("projects/ideal-8895/assets/anji/AJ");

// -------------------- PARAMETERS --------------------
var TARGET_NAME = 'AJ';

var YEAR_A = 2024;
var YEAR_B = 2025;
var REF_MONTH = 5;

var N_ON = 200;
var N_OFF = 200;
var N_OTHER_TRAIN = 200;

var SAMPLE_SEED_ON = 101;
var SAMPLE_SEED_OFF = 202;
var OTHER_TRAIN_SEED = 303;

var N_FIELD_PER_CLASS = 200;
var FIELD_SEED_ON = 2024;
var FIELD_SEED_OFF = 2025;
var FIELD_SEED_OTHER = 2026;

var RF_SEED = 42;
var N_TREES = 200;
var BAG_FRACTION = 0.7;

var NDVI_WIN_MIN = 0.50;
var EXPORT_FOLDER = 'POMM_TRANSFER_2026';

var geomDQ = roiDQ.geometry();
var geomTarget = roiTarget.geometry();

// ============================================================
// 1. SAMPLE PREPARATION
// ============================================================

function centroidClass(fc,classValue){
  return ee.FeatureCollection(fc).map(function(f){
    return ee.Feature(f.geometry().centroid(1)).copyProperties(f).set('class',classValue);
  });
}

function limitFC(fc,n,seed,columnName){
  return ee.FeatureCollection(fc).randomColumn(columnName,seed).sort(columnName).limit(n);
}

var dqOn = limitFC(centroidClass(dqGenOnRaw,0),N_ON,SAMPLE_SEED_ON,'dq_on_random');
var dqOff = limitFC(centroidClass(dqGenOffRaw,1),N_OFF,SAMPLE_SEED_OFF,'dq_off_random');
var dqOther = limitFC(centroidClass(dqOtherRaw,2),N_OTHER_TRAIN,OTHER_TRAIN_SEED,'dq_other_random');
var dqTrainPts = dqOn.merge(dqOff).merge(dqOther);

var targetOn = limitFC(centroidClass(targetOnRaw,0),N_FIELD_PER_CLASS,FIELD_SEED_ON,'field_on_random');
var targetOff = limitFC(centroidClass(targetOffRaw,1),N_FIELD_PER_CLASS,FIELD_SEED_OFF,'field_off_random');
var targetOther = limitFC(centroidClass(targetOtherRaw,2),N_FIELD_PER_CLASS,FIELD_SEED_OTHER,'field_other_random');
var targetTestPts = targetOn.merge(targetOff).merge(targetOther);

print('DQ train On / Off / Other:',dqOn.size(),dqOff.size(),dqOther.size());
print(TARGET_NAME+' test On / Off / Other:',targetOn.size(),targetOff.size(),targetOther.size());

// ============================================================
// 2. SENTINEL-2 FEATURES
// May 2024 + May 2025
// ============================================================

function maskS2(img){
  var scl = img.select('SCL');
  var clear = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
    .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(clear).copyProperties(img,['system:time_start']);
}

function prepS2(img){
  var sr10 = img.select(['B2','B3','B4','B8']).divide(10000)
    .rename(['B2_sr','B3_sr','B4_sr','B8_sr']);
  var sr20 = img.select(['B5','B6','B7','B8A','B11','B12']).divide(10000)
    .resample('bilinear')
    .rename(['B5_sr','B6_sr','B7_sr','B8A_sr','B11_sr','B12_sr']);
  return img.addBands(sr10).addBands(sr20);
}

function prefixBands(img,prefix){
  return img.rename(img.bandNames().map(function(name){
    return ee.String(prefix).cat('_').cat(ee.String(name));
  }));
}

function buildS2Data(geometry){
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(geometry)
    .filterDate(ee.Date.fromYMD(YEAR_A,4,1),ee.Date.fromYMD(YEAR_B,REF_MONTH,1).advance(1,'month'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
    .select(['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12','SCL'])
    .map(maskS2)
    .map(prepS2);

  var p10 = ee.Image(col.first()).select('B2_sr').projection();

  function monthComposite(year){
    var start = ee.Date.fromYMD(year,REF_MONTH,1);
    return col.filterDate(start,start.advance(1,'month')).median()
      .clip(geometry).setDefaultProjection(p10);
  }

  var bands = ['B4_sr','B5_sr','B6_sr','B7_sr','B8A_sr','B11_sr','B12_sr'];
  //var bands = ['B12_sr'];
  var monthA = monthComposite(YEAR_A);
  var monthB = monthComposite(YEAR_B);

  var image = prefixBands(monthA.select(bands),'S2A')
    .addBands(prefixBands(monthB.select(bands),'S2B'))
    .toFloat();

  return {image:image,projection:p10};
}

// ============================================================
// 3. AEF FEATURES
// 2024 + 2025 + dot + angle
// ============================================================

function getAEF(year,geometry){
  var start = ee.Date.fromYMD(year,1,1);
  return ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
    .filterDate(start,start.advance(1,'year'))
    .filterBounds(geometry)
    .mosaic()
    .clip(geometry)
    .toFloat();
}

function buildAEF(geometry){
  var embY = getAEF(YEAR_A,geometry);
  var embY1 = getAEF(YEAR_B,geometry);

  var Y = prefixBands(embY,'Y');
  var Y1 = prefixBands(embY1,'Y1');

  var dot = embY.multiply(embY1).reduce(ee.Reducer.sum()).rename('ae_dot');
  var angle = dot.max(-1).min(1).acos().rename('ae_ang');

  return Y.addBands(Y1).addBands(dot).addBands(angle).toFloat();
}

// ============================================================
// 4. BUILD FEATURE IMAGES
// ============================================================

var s2DQ = buildS2Data(geomDQ);
var s2Target = buildS2Data(geomTarget);

var aefImgDQ = buildAEF(geomDQ);
var aefImgTarget = buildAEF(geomTarget);

var featImgDQ = aefImgDQ.addBands(s2DQ.image).toFloat();
var featImgTarget = aefImgTarget.addBands(s2Target.image).toFloat();

var featBands = featImgDQ.bandNames();
var aefBands = aefImgDQ.bandNames();

print('Full feature number:',featBands.length());
print('AEF feature number:',aefBands.length());

print('Full DQ / Target bands identical:',
  ee.Algorithms.IsEqual(featImgDQ.bandNames(),featImgTarget.bandNames()));

print('AEF DQ / Target bands identical:',
  ee.Algorithms.IsEqual(aefImgDQ.bandNames(),aefImgTarget.bandNames()));

// ============================================================
// 5. EXTRACT DEQING TRAINING FEATURES
// ============================================================

var dqTrainingFull = featImgDQ.select(featBands).sampleRegions({
  collection:dqTrainPts,
  properties:['class'],
  scale:10,
  projection:s2DQ.projection,
  geometries:false,
  tileScale:16
}).filter(ee.Filter.notNull(featBands));

var dqTrainingAEF = aefImgDQ.select(aefBands).sampleRegions({
  collection:dqTrainPts,
  properties:['class'],
  scale:10,
  projection:s2DQ.projection,
  geometries:false,
  tileScale:16
}).filter(ee.Filter.notNull(aefBands));

print('Effective DQ training n:',dqTrainingFull.size());
print('DQ training histogram:',dqTrainingFull.aggregate_histogram('class'));

// ============================================================
// 6. TRAIN INTEGRATED CLASSIFICATION SYSTEM
// ============================================================

var rfFull = ee.Classifier.smileRandomForest({
  numberOfTrees:N_TREES,
  variablesPerSplit:null,
  minLeafPopulation:1,
  bagFraction:BAG_FRACTION,
  maxNodes:null,
  seed:RF_SEED
}).train({
  features:dqTrainingFull,
  classProperty:'class',
  inputProperties:featBands
});

var rfAEF = ee.Classifier.smileRandomForest({
  numberOfTrees:N_TREES,
  variablesPerSplit:null,
  minLeafPopulation:1,
  bagFraction:BAG_FRACTION,
  maxNodes:null,
  seed:RF_SEED
}).train({
  features:dqTrainingAEF,
  classProperty:'class',
  inputProperties:aefBands
});

print('DQ integrated RF training completed');

// ============================================================
// 7. INTEGRATED TARGET CLASSIFICATION
// Full AEF+S2 where available; AEF fills missing S2 areas.
// ============================================================

var fullRawMap = featImgTarget.select(featBands)
  .classify(rfFull)
  .rename('classification');

var aefRawMap = aefImgTarget.select(aefBands)
  .classify(rfAEF)
  .rename('classification');

var targetRawMap = fullRawMap
  .unmask(aefRawMap)
  .rename('classification')
  .clip(geomTarget);

// ============================================================
// 8. INDEPENDENT TARGET TEST
// ============================================================

var targetPred = fullRawMap.sampleRegions({
  collection:targetTestPts,
  properties:['class'],
  scale:10,
  projection:s2Target.projection,
  geometries:true,
  tileScale:16
}).filter(ee.Filter.notNull(['classification']));

print('Effective target test n:',targetPred.size());
print('Target test histogram:',targetPred.aggregate_histogram('class'));

// ============================================================
// 9. ACCURACY METRICS
// ============================================================

function safeDiv(a,b){
  a = ee.Number(a);
  b = ee.Number(b);
  return ee.Number(ee.Algorithms.If(b.gt(0),a.divide(b),0));
}

function f1(p,r){
  return ee.Number(ee.Algorithms.If(
    p.add(r).gt(0),
    p.multiply(r).multiply(2).divide(p.add(r)),
    0
  ));
}

var cm = targetPred.errorMatrix( 'class','classification',[0,1,2]);
var a = cm.array();

var c00 = ee.Number(a.get([0,0]));
var c01 = ee.Number(a.get([0,1]));
var c02 = ee.Number(a.get([0,2]));
var c10 = ee.Number(a.get([1,0]));
var c11 = ee.Number(a.get([1,1]));
var c12 = ee.Number(a.get([1,2]));
var c20 = ee.Number(a.get([2,0]));
var c21 = ee.Number(a.get([2,1]));
var c22 = ee.Number(a.get([2,2]));

var PA0 = safeDiv(c00,c00.add(c01).add(c02));
var PA1 = safeDiv(c11,c10.add(c11).add(c12));
var PA2 = safeDiv(c22,c20.add(c21).add(c22));

var UA0 = safeDiv(c00,c00.add(c10).add(c20));
var UA1 = safeDiv(c11,c01.add(c11).add(c21));
var UA2 = safeDiv(c22,c02.add(c12).add(c22));

var F10 = f1(PA0,UA0);
var F11 = f1(PA1,UA1);
var F12 = f1(PA2,UA2);
var macroF1 = F10.add(F11).add(F12).divide(3);

print('================================================');
print('MODEL TRANSFER: DEQING -> '+TARGET_NAME);
print('================================================');
print('Confusion matrix:',cm);
print('OA:',cm.accuracy());
print('Kappa:',cm.kappa());
print('PA On / Off / Other:',PA0,PA1,PA2);
print('UA On / Off / Other:',UA0,UA1,UA2);
print('F1 On / Off / Other:',F10,F11,F12);
print('Macro-F1:',macroF1);

// ============================================================
// 10. PAPER-READY RESULT
// ============================================================

var result = ee.Feature(null,{
  source_region:'Deqing',
  target_region:TARGET_NAME,
  transfer_mode:'Model_transfer',
  feature_mode:'AEF_S2',

  year_a:YEAR_A,
  year_b:YEAR_B,
  ref_month:REF_MONTH,

  rf_trees:N_TREES,
  rf_bag_fraction:BAG_FRACTION,
  rf_seed:RF_SEED,

  train_n:dqTrainingFull.size(),
  test_n:targetPred.size(),

  train_On_n:dqTrainingFull.filter(ee.Filter.eq('class',0)).size(),
  train_Off_n:dqTrainingFull.filter(ee.Filter.eq('class',1)).size(),
  train_Other_n:dqTrainingFull.filter(ee.Filter.eq('class',2)).size(),

  test_On_n:targetPred.filter(ee.Filter.eq('class',0)).size(),
  test_Off_n:targetPred.filter(ee.Filter.eq('class',1)).size(),
  test_Other_n:targetPred.filter(ee.Filter.eq('class',2)).size(),

  OA:cm.accuracy().multiply(100),
  Kappa:cm.kappa(),

  PA_On:PA0.multiply(100),
  PA_Off:PA1.multiply(100),
  PA_Other:PA2.multiply(100),

  UA_On:UA0.multiply(100),
  UA_Off:UA1.multiply(100),
  UA_Other:UA2.multiply(100),

  F1_On:F10.multiply(100),
  F1_Off:F11.multiply(100),
  F1_Other:F12.multiply(100),
  Macro_F1:macroF1.multiply(100),

  CM_00:c00,CM_01:c01,CM_02:c02,
  CM_10:c10,CM_11:c11,CM_12:c12,
  CM_20:c20,CM_21:c21,CM_22:c22
});

print('FINAL MODEL TRANSFER RESULT:',result);

// ============================================================
// 11. BROAD FOREST DOMAIN MASK
// Only used for final spatial map.
// ============================================================

function maskS2Forest(img){
  var scl = img.select('SCL');
  var clear = scl.neq(0).and(scl.neq(1)).and(scl.neq(3))
    .and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  return img.updateMask(clear);
}

function addNDVI(img){
  return img.addBands(img.normalizedDifference(['B8','B4']).rename('ndvi'));
}

var treeMask = ee.Image('ESA/WorldCover/v200/2021').select('Map').eq(10);

var winterTarget = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(geomTarget)
  .filterDate(ee.Date.fromYMD(YEAR_A,12,1),ee.Date.fromYMD(YEAR_B,3,1))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',30))
  .select(['B4','B8','SCL'])
  .map(maskS2Forest)
  .map(addNDVI)
  .median()
  .clip(geomTarget);

var winterNdviMask = winterTarget.select('ndvi').gte(NDVI_WIN_MIN);

var forestDomainMask = treeMask
  .and(winterNdviMask)
  .selfMask()
  .clip(geomTarget);

var targetMap = targetRawMap
  .updateMask(forestDomainMask)
  .clip(geomTarget);

// ============================================================
// 12. MAP DISPLAY
// ============================================================

Map.centerObject(roiTarget,10);

Map.addLayer(forestDomainMask,{
  palette:['006400']
},'Broad Forest Domain',false);

Map.addLayer(targetMap,{
  min:0,
  max:2,
  palette:['006400','90EE90','D2B48C']
},'DQ -> '+TARGET_NAME+' Model Transfer',true);

// ============================================================
// 13. EXPORT RESULT
// ============================================================

Export.table.toDrive({
  collection:ee.FeatureCollection([result]),
  description:'DQ_to_'+TARGET_NAME+'_ModelTransfer_Result',
  folder:EXPORT_FOLDER,
  fileNamePrefix:'DQ_to_'+TARGET_NAME+'_ModelTransfer_Result',
  fileFormat:'CSV'
});

// ============================================================
// 14. EXPORT FIELD PREDICTIONS
// ============================================================

var targetPredExport = targetPred.select(['class','classification']);

Export.table.toDrive({
  collection:targetPredExport,
  description:'DQ_to_'+TARGET_NAME+'_ModelTransfer_FieldPrediction',
  folder:EXPORT_FOLDER,
  fileNamePrefix:'DQ_to_'+TARGET_NAME+'_ModelTransfer_FieldPrediction',
  fileFormat:'CSV'
});

// ============================================================
// 15. EXPORT FINAL MAP
// 0 background; 1 On; 2 Off; 3 Other
// ============================================================

var targetMapExport = targetMap
  .add(1)
  .unmask(0)
  .toByte()
  .rename('classification');

Export.image.toDrive({
  image:targetMapExport,
  description:'DQ_to_'+TARGET_NAME+'_ModelTransfer_2024',
  folder:EXPORT_FOLDER,
  fileNamePrefix:'DQ_to_'+TARGET_NAME+'_ModelTransfer_2024',
  region:geomTarget,
  scale:10,
  maxPixels:1e13
});