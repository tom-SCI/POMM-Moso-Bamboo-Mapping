var roiTarget = ee.FeatureCollection("projects/ideal-8895/assets/longyou/LY"),
    fieldOnRaw = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lydnfield200"),
    fieldOffRaw = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lyxnfield_200"),
    otherTrainRaw = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lyqtfield200"),
    autoOnRaw = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lydn200gen"),
    autoOffRaw = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lyxn200gen"),
    targetOtherRaw = ee.FeatureCollection("projects/ideal-8895/assets/longyou/lyqttest200");
    
// -------------------- PARAMETERS --------------------
var TARGET_NAME = 'LY';

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
var EXPORT_FOLDER = 'POMM_WORKFLOW_2026';

var geomTarget = roiTarget.geometry();
Map.centerObject(roiTarget,10);

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

var trainOn = limitFC(centroidClass(autoOnRaw,0),N_ON,SAMPLE_SEED_ON,'train_on_random');
var trainOff = limitFC(centroidClass(autoOffRaw,1),N_OFF,SAMPLE_SEED_OFF,'train_off_random');
var trainOther = limitFC(centroidClass(otherTrainRaw,2),N_OTHER_TRAIN,OTHER_TRAIN_SEED,'train_other_random');
var trainSamples = trainOn.merge(trainOff).merge(trainOther);

var testOn = limitFC(centroidClass(fieldOnRaw,0),N_FIELD_PER_CLASS,FIELD_SEED_ON,'field_on_random');
var testOff = limitFC(centroidClass(fieldOffRaw,1),N_FIELD_PER_CLASS,FIELD_SEED_OFF,'field_off_random');
var testOther = limitFC(centroidClass(targetOtherRaw,2),N_FIELD_PER_CLASS,FIELD_SEED_OTHER,'field_other_random');
var independentTest = testOn.merge(testOff).merge(testOther);

print('Training On / Off / Other:',trainOn.size(),trainOff.size(),trainOther.size());
print('Independent test On / Off / Other:',testOn.size(),testOff.size(),testOther.size());

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
  var monthA = monthComposite(YEAR_A);
  var monthB = monthComposite(YEAR_B);

  var image = prefixBands(monthA.select(bands),'S2A')
    .addBands(prefixBands(monthB.select(bands),'S2B'))
    .toFloat();

  return {image:image,projection:p10};
}
// function nearCount(testFC,trainFC,distance){
//   var trainZone = trainFC.geometry().buffer(distance);
//   return testFC.filterBounds(trainZone).size();
// }

// print('================ SPATIAL LEAKAGE QC ================');

// print('ALL test within 10 m of ANY train:',
//   nearCount(independentTest,trainSamples,10));
// print('ALL test within 30 m of ANY train:',
//   nearCount(independentTest,trainSamples,30));
// print('ALL test within 100 m of ANY train:',
//   nearCount(independentTest,trainSamples,100));

// print('On test within 10 m of On train:',
//   nearCount(testOn,trainOn,10));
// print('On test within 30 m of On train:',
//   nearCount(testOn,trainOn,30));
// print('On test within 100 m of On train:',
//   nearCount(testOn,trainOn,100));

// print('Off test within 10 m of Off train:',
//   nearCount(testOff,trainOff,10));
// print('Off test within 30 m of Off train:',
//   nearCount(testOff,trainOff,30));
// print('Off test within 100 m of Off train:',
//   nearCount(testOff,trainOff,100));

// print('Other test within 10 m of Other train:',
//   nearCount(testOther,trainOther,10));
// print('Other test within 30 m of Other train:',
//   nearCount(testOther,trainOther,30));
// print('Other test within 100 m of Other train:',
//   nearCount(testOther,trainOther,100));
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
// 4. FEATURE IMAGES
// ============================================================

var s2Data = buildS2Data(geomTarget);
var aefImage = buildAEF(geomTarget);

var featureImage = aefImage.addBands(s2Data.image).toFloat();

var featureBands = featureImage.bandNames();
var aefBands = aefImage.bandNames();

print('Full feature number:',featureBands.length());
print('AEF feature number:',aefBands.length());

// ============================================================
// 5. EXTRACT TRAINING FEATURES
// ============================================================

var trainingFull = featureImage.select(featureBands).sampleRegions({
  collection:trainSamples,
  properties:['class'],
  scale:10,
  projection:s2Data.projection,
  geometries:false,
  tileScale:16
}).filter(ee.Filter.notNull(featureBands));

var trainingAEF = aefImage.select(aefBands).sampleRegions({
  collection:trainSamples,
  properties:['class'],
  scale:10,
  projection:s2Data.projection,
  geometries:false,
  tileScale:16
}).filter(ee.Filter.notNull(aefBands));

print('Effective training n:',trainingFull.size());
print('Training histogram:',trainingFull.aggregate_histogram('class'));

// ============================================================
// 6. TRAIN INTEGRATED CLASSIFICATION SYSTEM
// ============================================================

var localRFFull = ee.Classifier.smileRandomForest({
  numberOfTrees:N_TREES,
  variablesPerSplit:null,
  minLeafPopulation:1,
  bagFraction:BAG_FRACTION,
  maxNodes:null,
  seed:RF_SEED
}).train({
  features:trainingFull,
  classProperty:'class',
  inputProperties:featureBands
});

var localRFAEF = ee.Classifier.smileRandomForest({
  numberOfTrees:N_TREES,
  variablesPerSplit:null,
  minLeafPopulation:1,
  bagFraction:BAG_FRACTION,
  maxNodes:null,
  seed:RF_SEED
}).train({
  features:trainingAEF,
  classProperty:'class',
  inputProperties:aefBands
});

print(TARGET_NAME+' integrated RF training completed');

// ============================================================
// 7. INTEGRATED CLASSIFICATION
// ============================================================

var fullClassification = featureImage.select(featureBands)
  .classify(localRFFull)
  .rename('classification');

var aefClassification = aefImage.select(aefBands)
  .classify(localRFAEF)
  .rename('classification');

var rawClassification = fullClassification
  .unmask(aefClassification)
  .rename('classification')
  .clip(geomTarget);

// ============================================================
// 8. INDEPENDENT TEST
// ============================================================
var localPrediction = fullClassification.sampleRegions({
  collection:independentTest,
  properties:['class'],
  scale:10,
  projection:s2Data.projection,
  geometries:true,
  tileScale:16
}).filter(ee.Filter.notNull(['classification']));
// var localPrediction = rawClassification.sampleRegions({
//   collection:independentTest,
//   properties:['class'],
//   scale:10,
//   projection:s2Data.projection,
//   geometries:true,
//   tileScale:16
// }).filter(ee.Filter.notNull(['classification']));

print('================================================');
print('WORKFLOW TRANSFER — INDEPENDENT TEST');
print('================================================');
print('Original test n:',independentTest.size());
print('Original class histogram:',independentTest.aggregate_histogram('class'));
print('Effective test n:',localPrediction.size());
print('Effective class histogram:',localPrediction.aggregate_histogram('class'));
print('Dropped samples:',independentTest.size().subtract(localPrediction.size()));

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

var cm = localPrediction.errorMatrix('class','classification',[0,1,2]);
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
  target_region:TARGET_NAME,
  transfer_mode:'Workflow_transfer',
  feature_mode:'AEF_S2',

  year_a:YEAR_A,
  year_b:YEAR_B,
  ref_month:REF_MONTH,

  rf_trees:N_TREES,
  rf_bag_fraction:BAG_FRACTION,
  rf_seed:RF_SEED,

  train_n:trainingFull.size(),
  test_n:localPrediction.size(),

  train_On_n:trainingFull.filter(ee.Filter.eq('class',0)).size(),
  train_Off_n:trainingFull.filter(ee.Filter.eq('class',1)).size(),
  train_Other_n:trainingFull.filter(ee.Filter.eq('class',2)).size(),

  test_On_n:localPrediction.filter(ee.Filter.eq('class',0)).size(),
  test_Off_n:localPrediction.filter(ee.Filter.eq('class',1)).size(),
  test_Other_n:localPrediction.filter(ee.Filter.eq('class',2)).size(),

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

print('FINAL WORKFLOW TRANSFER RESULT:',result);

// ============================================================
// 11. BOOTSTRAP TABLE
// ============================================================

var bootstrapTable = localPrediction.select(
  ['class','classification'],
  ['truth','pred']
);

print('Bootstrap first:',bootstrapTable.first());
print('Bootstrap n:',bootstrapTable.size());

Export.table.toDrive({
  collection:bootstrapTable,
  description:TARGET_NAME+'_Workflow_Bootstrap',
  folder:'POMM_WORKFLOW_2026_num',
  fileNamePrefix:TARGET_NAME+'_Workflow_Bootstrap'+N_ON,
  fileFormat:'CSV',
  selectors:['truth','pred']
});
// ============================================================
// 12. BROAD FOREST DOMAIN MASK
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

var finalMap = rawClassification
  .updateMask(forestDomainMask)
  .clip(geomTarget);

// ============================================================
// 13. FINAL MAP VALIDATION
// ============================================================

var finalValidation = finalMap.unmask(3).sampleRegions({
  collection:independentTest,
  properties:['class'],
  scale:10,
  projection:s2Data.projection,
  geometries:true,
  tileScale:16
});

var maskedSamples = finalValidation.filter(ee.Filter.eq('classification',3));
var retainedSamples = finalValidation.filter(ee.Filter.neq('classification',3));

var cmRetained = retainedSamples.errorMatrix('class','classification',[0,1,2]);

var totalN = finalValidation.size();
var correctN = finalValidation.filter(ee.Filter.equals({
  leftField:'class',
  rightField:'classification'
})).size();

var endToEndOA = ee.Number(correctN).divide(totalN);

print('================================================');
print('FINAL WORKFLOW MAP VALIDATION');
print('================================================');
print('Total field samples:',totalN);
print('Retained samples:',retainedSamples.size());
print('Masked-out samples:',maskedSamples.size());
print('Masked-out histogram:',maskedSamples.aggregate_histogram('class'));
print('Retained CM:',cmRetained);
print('Retained OA:',cmRetained.accuracy());
print('Retained Kappa:',cmRetained.kappa());
print('Retained PA:',cmRetained.producersAccuracy());
print('Retained UA:',cmRetained.consumersAccuracy());
print('Retained F1:',cmRetained.fscore(1));
print('End-to-end OA:',endToEndOA);

// // ============================================================
// // 14. MAP DISPLAY
// // ============================================================

Map.addLayer(forestDomainMask,{
  palette:['006400']
},'Broad Forest Domain',false);

Map.addLayer(finalMap,{
  min:0,
  max:2,
  palette:['006400','90EE90','D2B48C']
},TARGET_NAME+' Workflow Three-Class Map',true);

// // ============================================================
// // 15. EXPORT RESULT
// // ============================================================

// Export.table.toDrive({
//   collection:ee.FeatureCollection([result]),
//   description:TARGET_NAME+'_WorkflowTransfer_Result',
//   folder:'POMM_WORKFLOW_2026_num',
//   fileNamePrefix:TARGET_NAME+'_WorkflowTransfer_Result'+N_ON,
//   fileFormat:'CSV'
// });

// ============================================================
// 16. EXPORT FIELD PREDICTION
// ============================================================

var predictionExport = localPrediction.select(['class','classification']);

// Export.table.toDrive({
//   collection:predictionExport,
//   description:TARGET_NAME+'_WorkflowTransfer_FieldPrediction',
//   folder:'POMM_WORKFLOW_2026',
//   fileNamePrefix:TARGET_NAME+'_WorkflowTransfer_FieldPrediction',
//   fileFormat:'CSV'
// });

// ============================================================
// 17. EXPORT FINAL MAP
// 0 background; 1 On; 2 Off; 3 Other
// ============================================================

var finalMapExport = finalMap.add(1).unmask(0).toByte().rename('classification');

Export.image.toDrive({
  image:finalMapExport,
  description:TARGET_NAME+'_Workflow_ThreeClass_2024',
  folder:'GEE_Workflow_Map',
  fileNamePrefix:TARGET_NAME+'_Workflow_ThreeClass_2024',
  region:geomTarget,
  scale:10,
  maxPixels:1e13
});

// ============================================================
// 18. EXPORT FINAL VALIDATION
// ============================================================

// Export.table.toDrive({
//   collection:finalValidation,
//   description:TARGET_NAME+'_Workflow_FinalValidation',
//   folder:'POMM_WORKFLOW_2026',
//   fileNamePrefix:TARGET_NAME+'_Workflow_FinalValidation',
//   fileFormat:'CSV'
// });