/**************** 0) ROI & 轻渲染设置 ****************/
var jjj = ee.FeatureCollection("projects/ee-l2892786691/assets/JJJ/jjj");
var region = jjj.geometry().transform('EPSG:4326', 1).buffer(0, 100);
Map.centerObject(region, 7);

var lczVis = {min:1, max:17, palette:[
  '8c0000','d10000','ff0000','bf4d00','ff6600','ff9955','faee05',
  'bcbcbc','ffccaa','555555','006a00','00aa00','648525','b9db79',
  '000000','fbf7ae','6a6aff'
]};

/**************** 1) 数据：GLCZ（高置信） & AEF(2018, 10m) ****************/
var GLCZ = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest')
  .filterBounds(region).mosaic().clip(region);
var lcz  = GLCZ.select('LCZ_Filter');        // 1..17
var prob = GLCZ.select('LCZ_Probability');   // 0..100
var TH = 90;
var lcz_highconf = lcz.updateMask(prob.gte(TH));  // 100m 标签（高置信）

// AEF 年度嵌入（10m）
var AEF2018 = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL')
  .filterDate('2018-01-01','2019-01-01')
  .filterBounds(jjj)   // 不 clip，避免边缘丢像元；后面推理再 clip(region)
  .mosaic();

var featBands = AEF2018.bandNames();
var aefProj   = AEF2018.projection();
print('AEF bands:', featBands);
print('AEF projection:', aefProj);

/**************** 2) 人工样本（Polygon）→ 轻度内缩打点（200/类） ****************/
var labelBand = 'Label';
var classes = ee.List.sequence(1, 17);
var N_MANUAL_PER_CLASS = 200;
var SHRINK_M = 20;
var SEED = 42;

var manualFC_raw = ee.FeatureCollection("projects/ee-l2892786691/assets/JJJ/2018jjjlow")
  .filterBounds(region)
  .filter(ee.Filter.gte('LCZ_label', 1))
  .filter(ee.Filter.lte('LCZ_label', 17));

// 轻度几何规范：向内收缩 20 m，避免边界混合像元
var manualPolys = manualFC_raw.map(function(f){
  var g3857 = f.geometry().transform('EPSG:3857', 1);
  var eroded = g3857.buffer(-SHRINK_M);
  var gFix   = ee.Geometry(ee.Algorithms.If(eroded.area(1).gt(0), eroded, g3857));
  var g4326  = gFix.transform('EPSG:4326', 1).intersection(region, 1);
  return ee.Feature(g4326).copyProperties(f, ['LCZ_label']);
}).map(function (f) { return f.set('area_m2', f.geometry().area(1)); })
  .filter(ee.Filter.gt('area_m2', 0));

var manualPts = ee.FeatureCollection(
  classes.map(function(c){
    c = ee.Number(c);
    var polys = manualPolys.filter(ee.Filter.eq('LCZ_label', c));
    var area  = polys.geometry().area(1);
    return ee.FeatureCollection(ee.Algorithms.If(
      area.gt(0),
      ee.FeatureCollection.randomPoints({
        region: polys.geometry(),
        points: N_MANUAL_PER_CLASS,
        seed: c.add(SEED).int(),
        maxError: 10
      }).map(function(p){ return ee.Feature(p.geometry()).set(labelBand, c); }),
      ee.FeatureCollection([])
    ));
  })
).flatten();
print('Manual points (<=200/class):', manualPts.size());

/**************** 3) GLCZ 高置信分层（100/类） ****************/
var proj100m = ee.Projection('EPSG:3857').atScale(100);
var labelHigh100 = lcz_highconf.rename(labelBand).toInt16().reproject(proj100m);

var N_HIGHCONF_PER_CLASS = 100;
var classValues = classes;
var classPoints = classes.map(function(_){ return ee.Number(N_HIGHCONF_PER_CLASS); });

var highconfPts = labelHigh100.stratifiedSample({
  numPoints: 0,
  classBand: labelBand,
  classValues: classValues,
  classPoints: classPoints,
  region: region,
  scale: 100,
  geometries: true,
  seed: SEED,
  tileScale: 2
});
print('High-conf points (<=100/class):', highconfPts.size());

/**************** 4) 合并点 → 在 AEF(10m) 上采样 ****************/
var allPts = manualPts.merge(highconfPts);
print('Total points (manual + highconf):', allPts.size());
Map.addLayer(allPts.style({pointSize: 2, color: 'cyan'}), {}, '样本点(预览)', false);

var samplesAll = AEF2018.sampleRegions({
  collection: allPts,
  properties: [labelBand],
  scale: 10,
  tileScale: 2,
  geometries: false
}).filter(ee.Filter.notNull([labelBand]));
print('Samples with AEF features:', samplesAll.size());

/**************** 5) 训练 / 留出验证（RF） ****************/
var TRAIN_FRAC = 0.7;
var withRand = samplesAll.randomColumn('rand', SEED);
var trainFC  = withRand.filter(ee.Filter.lt('rand', TRAIN_FRAC));
var testFC   = withRand.filter(ee.Filter.gte('rand', TRAIN_FRAC));
print('Train/Test sizes:', trainFC.size(), testFC.size());

var RF_TREES = 200;
var rf = ee.Classifier.smileRandomForest({numberOfTrees: RF_TREES, seed: SEED})
  .train({features: trainFC, classProperty: labelBand, inputProperties: featBands});

// 验证
var validated = testFC.classify(rf);
var cm = validated.errorMatrix(labelBand, 'classification');
print('Validation CM (AEF-only, 10m):', cm);
print('OA:', cm.accuracy(), 'Kappa:', cm.kappa());
print('UA:', cm.consumersAccuracy());
print('PA:', cm.producersAccuracy());

/**************** 6) 推理（10m） & 轻量显示 ****************/
var lcz2018_pred_10m = AEF2018.classify(rf)
  .rename('LCZ_2018_pred_AEF10m')
  .toInt16()
  .clip(region);

// 默认不开原始 10m 层（渲染重）；提供 120m 预览层以便浏览
var preview120 = lcz2018_pred_10m
  .toByte()
  .reproject({crs: aefProj.crs(), scale: 120});

Map.addLayer(preview120, lczVis, 'LCZ 2018 预测 (AEF-only, 120m 预览)', true);
Map.addLayer(lcz2018_pred_10m, lczVis, 'LCZ 2018 预测 (AEF-only, 10m 原图)', false);

/**************** 7) （可选）导出 ****************/
// Export.image.toAsset({
//   image: lcz2018_pred_10m.toByte(),
//   description: 'LCZ_2018_AEF10m_RF_manual200_highconf100',
//   assetId: 'users/你的用户名/LCZ_2018_AEF10m_RF_manual200_highconf100',
//   region: region,
//   scale: 10,
//   maxPixels: 1e13
// });
