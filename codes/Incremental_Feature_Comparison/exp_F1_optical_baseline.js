// ====================== 0) ROI：统一到 WGS84 + 安全修复 ======================
var jjj = ee.FeatureCollection("projects/ee-l2892786691/assets/JJJ/jjj");
var ROI4326 = jjj.geometry().transform('EPSG:4326', 1);
var region  = ROI4326.buffer(0, 100);   // 0 距离 + 100 m 容差，修复几何
Map.centerObject(region, 7);

// 统一米制投影（关键修复）
var proj3857_30  = ee.Projection('EPSG:3857').atScale(30);
var proj3857_100 = ee.Projection('EPSG:3857').atScale(100);

// ====================== 1) 读取 GLCZ 并裁剪 ======================
var GLCZ = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest')
  .filterBounds(region).mosaic().clip(region);

var lcz  = GLCZ.select('LCZ_Filter');        // 1..17
var prob = GLCZ.select('LCZ_Probability');   // 0..100（百分制）

var TH = 90;                                  // 置信度阈值（可试 80/90/95）
var lcz_highconf = lcz.updateMask(prob.gte(TH));

// 可视化
var lczVis = {min:1, max:17, palette:[
  '8c0000','d10000','ff0000','bf4d00','ff6600','ff9955','faee05',
  'bcbcbc','ffccaa','555555','006a00','00aa00','648525','b9db79',
  '000000','fbf7ae','6a6aff'
]};
Map.addLayer(lcz_highconf, lczVis, 'GLCZ 高置信(>=TH)', true);

// 统计（确认各类在高置信范围内有像元）
var lczHist = lcz_highconf.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(),
  geometry: region, scale: 100, maxPixels: 1e13, bestEffort: true
});
var histDict = ee.Dictionary(
  ee.Algorithms.If(lczHist.contains('LCZ_Filter'), lczHist.get('LCZ_Filter'), ee.Dictionary({}))
);
print('LCZ counts (>=TH, 100m):', histDict);

// ====================== 2) Landsat-8 2018 年度合成（30 m） ======================
function prepSrL8(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);

  // 官方系数缩放（反射率 & 亮温）
  var getFactorImg = function (factorNames) {
    var factorList = image.toDictionary().select(factorNames).values();
    return ee.Image.constant(factorList);
  };
  var scaleImg  = getFactorImg(['REFLECTANCE_MULT_BAND_.|TEMPERATURE_MULT_BAND_ST_B10']);
  var offsetImg = getFactorImg(['REFLECTANCE_ADD_BAND_.|TEMPERATURE_ADD_BAND_ST_B10']);
  var scaled = image.select('SR_B.|ST_B10').multiply(scaleImg).add(offsetImg);

  return image.addBands(scaled, null, true)
              .updateMask(qaMask).updateMask(saturationMask);
}

var L8_2018_raw = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(region)
  .filterDate('2018-01-01', '2019-01-01')
  .map(prepSrL8)
  .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10'])
  .median();

// // 指数（如需“纯 SR 对照”，可注释掉下面三行）
// var ndvi  = L8_2018_raw.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');
// var ndbi  = L8_2018_raw.normalizedDifference(['SR_B6','SR_B5']).rename('NDBI');
// var mndwi = L8_2018_raw.normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI');

// —— 关键：统一把 30 m 特征影像重投影到 EPSG:3857，并设为默认投影 ——
// （避免选到错误的 UTM 分带）
var L8_feat_30m = L8_2018_raw//.addBands([ndvi, ndbi, mndwi])
  .reproject(proj3857_30)
  .setDefaultProjection(proj3857_30);

var l8Bands = L8_feat_30m.bandNames();
print('L8 feature bands (30m, EPSG:3857):', l8Bands);

// ====================== 3) 用“标签(100 m)”分层抽样生成点（只决定点位） ======================
// ====================== 3) 手工样本(Polygon) → 每类200随机点 ======================
// 手工标注面：「projects/ee-l2892786691/assets/JJJ/2018jjjlow」字段名 LCZ_label (1..17)
var labelBand   = 'Label';
var manualFC_raw = ee.FeatureCollection("projects/ee-l2892786691/assets/JJJ/2018jjjlow")
  .filterBounds(region)
  .filter(ee.Filter.gte('LCZ_label', 1))
  .filter(ee.Filter.lte('LCZ_label', 17));

// 轻度几何清洗：向内收缩 & 简化，减少边界混合像元（可按需调整或关掉）
var MARGIN_M = 20, SIMP_TOL = 10;
var manualPolys = manualFC_raw.map(function(f){
  var g3857 = f.geometry().transform('EPSG:3857', 1);
  var eroded = g3857.buffer(-MARGIN_M);
  var gClean = ee.Geometry(ee.Algorithms.If(eroded.area(1).gt(0), eroded, g3857)).simplify(SIMP_TOL);
  var g4326 = gClean.transform('EPSG:4326', 1).intersection(region, 1);
  return ee.Feature(g4326).copyProperties(f, ['LCZ_label']);
}).map(function (f) { return f.set('area_m2', f.geometry().area(1)); })
  .filter(ee.Filter.gt('area_m2', 0));

// 按类随机打点（每类最多200点；类里面积为0则跳过）
var classes = ee.List.sequence(1, 17);
var N_MANUAL_PER_CLASS = 200;
var SEED = 42;

var manualPts = ee.FeatureCollection(
  classes.map(function(c){
    c = ee.Number(c);
    var polys = manualPolys.filter(ee.Filter.eq('LCZ_label', c));
    var area = polys.geometry().area(1);
    var pts = ee.FeatureCollection(ee.Algorithms.If(
      area.gt(0),
      ee.FeatureCollection.randomPoints({
        region: polys.geometry(),
        points: N_MANUAL_PER_CLASS,
        seed: c.add(SEED).int(),
        maxError: 10
      }).map(function(p){ return ee.Feature(p.geometry()).set('Label', c); }),
      ee.FeatureCollection([])
    ));
    return pts;
  })
).flatten();
print('Manual points by polygons (200/class max):', manualPts.size());

// ====================== 3b) 高置信 GLCZ(100 m) → 每类100分层点 ======================
var proj100m = ee.Projection('EPSG:3857').atScale(100);
var labelHigh100 = lcz_highconf.rename('Label').toInt16().reproject(proj100m);

// 每类上限100；若某类高置信像元少，则 stratifiedSample 会自动少取
var N_HIGHCONF_PER_CLASS = 100;
var classValues = classes; // 1..17
var classPoints = classes.map(function(_) { return ee.Number(N_HIGHCONF_PER_CLASS); });

var highconfPts = labelHigh100.stratifiedSample({
  numPoints: 0,                 // 完全由 classPoints 控制
  classBand: 'Label',
  classValues: classValues,
  classPoints: classPoints,
  region: region,
  scale: 100,
  geometries: true,
  seed: SEED,
  tileScale: 2
});
print('High-conf points from GLCZ (100/class max):', highconfPts.size());

// ====================== 4) 合并两批点 → 在 30 m L8 特征上取值 ======================
var allPts = manualPts.merge(highconfPts);
print('Total points (manual 200/class + highConf 100/class):', allPts.size());

var samplesAll = L8_feat_30m.sampleRegions({
  collection: allPts,
  properties: ['Label'],
  scale: 30,          // 与 L8_feat_30m 一致（EPSG:3857, 30m）
  tileScale: 2,
  geometries: false
}).filter(ee.Filter.notNull(['Label']));
print('Samples with L8(30m) features:', samplesAll.size());

// =========（后续与原来一致）训练/验证 =========
var TRAIN_FRAC = 0.7;
var withRand = samplesAll.randomColumn('rand', SEED);
var trainFC  = withRand.filter(ee.Filter.lt('rand', TRAIN_FRAC));
var testFC   = withRand.filter(ee.Filter.gte('rand', TRAIN_FRAC));
print('Train size:', trainFC.size(), 'Test size:', testFC.size());

var RF_TREES = 200;
var rf = ee.Classifier.smileRandomForest({numberOfTrees: RF_TREES, seed: SEED})
  .train({features: trainFC, classProperty: 'Label', inputProperties: l8Bands});

// ========= 留出评估 & 推理（保留你原来的第6、7节） =========


// ====================== 6) 留出法评估 ======================
var validated = testFC.classify(rf);
var cm = validated.errorMatrix(labelBand, 'classification');
print('Validation confusion matrix (L8, 30m)', cm);
print('Validation OA (L8, 30m)', cm.accuracy());
print('Validation Kappa (L8, 30m)', cm.kappa());
print('UA per class (L8, 30m)', cm.consumersAccuracy());
print('PA per class (L8, 30m)', cm.producersAccuracy());

// ====================== 7) 30 m 推理输出 ======================
var lcz2018_pred_30m = L8_feat_30m.classify(rf)
  .rename('LCZ_2018_pred_L8')
  .toInt16()
  .clip(region); // 可用 region 或 region3857，EE 会处理投影变换
Map.addLayer(lcz2018_pred_30m, lczVis, 'LCZ 2018 预测 (L8, 30 m)', true);


// 如需把结果导出成表（可选）：
// Export.table.toDrive({
//   collection: ee.FeatureCollection([ee.Feature(null, stats_L8_30m)]),
//   description: 'LCZ_MissingStats_L8_30m',
//   fileFormat: 'CSV'
// });

// ====================== 8) （可选）导出 ======================
// Export.image.toAsset({
//   image: lcz2018_pred_30m,
//   description: 'LCZ_2018_L8_RF_30m',
//   assetId: 'users/你的用户名/LCZ_2018_L8_RF_30m',
//   region: region,
//   scale: 30,
//   maxPixels: 1e13
// });
