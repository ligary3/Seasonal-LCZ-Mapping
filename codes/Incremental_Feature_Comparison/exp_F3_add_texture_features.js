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

// ===== 指数 =====
var ndvi  = L8_2018_raw.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');
var ndbi  = L8_2018_raw.normalizedDifference(['SR_B6','SR_B5']).rename('NDBI');
var mndwi = L8_2018_raw.normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI');

// ===== 纹理（GLCM）在 100 m 上计算 → 上到 30 m（靠 reproject） =====
// ===== 纹理（GLCM）在 30 m 上计算，窗口 5×5，64 灰阶 =====
var texSrc = L8_2018_raw.select(['SR_B5','SR_B6'])
  // 将反射率线性映射到 [0, 63]，再量化为 6-bit 整型
  .unitScale(0, 0.4).multiply(63).int16()   // 0.4 是经验上限，可按需改到 0.5
  .updateMask(L8_2018_raw.select(0).mask())
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

// 5×5 窗口（≈150 m 语境，通常对 LCZ 更友好）
var glcm30 = texSrc.glcmTexture({size: 5});

// 只取稳健的统计量（对噪声更稳）
var glcmSel30m = glcm30.select([
  'SR_B5_contrast','SR_B5_idm',
  'SR_B6_contrast','SR_B6_idm'
]).reproject(proj3857_30).setDefaultProjection(proj3857_30);

// 重新拼装特征栈（指数 + 改良纹理）
var L8_feat_30m = L8_2018_raw.addBands([ndvi, ndbi, mndwi])
  .addBands(glcmSel30m)
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

var l8Bands = L8_feat_30m.bandNames();
print('L8+指数+纹理(30m-GLCM-5x5-64gray) bands:', l8Bands);


// ====================== 3) 人工200/类（带轻度内缩） + 高置信100/类 ======================
var labelBand = 'Label';
var classes = ee.List.sequence(1, 17);
var N_MANUAL_PER_CLASS = 200;
var N_HIGHCONF_PER_CLASS = 100;
var SEED = 42;
var SHRINK_M = 20;   // 轻度向内收缩（米）

// 3a) 人工标注面（轻度内缩） → 随机点
var manualFC_raw = ee.FeatureCollection("projects/ee-l2892786691/assets/JJJ/2018jjjlow")
  .filterBounds(region)
  .filter(ee.Filter.gte('LCZ_label', 1))
  .filter(ee.Filter.lte('LCZ_label', 17));

var manualPts = ee.FeatureCollection(
  classes.map(function(c){
    c = ee.Number(c);
    var polysShrink = manualFC_raw
      .filter(ee.Filter.eq('LCZ_label', c))
      .map(function(f){
        var g = f.geometry().buffer(-SHRINK_M);                 // 内缩
        var gFixed = ee.Geometry(ee.Algorithms.If(g.area(1).gt(0), g, f.geometry()));
        return ee.Feature(gFixed);
      });
    var geom = ee.FeatureCollection(polysShrink).geometry();
    var area = geom.area(1);

    return ee.FeatureCollection(ee.Algorithms.If(
      area.gt(0),
      ee.FeatureCollection.randomPoints({
        region: geom,
        points: N_MANUAL_PER_CLASS,
        seed: c.add(SEED).int(),
        maxError: 10
      }).map(function(p){ return ee.Feature(p.geometry()).set(labelBand, c); }),
      ee.FeatureCollection([])
    ));
  })
).flatten();
print('Manual points (≤200/class, shrink='+SHRINK_M+'m):', manualPts.size());

// 3b) 高置信 GLCZ（100 m 分层）
var labelImg100 = lcz_highconf.rename(labelBand).reproject(proj3857_100);
var classValues = classes;
var classPoints = classes.map(function(_) { return ee.Number(N_HIGHCONF_PER_CLASS); });

var highconfPts = labelImg100.stratifiedSample({
  numPoints: 0,                 // 完全由 classPoints 控制
  classBand: labelBand,
  classValues: classValues,
  classPoints: classPoints,
  region: region,
  scale: 100,
  geometries: true,
  seed: SEED,
  tileScale: 2
});
print('High-conf points (≤100/class):', highconfPts.size());

// 合并两批点（只决定点位与标签）
var allPts = manualPts.merge(highconfPts);
print('Total labeled points (manual+highconf):', allPts.size());
Map.addLayer(allPts.style({pointSize: 2, color: 'cyan'}), {}, '样本点(预览)', false);

// ====================== 4) 在 30 m 特征上提取样本 ======================
var samplesAll = L8_feat_30m.sampleRegions({
  collection: allPts,
  properties: [labelBand],
  scale: 30,           // 与 L8_feat_30m 一致（EPSG:3857, 30m）
  tileScale: 2,
  geometries: false
}).filter(ee.Filter.notNull([labelBand]));
print('Samples with L8+指数+纹理(30m) features:', samplesAll.size());

// ====================== 5) 划分训练/验证，训练 RF ======================
var TRAIN_FRAC = 0.7;
var withRand = samplesAll.randomColumn('rand', SEED);
var trainFC  = withRand.filter(ee.Filter.lt('rand', TRAIN_FRAC));
var testFC   = withRand.filter(ee.Filter.gte('rand', TRAIN_FRAC));
print('Train size:', trainFC.size(), 'Test size:', testFC.size());

var RF_TREES = 200;
var rf = ee.Classifier.smileRandomForest({numberOfTrees: RF_TREES, seed: SEED})
  .train({features: trainFC, classProperty: labelBand, inputProperties: l8Bands});

// ====================== 6) 留出法评估 ======================
var validated = testFC.classify(rf);
var cm = validated.errorMatrix(labelBand, 'classification');
print('Validation confusion matrix (L8+指数+纹理, 30m)', cm);
print('Validation OA (L8+指数+纹理, 30m)', cm.accuracy());
print('Validation Kappa (L8+指数+纹理, 30m)', cm.kappa());
print('UA per class (L8+指数+纹理, 30m)', cm.consumersAccuracy());
print('PA per class (L8+指数+纹理, 30m)', cm.producersAccuracy());

// ====================== 7) 30 m 推理输出 ======================
var lcz2018_pred_30m = L8_feat_30m.classify(rf)
  .rename('LCZ_2018_pred_L8IdxTex')     // 标识为 L8+指数+纹理
  .toInt16()
  .clip(region);
Map.addLayer(lcz2018_pred_30m, lczVis, 'LCZ 2018 预测 (L8+指数+纹理, 30 m)', true);

// ====================== 8) （可选）导出 ======================
// Export.image.toAsset({
//   image: lcz2018_pred_30m,
//   description: 'LCZ_2018_L8IdxTex_RF_30m',
//   assetId: 'users/你的用户名/LCZ_2018_L8IdxTex_RF_30m',
//   region: region,
//   scale: 30,
//   maxPixels: 1e13
// });
