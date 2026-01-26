/******************** 0) ROI & 投影 *************************/
var jjj = ee.FeatureCollection("projects/ee-l2892786691/assets/JJJ/jjj");
var ROI4326 = jjj.geometry().transform('EPSG:4326', 1);
var region  = ROI4326.buffer(0, 100);
Map.centerObject(region, 7);

var proj3857_30  = ee.Projection('EPSG:3857').atScale(30);
var proj3857_100 = ee.Projection('EPSG:3857').atScale(100);

/******************** 1) GLCZ 高置信标签 (100 m) *************/
var GLCZ = ee.ImageCollection('RUB/RUBCLIM/LCZ/global_lcz_map/latest')
  .filterBounds(region).mosaic().clip(region);

var lcz  = GLCZ.select('LCZ_Filter');        // 1..17
var prob = GLCZ.select('LCZ_Probability');   // 0..100
var TH   = 90;
var lcz_highconf = lcz.updateMask(prob.gte(TH));

var lczVis = {min:1, max:17, palette:[
  '8c0000','d10000','ff0000','bf4d00','ff6600','ff9955','faee05',
  'bcbcbc','ffccaa','555555','006a00','00aa00','648525','b9db79',
  '000000','fbf7ae','6a6aff'
]};
Map.addLayer(lcz_highconf, lczVis, 'GLCZ 高置信(>=TH)', true);

// 统计直方图（确认各类有样本）
var lczHist = lcz_highconf.reduceRegion({
  reducer: ee.Reducer.frequencyHistogram(),
  geometry: region, scale: 100, maxPixels: 1e13, bestEffort: true
});
var histDict = ee.Dictionary(
  ee.Algorithms.If(lczHist.contains('LCZ_Filter'), lczHist.get('LCZ_Filter'), ee.Dictionary({}))
);
print('LCZ counts (>=TH, 100m):', histDict);

/******************** 2) L8 & S2 预处理 (2018) ****************/
var START = '2018-01-01', END = '2019-01-01';

// ---- L8 预处理：缩放到物理量 + 云/饱和掩膜
function prepSrL8(img) {
  var qaMask = img.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var satMask = img.select('QA_RADSAT').eq(0);
  var getFactorImg = function (names) {
    return ee.Image.constant(img.toDictionary().select(names).values());
  };
  var scale  = getFactorImg(['REFLECTANCE_MULT_BAND_.|TEMPERATURE_MULT_BAND_ST_B10']);
  var offset = getFactorImg(['REFLECTANCE_ADD_BAND_.|TEMPERATURE_ADD_BAND_ST_B10']);
  var scaled = img.select('SR_B.|ST_B10').multiply(scale).add(offset)
                  .updateMask(qaMask).updateMask(satMask);
  var sr = scaled.select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7']).clamp(0,1).toFloat();
  var st = scaled.select('ST_B10').toFloat();
  return sr.addBands(st);
}

// ---- S2 高质量云影掩膜：s2cloudless + SCL
function prepS2_cloudprob(img) {
  var prob = ee.Image(img.get('cloud_probability')); // from join
  var scl  = img.select('SCL');
  var cloud  = prob.gt(40);  // 阈值 35~60 可调
  var shadow = scl.eq(3);
  var cirrus = scl.eq(10);
  var snow   = scl.eq(11);
  var mask = cloud.or(shadow).or(cirrus).or(snow).not();

  var s2Bands = ['B2','B3','B4','B8','B11','B12'];
  var l8Names = ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7'];

  var sr = img.select(s2Bands).divide(10000).rename(l8Names)
              .updateMask(mask).clamp(0,1).toFloat();

  var stB10_placeholder = ee.Image(0).updateMask(ee.Image(0)).rename('ST_B10').toFloat();
  return sr.addBands(stB10_placeholder);
}

// L8 合成
var L8_SR = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(region).filterDate(START, END)
  .map(prepSrL8)
  .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10'])
  .median().clip(region).toFloat();

// S2 合成（join cloud prob）
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(region).filterDate(START, END)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80));

var s2prob = ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
  .filterBounds(region).filterDate(START, END);

var joined = ee.ImageCollection(ee.Join.saveFirst('cloud_probability').apply({
  primary: s2,
  secondary: s2prob,
  condition: ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'})
}));

var S2_SR = joined.map(prepS2_cloudprob)
  .select(['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','ST_B10'])
  .median().clip(region).toFloat();

/******************** 3) 光学特征：分传感器并排 ****************/
// L8 指数
var L8_NDVI = L8_SR.normalizedDifference(['SR_B5','SR_B4']).rename('L8_NDVI');
var L8_NDBI = L8_SR.normalizedDifference(['SR_B6','SR_B5']).rename('L8_NDBI');
var L8_MNDWI= L8_SR.normalizedDifference(['SR_B3','SR_B6']).rename('L8_MNDWI');

// L8 纹理 (30m, 5x5, 64灰阶，仅 B5/B6)
var texSrcL8 = L8_SR.select(['SR_B5','SR_B6'])
  .unitScale(0, 0.4).multiply(63).int16()
  .updateMask(L8_SR.select(0).mask())
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

var glcmL8 = texSrcL8.glcmTexture({size: 5}).select([
  'SR_B5_contrast','SR_B5_idm','SR_B6_contrast','SR_B6_idm'
]).rename(['L8_B5_contrast','L8_B5_idm','L8_B6_contrast','L8_B6_idm'])
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

// S2 指数（不做纹理）
var S2_NDVI = S2_SR.normalizedDifference(['SR_B5','SR_B4']).rename('S2_NDVI');
var S2_NDBI = S2_SR.normalizedDifference(['SR_B6','SR_B5']).rename('S2_NDBI');
var S2_MNDWI= S2_SR.normalizedDifference(['SR_B3','SR_B6']).rename('S2_MNDWI');

// 光学特征并排（保持 30m）
var Opt_feat_30m = ee.Image.cat(
  L8_SR.rename(['L8_B2','L8_B3','L8_B4','L8_B5','L8_B6','L8_B7','L8_ST'])
       .addBands([L8_NDVI, L8_NDBI, L8_MNDWI, glcmL8]),
  S2_SR.rename(['S2_B2','S2_B3','S2_B4','S2_B5','S2_B6','S2_B7','S2_ST_placeholder'])
       .select(['S2_B2','S2_B3','S2_B4','S2_B5','S2_B6','S2_B7'])
       .addBands([S2_NDVI, S2_NDBI, S2_MNDWI])
).reproject(proj3857_30).setDefaultProjection(proj3857_30);

var optMask = Opt_feat_30m.select(0).mask();

/******************** 4) Sentinel-1 特征 (30 m) ****************/
var eps = ee.Image.constant(1e-6);
var S1_2018 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(region)
  .filterDate(START, END)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) // 如需双轨：改 ee.Filter.or(...)
  .select(['VV','VH']);

var S1_lin = S1_2018.median().clip(region);
var VV_lin = S1_lin.select('VV').max(eps);
var VH_lin = S1_lin.select('VH').max(eps);

var VV_dB = VV_lin.log10().multiply(10).rename('VV_dB');
var VH_dB = VH_lin.log10().multiply(10).rename('VH_dB');
var VV_VH_ratio  = VV_lin.divide(VH_lin).rename('VV_VH_ratio');
var VV_VH_dBdiff = VV_dB.subtract(VH_dB).rename('VV_VH_dBdiff');

var s1_fill = ee.Image.cat([
  ee.Image.constant(-20).rename('VV_dB'),
  ee.Image.constant(-25).rename('VH_dB'),
  ee.Image.constant(2.0).rename('VV_VH_ratio'),
  ee.Image.constant(5.0).rename('VV_VH_dBdiff')
]);

var S1_feats_30m = VV_dB.addBands([VH_dB, VV_VH_ratio, VV_VH_dBdiff])
  .unmask(s1_fill)
  .reproject(proj3857_30).setDefaultProjection(proj3857_30)
  .updateMask(optMask);
print('S1 feature bands (30m):', S1_feats_30m.bandNames());

/******************** 5) 融合特征（全部 30 m） *****************/
var Feat_30m = Opt_feat_30m.addBands(S1_feats_30m)
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);
var inputBands = Feat_30m.bandNames();
print('All feature bands (L8 + S2 + 指数 + L8纹理 + S1, 30m):', inputBands);

/******************** 6) 采样点：人工200/类(20m内缩) + 高置信100/类 ****/
var labelBand = 'Label';
var classes = ee.List.sequence(1, 17);
var N_MANUAL_PER_CLASS = 200;
var N_HIGHCONF_PER_CLASS = 100;
var SEED = 42;
var SHRINK_M = 20;

// 人工标注面（轻度内缩） → 随机点
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
        var g = f.geometry().buffer(-SHRINK_M);
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

// 高置信分层（100 m）
var labelImg100 = lcz_highconf.rename(labelBand).reproject(proj3857_100);
var classValues = classes;
var classPoints = classes.map(function(_) { return ee.Number(N_HIGHCONF_PER_CLASS); });

var highconfPts = labelImg100.stratifiedSample({
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
print('High-conf points (≤100/class):', highconfPts.size());

// 合并两批点
var allPts = manualPts.merge(highconfPts);
print('Total labeled points (manual+highconf):', allPts.size());
Map.addLayer(allPts.style({pointSize: 2, color: 'cyan'}), {}, '样本点(预览)', false);

/******************** 7) 样本提取 (30 m)、训练 RF、评估 *********/
var samplesAll = Feat_30m.sampleRegions({
  collection: allPts,
  properties: [labelBand],
  scale: 30,
  tileScale: 2,
  geometries: false
}).filter(ee.Filter.notNull([labelBand]));
print('Samples with features (30m):', samplesAll.size());

var TRAIN_FRAC = 0.7;
var withRand = samplesAll.randomColumn('rand', SEED);
var trainFC  = withRand.filter(ee.Filter.lt('rand', TRAIN_FRAC));
var testFC   = withRand.filter(ee.Filter.gte('rand', TRAIN_FRAC));
print('Train size:', trainFC.size(), 'Test size:', testFC.size());

var RF_TREES = 200;
// 可按需：variablesPerSplit: Math.floor(inputBands.length/3)
var rf = ee.Classifier.smileRandomForest({numberOfTrees: RF_TREES, seed: SEED})
  .train({features: trainFC, classProperty: labelBand, inputProperties: inputBands});

var validated = testFC.classify(rf);
var cm = validated.errorMatrix(labelBand, 'classification');
print('Validation confusion matrix (L8+S2+索引+L8纹理+S1, 30m)', cm);
print('Validation OA', cm.accuracy());
print('Validation Kappa', cm.kappa());
print('UA per class', cm.consumersAccuracy());
print('PA per class', cm.producersAccuracy());

// （可选）看特征重要性
var expl = rf.explain();
print('RF importance:', ee.Dictionary(expl).get('importance'));

/******************** 8) 推理输出 (30 m) ************************/
var lcz2018_pred_30m = Feat_30m.classify(rf)
  .rename('LCZ_2018_pred_L8S2_L8Tex_S1')
  .toInt16()
  .clip(region);

Map.addLayer(lcz2018_pred_30m, lczVis, 'LCZ 2018 预测 (L8+S2+指数+L8纹理+S1, 30 m)', true);

/******************** 9) （可选）导出 ***************************/
// Export.image.toAsset({
//   image: lcz2018_pred_30m,
//   description: 'LCZ_2018_L8S2_Idx_L8Tex_S1_RF_30m',
//   assetId: 'users/你的用户名/LCZ_2018_L8S2_Idx_L8Tex_S1_RF_30m',
//   region: region,
//   scale: 30,
//   maxPixels: 1e13
// });
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
var TH = 90;                                 // 置信度阈值（可试 80/90/95）
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
  .median()
  .clip(region);

// ===== 指数 =====
var ndvi  = L8_2018_raw.normalizedDifference(['SR_B5','SR_B4']).rename('NDVI');
var ndbi  = L8_2018_raw.normalizedDifference(['SR_B6','SR_B5']).rename('NDBI');
var mndwi = L8_2018_raw.normalizedDifference(['SR_B3','SR_B6']).rename('MNDWI');

// ===== 纹理（GLCM）在 30 m 上计算，窗口 5×5，64 灰阶 =====
var texSrc = L8_2018_raw.select(['SR_B5','SR_B6'])
  .unitScale(0, 0.4).multiply(63).int16()    // 将反射率映射到 0..63
  .updateMask(L8_2018_raw.select(0).mask())
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

var glcm30 = texSrc.glcmTexture({size: 5});
var glcmSel30m = glcm30.select([
  'SR_B5_contrast','SR_B5_idm',
  'SR_B6_contrast','SR_B6_idm'
]).reproject(proj3857_30).setDefaultProjection(proj3857_30);

// ===== L8 特征栈（30 m）=====
var L8_feat_30m = L8_2018_raw.addBands([ndvi, ndbi, mndwi]).addBands(glcmSel30m)
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

var l8Mask = L8_feat_30m.select(0).mask();

// ====================== 3) Sentinel-1（VV/VH）2018 年特征 ======================
var eps = ee.Image.constant(1e-6); // 防止 log/除零
var S1_2018 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(region)
  .filterDate('2018-01-01', '2019-01-01')
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))  // ← 关键：用 eq 而不是 .or()
  .select(['VV','VH']);


var S1_lin = S1_2018.median().clip(region);
var VV_lin = S1_lin.select('VV').max(eps);
var VH_lin = S1_lin.select('VH').max(eps);

var VV_dB = VV_lin.log10().multiply(10).rename('VV_dB');
var VH_dB = VH_lin.log10().multiply(10).rename('VH_dB');
var VV_VH_ratio  = VV_lin.divide(VH_lin).rename('VV_VH_ratio');   // 线性域比值
var VV_VH_dBdiff = VV_dB.subtract(VH_dB).rename('VV_VH_dBdiff');  // dB差

// 用固定常数进行 unmask（避免稀疏导致的采样空洞），并对齐到 30 m
var s1_fill = ee.Image.cat([
  ee.Image.constant(-20).rename('VV_dB'),
  ee.Image.constant(-25).rename('VH_dB'),
  ee.Image.constant(2.0).rename('VV_VH_ratio'),
  ee.Image.constant(5.0).rename('VV_VH_dBdiff')
]);
var S1_feats_30m = VV_dB.addBands([VH_dB, VV_VH_ratio, VV_VH_dBdiff])
  .unmask(s1_fill)
  .reproject(proj3857_30).setDefaultProjection(proj3857_30)
  .updateMask(l8Mask); // 与 L8 有效掩膜对齐

print('S1 feature bands (30m):', S1_feats_30m.bandNames());

// ====================== 4) 融合特征（L8 + 指数 + 纹理 + S1，全部 30 m） ======================
var Feat_30m = L8_feat_30m.addBands(S1_feats_30m)
  .reproject(proj3857_30).setDefaultProjection(proj3857_30);

var inputBands = Feat_30m.bandNames();
print('All feature bands (L8+指数+纹理+S1, 30m):', inputBands);

// ====================== 5) 人工200/类（20 m 内缩） + 高置信100/类：采样点 ======================
var labelBand = 'Label';
var classes = ee.List.sequence(1, 17);
var N_MANUAL_PER_CLASS = 200;
var N_HIGHCONF_PER_CLASS = 100;
var SEED = 42;
var SHRINK_M = 20;   // 轻度向内收缩（米）

// 人工标注面（轻度内缩） → 随机点
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

// 高置信分层（100 m）
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

// 合并两批点
var allPts = manualPts.merge(highconfPts);
print('Total labeled points (manual+highconf):', allPts.size());
Map.addLayer(allPts.style({pointSize: 2, color: 'cyan'}), {}, '样本点(预览)', false);

// ====================== 6) 在 30 m 融合特征上提取样本、训练 RF ======================
var samplesAll = Feat_30m.sampleRegions({
  collection: allPts,
  properties: [labelBand],
  scale: 30,
  tileScale: 2,
  geometries: false
}).filter(ee.Filter.notNull([labelBand]));
print('Samples with (L8+指数+纹理+S1, 30m) features:', samplesAll.size());

var TRAIN_FRAC = 0.7;
var withRand = samplesAll.randomColumn('rand', SEED);
var trainFC  = withRand.filter(ee.Filter.lt('rand', TRAIN_FRAC));
var testFC   = withRand.filter(ee.Filter.gte('rand', TRAIN_FRAC));
print('Train size:', trainFC.size(), 'Test size:', testFC.size());

var RF_TREES = 200;
var rf = ee.Classifier.smileRandomForest({numberOfTrees: RF_TREES, seed: SEED})
  .train({features: trainFC, classProperty: labelBand, inputProperties: inputBands});

// ====================== 7) 留出法评估 ======================
var validated = testFC.classify(rf);
var cm = validated.errorMatrix(labelBand, 'classification');
print('Validation confusion matrix (L8+指数+纹理+S1, 30m)', cm);
print('Validation OA (L8+指数+纹理+S1, 30m)', cm.accuracy());
print('Validation Kappa (L8+指数+纹理+S1, 30m)', cm.kappa());
print('UA per class (L8+指数+纹理+S1, 30m)', cm.consumersAccuracy());
print('PA per class (L8+指数+纹理+S1, 30m)', cm.producersAccuracy());

// ====================== 8) 30 m 推理输出 ======================
var lcz2018_pred_30m = Feat_30m.classify(rf)
  .rename('LCZ_2018_pred_L8IdxTexS1')
  .toInt16()
  .clip(region);
Map.addLayer(lcz2018_pred_30m, lczVis, 'LCZ 2018 预测 (L8+指数+纹理+S1, 30 m)', true);

// ====================== 9) （可选）导出 ======================
// Export.image.toAsset({
//   image: lcz2018_pred_30m,
//   description: 'LCZ_2018_L8IdxTexS1_RF_30m',
//   assetId: 'users/你的用户名/LCZ_2018_L8IdxTexS1_RF_30m',
//   region: region,
//   scale: 30,
//   maxPixels: 1e13
// });
