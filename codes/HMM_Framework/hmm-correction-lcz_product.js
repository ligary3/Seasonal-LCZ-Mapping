/****************************************************************
 * 【V306 Final Perfect Logic - Anti-Noise & Anti-FalseWater】
 * 1. 引擎：HMM V302 Strict 1D (稳定版)
 * 2. 逻辑：双层时空滤波 (Two-stage Temporal Filter)
 * 3. 2017逆推：增加了“去伪水体”逻辑 (Anti-FalseWater)
 ****************************************************************/

// ==============================================================
// 1. 初始化
// ==============================================================
var ROI_ASSET = "projects/ee-l2892786691/assets/JJJ/jjj";
var LABEL_ASSET = "projects/ee-l2892786691/assets/JJJ/2018seasonlabel";
var REF_IMG_ASSET = "projects/ee-l2892786691/assets/JJJ/2018_ref";

var jjj_vec = ee.FeatureCollection(ROI_ASSET);
//var baselineRef = ee.Image(REF_IMG_ASSET);
var seasonal_samples = ee.FeatureCollection(LABEL_ASSET);

Map.centerObject(jjj_vec, 7);
var ROI = jjj_vec.geometry();

var SEASONS_ORDER = ['SPRING', 'SUMMER', 'AUTUMN', 'WINTER'];
var NUM_CLASSES = 17; 
var EPSILON = 1e-6; 
var POINTS_PER_CLASS = 500; 

// 【定制逻辑】
function isBuiltIndex(idx) { // 注意是减1后的值
  var i = ee.Number(idx);
  // Built: 0-7(LCZ1-8), 9(LCZ10), 14(LCZ15/E)
  return (i.lte(7)).or(i.eq(9)).or(i.eq(14));
}

function getBuiltMask(img) { // 真实LCZ类别
  // Built: 1-8, 10, 15
  return img.lte(8).or(img.eq(10)).or(img.eq(15));
}

var YEARLY_DATES = {
  '2017': { 'SPRING': {start: '2017-04-04', end: '2017-05-18'}, 'SUMMER': {start: '2017-05-19', end: '2017-08-26'}, 'AUTUMN': {start: '2017-08-27', end: '2017-10-11'}, 'WINTER': {start: '2017-10-12', end: '2018-03-25'} },
  '2018': { 'SPRING': {start: '2018-03-26', end: '2018-06-01'}, 'SUMMER': {start: '2018-06-02', end: '2018-08-30'}, 'AUTUMN': {start: '2018-08-31', end: '2018-10-09'}, 'WINTER': {start: '2018-10-10', end: '2019-03-19'} },
  '2019': { 'SPRING': {start: '2019-03-20', end: '2019-05-24'}, 'SUMMER': {start: '2019-05-25', end: '2019-08-16'}, 'AUTUMN': {start: '2019-08-17', end: '2019-10-14'}, 'WINTER': {start: '2019-10-15', end: '2020-03-24'} },
  '2020': { 'SPRING': {start: '2020-03-25', end: '2020-06-05'}, 'SUMMER': {start: '2020-06-06', end: '2020-08-18'}, 'AUTUMN': {start: '2020-08-19', end: '2020-10-14'}, 'WINTER': {start: '2020-10-15', end: '2021-03-26'} },
  '2021': { 'SPRING': {start: '2021-03-27', end: '2021-06-11'}, 'SUMMER': {start: '2021-06-12', end: '2021-08-14'}, 'AUTUMN': {start: '2021-08-15', end: '2021-10-15'}, 'WINTER': {start: '2021-10-16', end: '2022-04-06'} },
  '2022': { 'SPRING': {start: '2022-04-07', end: '2022-05-22'}, 'SUMMER': {start: '2022-05-23', end: '2022-08-23'}, 'AUTUMN': {start: '2022-08-24', end: '2022-10-05'}, 'WINTER': {start: '2022-10-06', end: '2023-03-09'} },
  '2023': { 'SPRING': {start: '2023-03-10', end: '2023-06-06'}, 'SUMMER': {start: '2023-06-07', end: '2023-08-24'}, 'AUTUMN': {start: '2023-08-25', end: '2023-11-05'}, 'WINTER': {start: '2023-11-06', end: '2024-04-05'} },
  '2024': { 'SPRING': {start: '2024-04-06', end: '2024-06-04'}, 'SUMMER': {start: '2024-06-05', end: '2024-08-28'}, 'AUTUMN': {start: '2024-08-29', end: '2024-10-19'}, 'WINTER': {start: '2024-10-20', end: '2025-03-22'} }
};

// ==============================================================
// 2. 特征工程
// ==============================================================
var S2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var S1 = ee.ImageCollection('COPERNICUS/S1_GRD');
var AEF_COL = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL');

function maskS2(img){ return img.updateMask(img.select('QA60').bitwiseAnd(1<<10).eq(0).and(img.select('QA60').bitwiseAnd(1<<11).eq(0))).divide(10000); }
function addIndices(img){ return img.addBands([img.normalizedDifference(['B8','B4']).rename('NDVI'), img.normalizedDifference(['B3','B11']).rename('MNDWI'), img.normalizedDifference(['B11','B8']).rename('NDBI')]); }
function fillGaps(img, ann) { return img.unmask(ann).unmask(0); }

function getFeatures_Full(year, seasonName) {
  var yStr = ee.Number(year).format('%d');
  var nextYStr = ee.Number(year).add(1).format('%d');
  var startYear = ee.String(yStr).cat('-01-01');
  var endYear = ee.String(nextYStr).cat('-01-01');
  var dates = YEARLY_DATES[year] || YEARLY_DATES['2018'];
  var range = ee.Dictionary(ee.Dictionary(dates).get(seasonName));
  var s2Bands = ['B2','B3','B4','B8','B11','B12'];

  var s2Col = S2.filterBounds(ROI).filterDate(range.get('start'), range.get('end')).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80)).map(maskS2);
  var s2_median = s2Col.median();
  var s2_safe = ee.Algorithms.If(s2_median.bandNames().size().gt(0), s2_median.select(s2Bands), ee.Image.constant([0,0,0,0,0,0]).rename(s2Bands).toFloat());
  var s2_sea = ee.Image(s2_safe);
   
  var s2AnnCol = S2.filterBounds(ROI).filterDate(startYear, endYear).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80)).map(maskS2);
  var s2_ann = ee.Algorithms.If(s2AnnCol.median().bandNames().size().gt(0), s2AnnCol.median().select(s2Bands), ee.Image.constant([0,0,0,0,0,0]).rename(s2Bands).toFloat());
  s2_ann = ee.Image(s2_ann);

  var s2_fin = addIndices(fillGaps(s2_sea, s2_ann)).select(['B2','B3','B4','B8','B11','B12','NDVI','MNDWI','NDBI']);
  var s2_base = addIndices(fillGaps(s2_ann, s2_sea)).select(['B2','B3','B4','B8','B11','B12','NDVI','MNDWI','NDBI']);
   
  var deltaIdx = s2_fin.select(['NDVI','MNDWI','NDBI']).subtract(s2_base.select(['NDVI','MNDWI','NDBI'])).rename(['Delta_NDVI','Delta_MNDWI','Delta_NDBI']);

  var s1Col = S1.filterBounds(ROI).filterDate(range.get('start'), range.get('end')).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'));
  var s1_sea = s1Col.select(['VV','VH']).mean();
  var s1_safe = ee.Algorithms.If(s1_sea.bandNames().size().gt(0), s1_sea, ee.Image.constant([0,0]).rename(['VV','VH']).toFloat());
  var s1_fin = ee.Image(s1_safe).rename(['VV','VH']);

  var s1AnnCol = S1.filterBounds(ROI).filterDate(startYear, endYear).filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'));
  var s1_ann = s1AnnCol.select(['VV','VH']).mean();
  var s1_ann_safe = ee.Algorithms.If(s1_ann.bandNames().size().gt(0), s1_ann, ee.Image.constant([0,0]).rename(['VV','VH']).toFloat());
  s1_ann = ee.Image(s1_ann_safe);

  var deltaSar = s1_fin.subtract(s1_ann).rename(['Delta_VV','Delta_VH']);

  var aefColYear = AEF_COL.filterDate(startYear, endYear);
  var aef = ee.Algorithms.If(aefColYear.size().gt(0), aefColYear.filterBounds(ROI).mosaic().unmask(0), AEF_COL.filterDate('2018-01-01', '2019-01-01').filterBounds(ROI).mosaic().unmask(0));

  return s2_fin.addBands(s1_fin).addBands(deltaIdx).addBands(deltaSar).addBands(ee.Image(aef)).clip(ROI).toFloat();
}

// ==============================================================
// 3. 采样 (Stratified)
// ==============================================================
var samplesProcessed = seasonal_samples.map(function(f){
  var sp = f.get('spring') || f.get('SPRING');
  var sm = f.get('summer') || f.get('SUMMER');
  var au = f.get('autumn') || f.get('AUTUMN');
  var wi = f.get('winter') || f.get('WINTER');
  return f.set({
    'label_sp': ee.Number(sp).subtract(1).toInt(),
    'label_su': ee.Number(sm).subtract(1).toInt(), 
    'label_au': ee.Number(au).subtract(1).toInt(),
    'label_wi': ee.Number(wi).subtract(1).toInt()
  });
}).filter(ee.Filter.and(ee.Filter.notNull(['label_sp']), ee.Filter.notNull(['label_su']), ee.Filter.notNull(['label_au']), ee.Filter.notNull(['label_wi'])));

var empty = ee.Image().byte();
var labelImg = empty.paint(samplesProcessed, 'label_sp').rename('label_sp')
  .addBands(empty.paint(samplesProcessed, 'label_su').rename('label_su'))
  .addBands(empty.paint(samplesProcessed, 'label_au').rename('label_au'))
  .addBands(empty.paint(samplesProcessed, 'label_wi').rename('label_wi'))
  .int32();

var stratifiedPoints = labelImg.stratifiedSample({
  numPoints: POINTS_PER_CLASS, classBand: 'label_sp', region: samplesProcessed.geometry(), scale: 30, geometries: true, dropNulls: true, seed: 42
});

var pointsWithRandom = stratifiedPoints.randomColumn('random', 42);
var trainSet = pointsWithRandom.filter(ee.Filter.lt('random', 0.7));
var testSet = pointsWithRandom.filter(ee.Filter.gte('random', 0.7));

print('Training:', trainSet.size(), 'Testing:', testSet.size());

// ==============================================================
// 4. 转移矩阵 (应用定制逻辑)
// ==============================================================
function computeTransitionMatrix(collection, fromCol, toCol) {
  var matrixList = ee.List.sequence(0, 16).map(function(fromClass) {
    var fromSubset = collection.filter(ee.Filter.eq(fromCol, fromClass));
    var counts = ee.List.sequence(0, 16).map(function(toClass) {
      // 【关键逻辑】：使用 isBuiltIndex
      var isBuilding = isBuiltIndex(fromClass); 
      var isDiagonal = ee.Number(fromClass).eq(toClass);
      var n = fromSubset.filter(ee.Filter.eq(toCol, toClass)).size();
      return ee.Algorithms.If(isBuilding, ee.Algorithms.If(isDiagonal, 10000, EPSILON), ee.Number(n).max(EPSILON));
    });
    return counts;
  });
  var matrix = ee.Array(matrixList);
  var rowSums = matrix.reduce(ee.Reducer.sum(), [1]);
  var safeSums = rowSums.add(rowSums.eq(0)); 
  return matrix.divide(safeSums.repeat(1, NUM_CLASSES)).max(EPSILON).log(); 
}

var T_Sp_Sm = computeTransitionMatrix(trainSet, 'label_sp', 'label_su');
var T_Sm_Au = computeTransitionMatrix(trainSet, 'label_su', 'label_au');
var T_Au_Wi = computeTransitionMatrix(trainSet, 'label_au', 'label_wi');
var T_Wi_Sp = computeTransitionMatrix(trainSet, 'label_wi', 'label_sp');

// ==============================================================
// 5. 训练
// ==============================================================
var trainedClassifiers = {};
var accuracyMetrics = [];

print('---- Training ----');
SEASONS_ORDER.forEach(function(season){
  var seasonLabel = 'label_' + season.substring(0, 2).toLowerCase();
  var feat2018 = getFeatures_Full(2018, season);
  var trainingData = feat2018.sampleRegions({collection: trainSet, properties: [seasonLabel], scale: 30, tileScale: 16, geometries: false});
  var rf = ee.Classifier.smileRandomForest({numberOfTrees: 200}).setOutputMode('MULTIPROBABILITY')
      .train({features: trainingData, classProperty: seasonLabel, inputProperties: feat2018.bandNames()});
  trainedClassifiers[season] = rf;
   
  var validationData = feat2018.sampleRegions({collection: testSet, properties: [seasonLabel], scale: 30, tileScale: 16, geometries: false});
  var rfVal = ee.Classifier.smileRandomForest({numberOfTrees: 200}).setOutputMode('CLASSIFICATION')
      .train({features: trainingData, classProperty: seasonLabel, inputProperties: feat2018.bandNames()});
  var validated = validationData.classify(rfVal);
  var acc = validated.errorMatrix(seasonLabel, 'classification');
  print(season, 'OA:', acc.accuracy(), 'Kappa:', acc.kappa());
  accuracyMetrics.push(ee.Feature(null, {'Season': season, 'OA': acc.accuracy(), 'Kappa': acc.kappa()}));
});

// ==============================================================
// 6. HMM 引擎 (Strict 1D)
// ==============================================================
function runHMM_Final(year) {
  var obsLogs = SEASONS_ORDER.map(function(season) {
    var feat = getFeatures_Full(year, season);
    var rawProb = feat.classify(trainedClassifiers[season]);
    var learnedClasses = ee.List(trainedClassifiers[season].explain().get('classes'));
    var rawBandNames = learnedClasses.map(function(n){ return ee.String('P_').cat(ee.Number(n).format('%d')) });
    var flatProb = rawProb.arrayFlatten([rawBandNames]);
    var fullBands = ee.List.sequence(0, 16).map(function(i){
      var targetName = ee.String('P_').cat(ee.Number(i).format('%d'));
      return ee.Image(ee.Algorithms.If(flatProb.bandNames().contains(targetName), flatProb.select(targetName), ee.Image(1e-6))).unmask(1e-6);
    });
    return ee.ImageCollection(fullBands).toBands().toArray().arrayProject([0]).max(1e-6).log();
  });

  function vStep(prev, trans, emiss) {
    var transImg = ee.Image(trans);
    var pathProb = prev.arrayReshape(ee.Image(ee.Array([17, 1])), 2).arrayRepeat(1, 17).add(transImg);
    var maxVal2D = pathProb.arrayReduce(ee.Reducer.max(), [0]);
    var idxM = ee.Image(ee.Array(ee.List.sequence(0, 16)).reshape([17, 1])).arrayRepeat(1, 17);
    var isMax = pathProb.eq(maxVal2D.arrayRepeat(0, 17));
    var bestPrevIndex2D = idxM.multiply(isMax).arrayReduce(ee.Reducer.max(), [0]);
    var dp = maxVal2D.arrayProject([1]).add(emiss); 
    var bp = bestPrevIndex2D.arrayProject([1]);      
    return {dp: dp, bp: bp};
  }

  var s0 = obsLogs[0]; 
  var s1 = vStep(s0, T_Sp_Sm, obsLogs[1]); 
  var s2 = vStep(s1.dp, T_Sm_Au, obsLogs[2]); 
  var s3 = vStep(s2.dp, T_Au_Wi, obsLogs[3]); 
   
  var w_idx = s3.dp.arrayArgmax().arrayGet([0]).toInt(); 

  function getP(bp, c) {
    var bpFlat = bp.arrayFlatten([ee.List.sequence(0, 16).map(function(i){ return ee.String('b').cat(ee.Number(i).int()) })]);
    var mask = ee.Image.constant(ee.List.sequence(0, 16)).eq(c);
    return bpFlat.multiply(mask).reduce(ee.Reducer.sum()).toInt();
  }
   
  var a = getP(s3.bp, w_idx), sm = getP(s2.bp, a), sp = getP(s1.bp, sm);
   
  return {
    SPRING: sp.add(1).toByte().rename('LCZ'), 
    SUMMER: sm.add(1).toByte().rename('LCZ'), 
    AUTUMN: a.add(1).toByte().rename('LCZ'), 
    WINTER: w_idx.add(1).toByte().rename('LCZ')
  };
}

// ==============================================================
// 7. 年际后处理 (Perfect Logic: Two-Stage Filter + Anti-FalseWater)
// ==============================================================
var rawResults = {};
var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

years.forEach(function(y){
  rawResults[y] = runHMM_Final(y);
});

var refined_Au = {};

// 辅助函数
function isNature(img) { return getBuiltMask(img).not(); }
function isBuilt(img) { return getBuiltMask(img); }
// 假设 LCZ G 是第 17 类 (根据 WUDAPT 标准)
function isWater(img) { return img.eq(17); } 

// --- 步骤 A: 2018 基准年 ---
refined_Au['2018'] = rawResults[2018].AUTUMN;

// --- 步骤 B: 2017年 (强逻辑逆推：抗噪+抗伪水体) ---
// 目标：解决“伪拆迁”和“植被误判为水(LCZG)”的问题

var raw17 = rawResults[2017].AUTUMN;
var ref18 = refined_Au['2018']; // 2018 是真值锚点

// 1. 初始化 2017 结果为 HMM 原始结果
var fixed17 = raw17;

// 修正 1: 禁止“伪拆迁” (Built -> Nature)
// 逻辑：如果 2018 是自然，2017 却算成了建筑，说明 2017 是噪声。
// (因为不可能 2017 有楼，2018 拆了变农田)
var isFalseDemolition17 = isNature(ref18).and(isBuilt(raw17));
fixed17 = fixed17.where(isFalseDemolition17, ref18);

// 修正 2: 禁止“伪水体” (Land -> Water -> Land) - 【解决你的新问题】
// 逻辑：如果 2018 不是水 (是陆地/植被/建筑)，但 2017 算成了水 (LCZ G)。
// 解释：植被阴影或暗色土壤极易被误判为水。湖泊变农田极少见。
// 操作：信任 2018 的陆地类别。
var isFalseWater = isWater(fixed17).and(isWater(ref18).not());
fixed17 = fixed17.where(isFalseWater, ref18);

// 修正 3: (可选) 保持建筑类别一致性
// 逻辑：如果 2017 和 2018 都是建筑，但类别不同 (如 1 -> 2)，保持稳定。
var isUnstableBuilt = isBuilt(ref18).and(isBuilt(fixed17)).and(ref18.neq(fixed17));
fixed17 = fixed17.where(isUnstableBuilt, ref18);

// 更新结果
refined_Au['2017'] = fixed17;


// --- 步骤 C: 2019 - 2023 滑动窗口修正 ---
var years_forward = [2019, 2020, 2021, 2022, 2023];

years_forward.forEach(function(y) {
  var prev = refined_Au[(y - 1).toString()]; // 上一年 (已修正)
  var curr = rawResults[y].AUTUMN;           // 今年 (原始)
  var next = rawResults[y + 1].AUTUMN;       // 下一年 (原始)
   
  // 1. 默认底图
  var finalImg = curr;

  // --- Layer 1: 大类物理约束 (Binary Physics) ---
  // 解决 "Nature -> Built -> Nature" (假扩张)
  var isFalseExpansion = isNature(prev).and(isBuilt(curr)).and(isNature(next));
  finalImg = finalImg.where(isFalseExpansion, prev);

  // 解决 "Built -> Nature -> Built" (假拆迁)
  var isFalseDemolition = isBuilt(prev).and(isNature(curr)).and(isBuilt(next));
  finalImg = finalImg.where(isFalseDemolition, prev);

  // --- Layer 2: 同类统计平滑 (Class-specific Sandwich Filter) ---
  // 解决 "LCZ 1 -> LCZ 2 -> LCZ 1" 这种同类内部抖动
  // 逻辑：如果 T-1 == T+1，但 T != T-1，则强制 T = T-1
  // 注意：此时 finalImg 已经经过了 Layer 1 处理，所以不会和 Layer 1 冲突
  var isSandwichNoise = prev.eq(next).and(finalImg.neq(prev));
  finalImg = finalImg.where(isSandwichNoise, prev);
   
  refined_Au[y.toString()] = finalImg;
});

// --- 步骤 D: 2024 边界处理 ---
refined_Au['2024'] = rawResults[2024].AUTUMN;


// ==============================================================
// 8. 导出
// ==============================================================
var exportYears = [2024]; // 【提示】建议一次只跑一年，避免内存溢出

var targetChunks = 8; 
var exportFolder = 'LCZ_V306_Product'; 
var exportScale = 10;

print("Preparing Grid...");
var combinedROI = jjj_vec.union(1).geometry(); 
var roiBounds = combinedROI.bounds(1).transform('EPSG:4326', 1);
var coords = ee.List(roiBounds.coordinates().get(0));
var xs = coords.map(function(p) { return ee.List(p).get(0); });
var ys = coords.map(function(p) { return ee.List(p).get(1); });
var minX = ee.Number(xs.reduce(ee.Reducer.min()));
var maxX = ee.Number(xs.reduce(ee.Reducer.max()));
var minY = ee.Number(ys.reduce(ee.Reducer.min()));
var maxY = ee.Number(ys.reduce(ee.Reducer.max()));

var width = maxX.subtract(minX).max(0.001);
var height = maxY.subtract(minY).max(0.001);
var ratio = width.divide(height);
var cols = ratio.sqrt().multiply(ee.Number(targetChunks).sqrt()).round().max(1);
var rows = ee.Number(targetChunks).divide(cols).ceil().max(1);

var latStep = height.divide(rows);
var lonStep = width.divide(cols);

var gridPolys = ee.List.sequence(0, rows.subtract(1)).map(function(r) {
  var y1 = minY.add(ee.Number(r).multiply(latStep));
  var y2 = y1.add(latStep);
  return ee.List.sequence(0, cols.subtract(1)).map(function(c) {
    var x1 = minX.add(ee.Number(c).multiply(lonStep));
    var x2 = x1.add(lonStep);
    return ee.Feature(ee.Geometry.Rectangle([x1, y1, x2, y2], 'EPSG:4326', false));
  });
}).flatten();

var grid = ee.FeatureCollection(gridPolys);
var validGrid = grid.filterBounds(combinedROI);
var gridList = validGrid.toList(validGrid.size());
var numTilesClient = gridList.size().getInfo(); 

print('Total Export Tiles:', numTilesClient);

if (numTilesClient === 0) {
  validGrid = ee.FeatureCollection([ee.Feature(roiBounds)]);
  gridList = validGrid.toList(1);
  numTilesClient = 1;
}

exportYears.forEach(function(y) {
  if (!rawResults[y] || !refined_Au[y.toString()]) {
    print('Error: Year ' + y + ' results not found.');
    return;
  }
  var raw = rawResults[y];
  var refinedAutumn = refined_Au[y.toString()];
   
  var isAuBuilt = getBuiltMask(refinedAutumn);// autumn to other season
  var aOut = refinedAutumn; 
  var sOut = raw.SPRING.where(isAuBuilt, refinedAutumn);
  var smOut = raw.SUMMER.where(isAuBuilt, refinedAutumn);
  var wOut = raw.WINTER.where(isAuBuilt, refinedAutumn);
   
  aOut = aOut.clip(combinedROI);
  sOut = sOut.clip(combinedROI);
  smOut = smOut.clip(combinedROI);
  wOut = wOut.clip(combinedROI);

  for (var i = 0; i < numTilesClient; i++) {
    var tileFeat = ee.Feature(gridList.get(i));
    var tileGeom = tileFeat.geometry(); 
    var taskSuffix = '_' + y + '_Tile' + (i+1); 
    
    // Export.image.toDrive({
    //   image: aOut, 
    //   description: 'Autumn_LCZ' + taskSuffix, 
    //   region: tileGeom, 
    //   scale: exportScale, 
    //   maxPixels: 1e13, 
    //   folder: exportFolder,
    //   crs: 'EPSG:3857' 
    // });
       Export.image.toDrive({ image: sOut, description: 'Spring_LCZ' + taskSuffix, region: tileGeom, scale: exportScale, maxPixels: 1e13, folder: exportFolder, crs: 'EPSG:4326' });
    Export.image.toDrive({ image: smOut, description: 'Summer_LCZ' + taskSuffix, region: tileGeom, scale: exportScale, maxPixels: 1e13, folder: exportFolder, crs: 'EPSG:4326' });
    Export.image.toDrive({ image: wOut, description: 'Winter_LCZ' + taskSuffix, region: tileGeom, scale: exportScale, maxPixels: 1e13, folder: exportFolder, crs: 'EPSG:4326' });
  }
});