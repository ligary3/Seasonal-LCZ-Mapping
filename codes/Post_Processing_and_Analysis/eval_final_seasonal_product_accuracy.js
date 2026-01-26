/****************************************************************
 * 【V304 Final Stable Custom】V302引擎 + V303逻辑
 * 1. 引擎：回滚到 V302 (Strict 1D)，解决 Error code 3 维度报错。
 * 2. 逻辑：保留 V303 的定制分类 (LCZ 15=Built, LCZ 9=Natural)。
 ****************************************************************/

// ==============================================================
// 1. 初始化
// ==============================================================
var ROI_ASSET = "projects/ee-l2892786691/assets/JJJ/jjj";
var LABEL_ASSET = "projects/ee-l2892786691/assets/JJJ/2018seasonlabel";
var REF_IMG_ASSET = "projects/ee-l2892786691/assets/JJJ/2018_ref";

var jjj_vec = ee.FeatureCollection(ROI_ASSET);
var baselineRef = ee.Image(REF_IMG_ASSET);
var seasonal_samples = ee.FeatureCollection(LABEL_ASSET);

Map.centerObject(jjj_vec, 7);
var ROI = jjj_vec.geometry();

var SEASONS_ORDER = ['SPRING', 'SUMMER', 'AUTUMN', 'WINTER'];
var NUM_CLASSES = 17; 
var EPSILON = 1e-6; 
var POINTS_PER_CLASS = 500; 

// 【定制逻辑】
function isBuiltIndex(idx) {//注意是减1后的值
  var i = ee.Number(idx);
  // Built: 0-7(LCZ1-8), 9(LCZ10), 14(LCZ15/E)
  return (i.lte(7)).or(i.eq(9)).or(i.eq(14));
}

function getBuiltMask(img) {//真实LCZ类别。
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
function fillGaps(img, ann) { return img.unmask(ann).unmask(0); }//年均填充，不足的用0填充。

function getFeatures_Full(year, seasonName) {
  var yStr = ee.Number(year).format('%d');
  var nextYStr = ee.Number(year).add(1).format('%d');
  var startYear = ee.String(yStr).cat('-01-01');
  var endYear = ee.String(nextYStr).cat('-01-01');
  var dates = YEARLY_DATES[year] || YEARLY_DATES['2018'];
  var range = ee.Dictionary(ee.Dictionary(dates).get(seasonName));
  var s2Bands = ['B2','B3','B4','B8','B11','B12'];

  var s2Col = S2.filterBounds(ROI).filterDate(range.get('start'), range.get('end')).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80)).map(maskS2);
  var s2_median = s2Col.median();//对应年份季节的中位数图像
  var s2_safe = ee.Algorithms.If(s2_median.bandNames().size().gt(0), s2_median.select(s2Bands), ee.Image.constant([0,0,0,0,0,0]).rename(s2Bands).toFloat());//无数据时用0填充
  var s2_sea = ee.Image(s2_safe);
  
  var s2AnnCol = S2.filterBounds(ROI).filterDate(startYear, endYear).filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80)).map(maskS2);
  var s2_ann = ee.Algorithms.If(s2AnnCol.median().bandNames().size().gt(0), s2AnnCol.median().select(s2Bands), ee.Image.constant([0,0,0,0,0,0]).rename(s2Bands).toFloat());
  s2_ann = ee.Image(s2_ann);//对应年份的年均图像

  var s2_fin = addIndices(fillGaps(s2_sea, s2_ann)).select(['B2','B3','B4','B8','B11','B12','NDVI','MNDWI','NDBI']);
  var s2_base = addIndices(fillGaps(s2_ann, s2_sea)).select(['B2','B3','B4','B8','B11','B12','NDVI','MNDWI','NDBI']);//用季节中位数填充年均，保证有数据。
  
  var deltaIdx = s2_fin.select(['NDVI','MNDWI','NDBI']).subtract(s2_base.select(['NDVI','MNDWI','NDBI'])).rename(['Delta_NDVI','Delta_MNDWI','Delta_NDBI']);

  //同样的处理 Sentinel-1，SAR
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
var labelImg = empty.paint(samplesProcessed, 'label_sp').rename('label_sp')//四个波段图像
  .addBands(empty.paint(samplesProcessed, 'label_su').rename('label_su'))
  .addBands(empty.paint(samplesProcessed, 'label_au').rename('label_au'))
  .addBands(empty.paint(samplesProcessed, 'label_wi').rename('label_wi'))
  .int32();

var stratifiedPoints = labelImg.stratifiedSample({//分层采样，每个类采样POINTS_PER_CLASS个点，stratifiedSample 会采集 labelImg 这个影像中包含的所有波段
  numPoints: POINTS_PER_CLASS, classBand: 'label_sp', region: samplesProcessed.geometry(), scale: 30, geometries: true, dropNulls: true, seed: 42
});//// label_sp告诉GEE用哪个波段来做分层（分类依据）

var pointsWithRandom = stratifiedPoints.randomColumn('random', 42);//添加随机列，用于划分训练测试集0-1
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

function exportMatrix(matrix, name) {
  var colNames = ee.List.sequence(1, 17).map(function(i){ return ee.String('To_LCZ_').cat(ee.Number(i).format('%d')) });
  var rows = ee.List.sequence(0, 16).map(function(r){
    var rowVals = matrix.exp().slice(0, r, ee.Number(r).add(1)).project([1]).toList();
    var dict = ee.Dictionary.fromLists(colNames, rowVals).set('From_LCZ', ee.Number(r).add(1));
    return ee.Feature(null, dict);
  });
  Export.table.toDrive({collection: ee.FeatureCollection(rows), description: name, fileFormat: 'CSV', folder: 'LCZ_V304_Stats'});
}
// exportMatrix(T_Sp_Sm, 'TransMatrix_Spring_Summer');
// exportMatrix(T_Sm_Au, 'TransMatrix_Summer_Autumn');
// exportMatrix(T_Au_Wi, 'TransMatrix_Autumn_Winter');
// exportMatrix(T_Wi_Sp, 'TransMatrix_Winter_Spring');

// ==============================================================
// 5. 训练
// ==============================================================
var trainedClassifiers = {};
var accuracyMetrics = [];

print('---- Training ----');
SEASONS_ORDER.forEach(function(season){//训练2018各季节的RF
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
//Export.table.toDrive({collection: ee.FeatureCollection(accuracyMetrics), description: 'LCZ_Accuracy_Report', fileFormat: 'CSV', folder: 'LCZ_V304_Stats'});

// ==============================================================
// 6. HMM 引擎 (回滚 V302 引擎 - Strict 1D)
// ==============================================================
function runHMM_Final(year) {
  var obsLogs = SEASONS_ORDER.map(function(season) {//计算观测概率，也就是发射概率
    var feat = getFeatures_Full(year, season);
    var rawProb = feat.classify(trainedClassifiers[season]);
    var learnedClasses = ee.List(trainedClassifiers[season].explain().get('classes'));
    var rawBandNames = learnedClasses.map(function(n){ return ee.String('P_').cat(ee.Number(n).format('%d')) });
    var flatProb = rawProb.arrayFlatten([rawBandNames]);
    var fullBands = ee.List.sequence(0, 16).map(function(i){
      var targetName = ee.String('P_').cat(ee.Number(i).format('%d'));
      return ee.Image(ee.Algorithms.If(flatProb.bandNames().contains(targetName), flatProb.select(targetName), ee.Image(1e-6))).unmask(1e-6);
    });
    // [Fix]: Strictly 1D [17]
    return ee.ImageCollection(fullBands).toBands().toArray().arrayProject([0]).max(1e-6).log();
  });

  /**
   * Performs the Viterbi step for Hidden Markov Model (HMM) calculations.
   * 
   * @param {ee.Image} prev - Previous state probabilities (1D array of length 17)
   * @param {ee.Image} trans - Transition probability matrix (17x17)
   * @param {ee.Image} emiss - Emission probabilities for current observation (1D array of length 17)
   * @returns {Object} An object containing:
   *   - dp: Updated state probabilities (1D array of length 17)
   *   - bp: Backpointer indices for best previous state (1D array of length 17)
   */
  function vStep(prev, trans, emiss) {
    var transImg = ee.Image(trans);
    
    // 1D [17] -> 2D [17, 17] for calculation
    var pathProb = prev.arrayReshape(ee.Image(ee.Array([17, 1])), 2).arrayRepeat(1, 17).add(transImg);
    var maxVal2D = pathProb.arrayReduce(ee.Reducer.max(), [0]);
    
    var idxM = ee.Image(ee.Array(ee.List.sequence(0, 16)).reshape([17, 1])).arrayRepeat(1, 17);
    var isMax = pathProb.eq(maxVal2D.arrayRepeat(0, 17));
    var bestPrevIndex2D = idxM.multiply(isMax).arrayReduce(ee.Reducer.max(), [0]);
    
    // Project back to 1D [17] IMMEDIATELY
    var dp = maxVal2D.arrayProject([1]).add(emiss); 
    var bp = bestPrevIndex2D.arrayProject([1]);     
    return {dp: dp, bp: bp};
  }

  var s0 = obsLogs[0]; // 1D
  var s1 = vStep(s0, T_Sp_Sm, obsLogs[1]); 
  var s2 = vStep(s1.dp, T_Sm_Au, obsLogs[2]); 
  var s3 = vStep(s2.dp, T_Au_Wi, obsLogs[3]); 
  
  // 1D argmax
  var w_idx = s3.dp.arrayArgmax().arrayGet([0]).toInt(); 

  // Flatten Backtrace
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
// 7. 年际后处理 (双向时空滤波 - 允许拆迁，过滤噪声)
// ==============================================================
var rawResults = {};
// 注意：确保 years 数组包含你所有需要处理的年份
var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

// 运行 HMM 获取原始结果
years.forEach(function(y){
  rawResults[y] = runHMM_Final(y);
});

// 定义存储修正结果的字典
var refined_Au = {};

// 辅助函数：判断是否为自然地表 (非建筑)
function isNature(img) { return getBuiltMask(img).not(); }
function isBuilt(img) { return getBuiltMask(img); }

// --- 步骤 A: 确定基准年份 (2018) ---
// 2018 作为基准年，不做时空修正，直接使用原始结果
refined_Au['2018'] = rawResults[2018].AUTUMN;

// --- 步骤 B: 2017年 (逆推) ---
// 由于没有 2016 数据参考，我们无法判断 2017 的突变是噪声还是真实变化。
// 策略：保守处理，信任 2017 原始数据，或者仅做非常轻微的逻辑修正。
// 这里选择：信任原始数据 (允许 2017 是建筑，2018 变成了自然，即允许拆迁)
refined_Au['2017'] = rawResults[2017].AUTUMN; 


// --- 步骤 C: 2019 - 2023 (三期滑动窗口逻辑) ---
// 核心逻辑：检查 T-1, T, T+1 三年
var years_forward = [2019, 2020, 2021, 2022, 2023];

years_forward.forEach(function(y) {
  var prev = refined_Au[(y - 1).toString()]; // 上一年 (已修正)
  var curr = rawResults[y].AUTUMN;           // 今年 (原始)
  var next = rawResults[y + 1].AUTUMN;       // 下一年 (原始 - 用于验证)
  
  // 1. 默认底图：使用今年的原始分类
  var finalImg = curr;

  // 2. 过滤“假扩张”噪声 (Nature -> Built -> Nature)
  // 如果去年是自然，今年突变成建筑，明年又变回自然 -> 认为是噪声，修回自然
  // 逻辑：Prev=Nature AND Curr=Built AND Next=Nature ==> Fix to Prev
  var isFalseExpansion = isNature(prev).and(isBuilt(curr)).and(isNature(next));
  finalImg = finalImg.where(isFalseExpansion, prev);

  // 3. 过滤“假拆迁”噪声 (Built -> Nature -> Built)
  // 如果去年是建筑，今年突变成自然，明年又变回建筑 -> 认为是噪声，修回建筑
  // 逻辑：Prev=Built AND Curr=Nature AND Next=Built ==> Fix to Prev
  var isFalseDemolition = isBuilt(prev).and(isNature(curr)).and(isBuilt(next));
  finalImg = finalImg.where(isFalseDemolition, prev);

  // --- 关键点解释 ---
  // 情况：去年=建筑(Built), 今年=自然(Nature), 明年=自然(Nature)
  // 此时 "isFalseDemolition" 判断为 false (因为 Next 不是 Built)。
  // 所以 where 函数不执行，finalImg 保持为 curr (Nature)。
  // 结果：成功保留了真实的拆迁变化！
  
  refined_Au[y.toString()] = finalImg;
});

// --- 步骤 D: 2024年 (边界处理) ---
// 2024 没有下一年来验证，无法区分“真拆迁”还是“假拆迁”。
// 策略：信任原始数据 (或者你可以选择严格继承 2023，取决于你觉得 2024 噪声大不大)
// 这里选择：信任原始数据，但防止极其离谱的突变 (可选)
refined_Au['2024'] = rawResults[2024].AUTUMN;


// ==============================================================
// 8. 导出
// ==============================================================
function fuse(hmm, ref) { return hmm.where(getBuiltMask(ref), ref); }

var exportFolder = 'LCZ_V304_Product';
var exportScale = 10;

// years.forEach(function(y) {
//   var raw = rawResults[y];
//   var refinedAutumn = refined_Au[y.toString()];
//   var isAuBuilt = getBuiltMask(refinedAutumn);
  
//   var sOut = raw.SPRING.where(isAuBuilt, refinedAutumn);
//   var smOut = raw.SUMMER.where(isAuBuilt, refinedAutumn);
//   var wOut = raw.WINTER.where(isAuBuilt, refinedAutumn);
//   var aOut = refinedAutumn; 
  
//   // if (y === 2018) {//这个可删除了，只有2018的不对
//   //   sOut = fuse(sOut, baselineRef);
//   //   smOut = fuse(smOut, baselineRef);
//   //   aOut = fuse(aOut, baselineRef);
//   //   wOut = fuse(wOut, baselineRef);
//   // }

//   sOut = sOut.clip(jjj_vec.geometry());
//   smOut = smOut.clip(jjj_vec.geometry());
//   aOut = aOut.clip(jjj_vec.geometry());
//   wOut = wOut.clip(jjj_vec.geometry());

//   Export.image.toDrive({image: sOut, description: y + '_Spring_LCZ', region: jjj_vec.geometry(), scale: exportScale, maxPixels: 1e13, fileFormat: 'GeoTIFF', folder: exportFolder});
//   Export.image.toDrive({image: smOut, description: y + '_Summer_LCZ', region: jjj_vec.geometry(), scale: exportScale, maxPixels: 1e13, fileFormat: 'GeoTIFF', folder: exportFolder});
//   Export.image.toDrive({image: aOut, description: y + '_Autumn_LCZ', region: jjj_vec.geometry(), scale: exportScale, maxPixels: 1e13, fileFormat: 'GeoTIFF', folder: exportFolder});
//   Export.image.toDrive({image: wOut, description: y + '_Winter_LCZ', region: jjj_vec.geometry(), scale: exportScale, maxPixels: 1e13, fileFormat: 'GeoTIFF', folder: exportFolder});
// });

print("V304: V302 Engine + V303 Custom Logic. Confirmed.");
// ==============================================================
// ==============================================================
// 9. 最终产品精度验证 (2018年) - [修复 .toList() 报错版]
// ==============================================================
print('==== 最终产品精度评估 (2018 Final Product Validation) ====');

// --- 1. 重建 2018 最终产品 ---
var raw2018 = rawResults[2018];
var final_Au_2018 = rawResults[2018].AUTUMN; 
var isAuBuilt_2018 = getBuiltMask(final_Au_2018);

var finalProducts2018 = {
  'SPRING': raw2018.SPRING.where(isAuBuilt_2018, final_Au_2018),
  'SUMMER': raw2018.SUMMER.where(isAuBuilt_2018, final_Au_2018),
  'AUTUMN': final_Au_2018,
  'WINTER': raw2018.WINTER.where(isAuBuilt_2018, final_Au_2018)
};

// --- 准备工作 ---
var fullValidationFeatures = []; 
var classIndices = ee.List.sequence(0, 16); // 0-16 代表 LCZ 1-17
// 准备矩阵的列名: [To_LCZ_1, ... To_LCZ_17]
var matrixColNames = classIndices.map(function(i){ 
  return ee.String('To_LCZ_').cat(ee.Number(i).add(1).format('%d')) 
});

var seasonConfig = [
  {name: 'SPRING', labelCol: 'label_sp'},
  {name: 'SUMMER', labelCol: 'label_su'},
  {name: 'AUTUMN', labelCol: 'label_au'},
  {name: 'WINTER', labelCol: 'label_wi'}
];

// --- 2. 循环处理每个季节 ---
seasonConfig.forEach(function(cfg){
  var seasonName = cfg.name;
  var labelCol = cfg.labelCol;
  
  // 预测结果 (1-17) -> 索引 (0-16)
  var predictedImg = finalProducts2018[seasonName].subtract(1).rename('prediction').toInt();
  
  // 采样
  var validation = predictedImg.sampleRegions({
    collection: testSet,
    properties: [labelCol],
    scale: 30,
    tileScale: 16,
    geometries: false
  });
  
  // 计算混淆矩阵
  var errorM = validation.errorMatrix(labelCol, 'prediction', classIndices);
  
  // ==========================================
  //  Part A: 计算汇总指标 (OA, Kappa, PA, UA)
  // ==========================================
  var oa = errorM.accuracy();
  var kappa = errorM.kappa();
  var paList = errorM.producersAccuracy().project([0]).toList(); 
  var uaList = errorM.consumersAccuracy().project([1]).toList(); 
  
  var stats = ee.Dictionary({
    'Season': seasonName,
    'OA': oa,
    'Kappa': kappa
  });

  var classStats = classIndices.map(function(idx){
    var i = ee.Number(idx);
    var lczNum = i.add(1).format('%d');
    return [
      ee.String('PA_LCZ_').cat(lczNum), paList.get(i),
      ee.String('UA_LCZ_').cat(lczNum), uaList.get(i)
    ];
  }).flatten();
  
  fullValidationFeatures.push(ee.Feature(null, stats.combine(ee.Dictionary(classStats))));
  
  print('--- ' + seasonName + ' ---', 'OA:', oa, 'Kappa:', kappa);

  // ==========================================
  //  Part B: 导出完整混淆矩阵 (修复报错部分)
  // ==========================================
  // 【修复点】：先转 array() 再转 toList()
  var matrixList = errorM.array().toList(); 
  
  // 使用 map 索引的方式构建 Feature，比 indexOf 更稳定
  var matrixRows = classIndices.map(function(i){
    var idx = ee.Number(i);
    var rowList = ee.List(matrixList.get(idx)); // 获取第i行
    var rowDict = ee.Dictionary.fromLists(matrixColNames, rowList);
    // 设置 From_LCZ 属性
    return ee.Feature(null, rowDict).set('From_LCZ', idx.add(1));
  });

  // 导出该季节的矩阵 CSV
  Export.table.toDrive({
    collection: ee.FeatureCollection(matrixRows),
    description: 'Matrix_2018_' + seasonName,
    fileFormat: 'CSV',
    folder: 'LCZ_V304_Validation'
  });
});

// --- 3. 导出汇总指标 CSV ---
Export.table.toDrive({
  collection: ee.FeatureCollection(fullValidationFeatures),
  description: 'Summary_Metrics_2018_Full',
  fileFormat: 'CSV',
  folder: 'LCZ_V304_Validation'
});