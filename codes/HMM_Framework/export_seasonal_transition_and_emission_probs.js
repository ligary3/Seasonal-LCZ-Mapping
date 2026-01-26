/****************************************************************
 * 【V304 Final Stable Custom + Validation Exports】
 * 1. 引擎：V302 (Strict 1D) 解决维度报错。
 * 2. 逻辑：保留 V303 定制分类 (LCZ 15=Built, LCZ 9=Natural)。
 * 3. 新增功能：导出论文图表所需的所有统计数据 (Step 9-11)。
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
exportMatrix(T_Sp_Sm, 'TransMatrix_Spring_Summer');
exportMatrix(T_Sm_Au, 'TransMatrix_Summer_Autumn');
exportMatrix(T_Au_Wi, 'TransMatrix_Autumn_Winter');
exportMatrix(T_Wi_Sp, 'TransMatrix_Winter_Spring');

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
Export.table.toDrive({collection: ee.FeatureCollection(accuracyMetrics), description: 'LCZ_Accuracy_Report', fileFormat: 'CSV', folder: 'LCZ_V304_Stats'});

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
// 7. 年际后处理 (使用 Custom Logic)
// ==============================================================
var rawResults = {};
var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

years.forEach(function(y){
  rawResults[y] = runHMM_Final(y);
});

var refined_Au = {};
refined_Au['2018'] = rawResults[2018].AUTUMN;

// 2017 逆推
var au17 = rawResults[2017].AUTUMN;
var au18 = refined_Au['2018'];
// 【应用定制逻辑】
var is17Built = getBuiltMask(au17);
var is18Nature = getBuiltMask(au18).not();
refined_Au['2017'] = au17.where(is17Built.and(is18Nature), au18);

// 2019-2024 顺推
refined_Au['2019'] = refined_Au['2018'].where(rawResults[2019].AUTUMN.neq(refined_Au['2018']).and(rawResults[2019].AUTUMN.eq(rawResults[2020].AUTUMN)), rawResults[2019].AUTUMN);
refined_Au['2020'] = refined_Au['2019'].where(rawResults[2020].AUTUMN.neq(refined_Au['2019']).and(rawResults[2020].AUTUMN.eq(rawResults[2021].AUTUMN)), rawResults[2020].AUTUMN);
refined_Au['2021'] = refined_Au['2020'].where(rawResults[2021].AUTUMN.neq(refined_Au['2020']).and(rawResults[2021].AUTUMN.eq(rawResults[2022].AUTUMN)), rawResults[2021].AUTUMN);
refined_Au['2022'] = refined_Au['2021'].where(rawResults[2022].AUTUMN.neq(refined_Au['2021']).and(rawResults[2022].AUTUMN.eq(rawResults[2023].AUTUMN)), rawResults[2022].AUTUMN);
refined_Au['2023'] = refined_Au['2022'].where(rawResults[2023].AUTUMN.neq(refined_Au['2022']).and(rawResults[2023].AUTUMN.eq(rawResults[2024].AUTUMN)), rawResults[2023].AUTUMN);
refined_Au['2024'] = rawResults[2024].AUTUMN; 

// ==============================================================
// 8. 导出地图产品
// ==============================================================
function fuse(hmm, ref) { return hmm.where(getBuiltMask(ref), ref); }

var exportFolder = 'LCZ_V304_Product';
var exportScale = 10;


// ==============================================================
// 9. [新增] 2.3.1 季节性 RF 性能详细评估 (Seasonal Metrics Export)
//     用于导出 OA, Kappa 以及每个类别的 PA, UA, F1
// ==============================================================
// ==============================================================
// 9. [终极稳健版] 2.3.1 季节性 RF 性能详细评估
//     修复策略：List.flatten() 自动适应任何维度，杜绝报错
// ==============================================================
print('正在计算季节性 RF 精度指标 (OA/PA/UA/F1)...');

// 使用客户端 map 遍历季节
var seasonalFeaturesList = SEASONS_ORDER.map(function(season) {
  var seasonLabel = 'label_' + season.substring(0, 2).toLowerCase();
  var feat = getFeatures_Full(2018, season);
  
  // 1. 准备验证数据
  var validationData = feat.sampleRegions({
    collection: testSet, 
    properties: [seasonLabel], 
    scale: 30, 
    tileScale: 16, 
    geometries: false
  });
  
  // 2. 原地重新训练 RF (防止模型对象丢失)
  var rf = ee.Classifier.smileRandomForest({numberOfTrees: 200})
      .setOutputMode('CLASSIFICATION')
      .train({
        features: feat.sampleRegions({
            collection: trainSet, 
            properties: [seasonLabel], 
            scale: 30, tileScale: 16, geometries: false
        }), 
        classProperty: seasonLabel, 
        inputProperties: feat.bandNames()
      });
      
  // 3. 分类与混淆矩阵
  var validated = validationData.classify(rf);
  var order = ee.List.sequence(0, 16).map(function(n) { return ee.Number(n).toInt(); });
  var matrix = validated.errorMatrix(seasonLabel, 'classification', order);
  
  // 4. 提取指标
  var oa = matrix.accuracy();
  var kappa = matrix.kappa();
  
  // [核心修复] 不管是 Nx1 还是 1xN，统统转 List 再 flatten
  // 结果必然是包含 17 个数字的一维列表
  var paList = matrix.producersAccuracy().toList().flatten(); 
  var uaList = matrix.consumersAccuracy().toList().flatten(); 
  
  // 5. 遍历每一类 (0-16)
  return order.map(function(i) {
    var classIdx = ee.Number(i);
    var classId = classIdx.add(1); // 1-17
    
    // 从扁平化列表中安全取值
    var pa = ee.Number(paList.get(classIdx));
    var ua = ee.Number(uaList.get(classIdx));
    
    // 计算 F1
    var f1 = ee.Number(2).multiply(pa).multiply(ua)
              .divide(pa.add(ua).max(1e-6));
    
    return ee.Feature(null, {
      'Season': season,
      'Class_ID': classId,
      'OA': oa,       
      'Kappa': kappa, 
      'PA': pa,
      'UA': ua,
      'F1': f1
    });
  });
});

// Flatten the feature list
var exportFeatures = ee.FeatureCollection(ee.List(seasonalFeaturesList).flatten());

// 6. 导出 CSV
Export.table.toDrive({
  collection: exportFeatures,
  description: 'Seasonal_RF_Performance_Metrics',
  fileFormat: 'CSV',
  folder: 'LCZ_Paper_Figs',
  selectors: ['Season', 'Class_ID', 'OA', 'Kappa', 'PA', 'UA', 'F1']
});


