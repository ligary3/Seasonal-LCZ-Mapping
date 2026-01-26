/**
 * ============================================================
 * LCZ 季节时间生成器 (适配 V304 代码结构)
 * ============================================================
 * 逻辑调整：
 * 1. 顺序：Spring -> Summer -> Autumn -> Winter
 * 2. Winter：从当年年底开始，延续到下一年春天开始前。
 * 3. 严格匹配你的 YEARLY_DATES 字典格式。
 */

// 1. 城市配置
var cityConfig = [
  // --- A. 标准自适应组 (JJJ, Berlin, NYC, Murmansk) ---
  {name: 'JJJ_Region', geom: ee.Geometry.Point([116.40, 39.90]), method: 'adaptive', t_spring: 10, t_summer: 22, t_winter: 10},
  {name: 'Berlin',     geom: ee.Geometry.Point([13.40, 52.52]),  method: 'adaptive', t_spring: 10, t_summer: 22, t_winter: 10},
  {name: 'New_York',   geom: ee.Geometry.Point([-74.00, 40.71]), method: 'adaptive', t_spring: 10, t_summer: 22, t_winter: 10},
  {name: 'Murmansk',   geom: ee.Geometry.Point([33.08, 68.97]),  method: 'adaptive', t_spring: 0,  t_summer: 10, t_winter: 0},
  
  // --- B. 固定组 (Hong Kong, Cairo, Singapore, Sao Paulo) ---
  // 统一按 WMO 3月1日起算，确保 Spring 在前，Winter 在最后跨年
  {name: 'Hong_Kong',  geom: ee.Geometry.Point([114.17, 22.32]), method: 'fixed'},
  {name: 'Cairo',      geom: ee.Geometry.Point([31.23, 30.04]),  method: 'fixed'},
  {name: 'Singapore',  geom: ee.Geometry.Point([103.81, 1.35]),  method: 'fixed'},
  {name: 'Sao_Paulo',  geom: ee.Geometry.Point([-46.63, -23.55]),method: 'fixed'} 
];

var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR").select('temperature_2m');
var years = ee.List.sequence(2017, 2024);

// 5日滑动平均
var applyMovingAverage = function(fc) {
  var dayWindow = 2.5 * 24 * 3600 * 1000; 
  var filter = ee.Filter.maxDifference({difference: dayWindow, leftField: 'millis', rightField: 'millis'});
  var join = ee.Join.saveAll({matchesKey: 'window'});
  var joined = join.apply(fc, fc, filter);
  return joined.map(function(f) {
    var win = ee.List(f.get('window'));
    var winFC = ee.FeatureCollection(win);
    return f.set('temp_5d', winFC.aggregate_mean('temp'));
  });
};

// 核心逻辑
var getDatesForCode = function(cfg, year) {
  var yStr = ee.Number(year).format('%d');
  var nextYStr = ee.Number(year).add(1).format('%d');
  
  // === 逻辑 1: 固定时间 (WMO) ===
  // 适用于 HK, Cairo, Singapore, Sao Paulo
  // Spring: 03-01 ~ 05-31
  // Summer: 06-01 ~ 08-31
  // Autumn: 09-01 ~ 11-30
  // Winter: 12-01 ~ 下一年02-28
  if (cfg.method === 'fixed') {
    var nextFebLast = ee.Date.fromYMD(ee.Number(year).add(1), 3, 1).advance(-1, 'day');
    return ee.Feature(null, {
      'City': cfg.name, 'Year': year,
      'Spring_Start': yStr.cat('-03-01'), 'Spring_End': yStr.cat('-05-31'),
      'Summer_Start': yStr.cat('-06-01'), 'Summer_End': yStr.cat('-08-31'),
      'Autumn_Start': yStr.cat('-09-01'), 'Autumn_End': yStr.cat('-11-30'),
      'Winter_Start': yStr.cat('-12-01'), 'Winter_End': nextFebLast.format('YYYY-MM-dd')
    });
  }

  // === 逻辑 2: 自适应 (JJJ, Berlin, etc) ===
  // 获取数据范围：当年1月 - 下一年4月 (为了计算 Winter End)
  var startDate = yStr.cat('-01-01');
  var endDate   = nextYStr.cat('-04-30');
  
  var dailyTemps = era5.filterDate(startDate, endDate)
    .map(function(img) {
      var t = img.reduceRegion({
        reducer: ee.Reducer.mean(), geometry: cfg.geom, scale: 11132
      }).get('temperature_2m');
      var tC = ee.Algorithms.If(t, ee.Number(t).subtract(273.15), -999);
      return ee.Feature(null, {'temp': tC, 'millis': img.date().millis(), 'date': img.date().format('YYYY-MM-dd')});
    });
  
  var fcSmooth = applyMovingAverage(ee.FeatureCollection(dailyTemps).sort('millis'));
  
  // 1. Spring Start: 2月15日后第一个 > 阈值
  var t_search_spring = ee.Date.fromYMD(year, 2, 15).millis();
  var springObj = fcSmooth.filter(ee.Filter.and(ee.Filter.gt('millis', t_search_spring), ee.Filter.gte('temp_5d', cfg.t_spring))).sort('millis').first();
  var springStart = ee.Date(ee.Algorithms.If(springObj, springObj.get('millis'), ee.Date.fromYMD(year, 3, 15).millis()));
  
  // 2. Summer Start: Spring后第一个 > 阈值
  var summerObj = fcSmooth.filter(ee.Filter.and(ee.Filter.gt('millis', springStart.millis()), ee.Filter.gte('temp_5d', cfg.t_summer))).sort('millis').first();
  var defaultSummer = (cfg.name === 'Murmansk') ? ee.Date.fromYMD(year, 6, 15) : ee.Date.fromYMD(year, 6, 1);
  var summerStart = ee.Date(ee.Algorithms.If(summerObj, summerObj.get('millis'), defaultSummer.millis()));
  
  // 3. Autumn Start: Summer后第一个 < 阈值 (入秋)
  // 搜索范围：7月15日后
  var t_search_autumn = ee.Date.fromYMD(year, 7, 15).millis();
  var autumnObj = fcSmooth.filter(ee.Filter.and(ee.Filter.gt('millis', t_search_autumn), ee.Filter.lt('temp_5d', cfg.t_summer))).sort('millis').first();
  var autumnStart = ee.Date(ee.Algorithms.If(autumnObj, autumnObj.get('millis'), ee.Date.fromYMD(year, 9, 1).millis()));
  
  // 4. Winter Start: Autumn后第一个 < 阈值 (入冬)
  var winterStartObj = fcSmooth.filter(ee.Filter.and(ee.Filter.gt('millis', autumnStart.millis()), ee.Filter.lt('temp_5d', cfg.t_winter))).sort('millis').first();
  var winterStart = ee.Date(ee.Algorithms.If(winterStartObj, winterStartObj.get('millis'), ee.Date.fromYMD(year, 11, 1).millis()));
  
  // 5. Winter End (Next Year Spring Start - 1 day)
  // 我们需要找下一年的春天开始时间
  var nextSpringSearch = ee.Date.fromYMD(ee.Number(year).add(1), 2, 15).millis();
  var nextSpringObj = fcSmooth.filter(ee.Filter.and(ee.Filter.gt('millis', nextSpringSearch), ee.Filter.gte('temp_5d', cfg.t_spring))).sort('millis').first();
  var nextSpringStart = ee.Date(ee.Algorithms.If(nextSpringObj, nextSpringObj.get('millis'), ee.Date.fromYMD(ee.Number(year).add(1), 3, 15).millis()));
  var winterEnd = nextSpringStart.advance(-1, 'day');

  // 计算中间的 End 日期
  var springEnd = summerStart.advance(-1, 'day');
  var summerEnd = autumnStart.advance(-1, 'day');
  var autumnEnd = winterStart.advance(-1, 'day');

  return ee.Feature(null, {
    'City': cfg.name, 'Year': year,
    'Spring_Start': springStart.format('YYYY-MM-dd'), 'Spring_End': springEnd.format('YYYY-MM-dd'),
    'Summer_Start': summerStart.format('YYYY-MM-dd'), 'Summer_End': summerEnd.format('YYYY-MM-dd'),
    'Autumn_Start': autumnStart.format('YYYY-MM-dd'), 'Autumn_End': autumnEnd.format('YYYY-MM-dd'),
    'Winter_Start': winterStart.format('YYYY-MM-dd'), 'Winter_End': winterEnd.format('YYYY-MM-dd')
  });
};

var resultList = ee.List([]);
cityConfig.forEach(function(cfg) {
  var cityRes = years.map(function(y) { return getDatesForCode(cfg, y); });
  resultList = resultList.cat(cityRes);
});

Export.table.toDrive({
  collection: ee.FeatureCollection(resultList),
  description: 'LCZ_Dates_For_V304_Code',
  fileFormat: 'CSV'
});