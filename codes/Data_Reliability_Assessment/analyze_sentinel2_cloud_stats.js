// ==============================================================
// 独立脚本：专门用于统计 Sentinel-2 历史云量 (修正版 - L1C)
// ==============================================================

var ROI_ASSET = "projects/ee-l2892786691/assets/JJJ/jjj"; // 你的 ROI
var jjj_vec = ee.FeatureCollection(ROI_ASSET);
var ROI = jjj_vec.geometry();
Map.centerObject(ROI, 7);

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

// 【关键修改】使用 L1C 集合以获取 2017-2018 的完整历史数据
// SR (L2A) 在 2019 年之前在中国区域严重缺失，会导致云量统计错误 (接近0)
var S2_HISTORY = ee.ImageCollection('COPERNICUS/S2_HARMONIZED'); 

var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];
var SEASONS_ORDER = ['SPRING', 'SUMMER', 'AUTUMN', 'WINTER'];

var statsList = [];

print('开始统计... (请稍候 10-20 秒)');

years.forEach(function(year) {
  SEASONS_ORDER.forEach(function(season) {
    var dates = YEARLY_DATES[year];
    var range = ee.Dictionary(ee.Dictionary(dates).get(season));
    
    // 关键修正：确保没有任何 map(mask) 操作
    // 只过滤掉覆盖率 > 80% 的太差影像，保留其他的用于统计平均水平
    var s2Col = S2_HISTORY.filterBounds(ROI)
                          .filterDate(range.get('start'), range.get('end'))
                          .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 80)); 

    // 计算平均云量 (属性级聚合 - 速度快且包含被Mask掉的云)
    var mean_cloud = s2Col.aggregate_mean('CLOUDY_PIXEL_PERCENTAGE');
    
    // 统计数量，验证数据完整性
    var count = s2Col.size();
    
    // 处理空值
    var val = ee.Algorithms.If(count.gt(0), mean_cloud, -1);
    
    statsList.push(ee.Feature(null, {
      'Year': year,
      'Season': season,
      'Cloud_Rate': val, // 这次出来的数应该是 30-60 左右
      'Count': count     // 2017年也应该有几百上千张
    }));
  });
});

var statsFC = ee.FeatureCollection(statsList);

// 1. 直接打印数值检查
print('统计结果列表 (Check Values):', statsFC);

// 2. 导出 CSV
Export.table.toDrive({
  collection: statsFC,
  description: 'JJJ_Cloud_Stats_L1C_Corrected',
  folder: 'LCZ_Stats',
  fileFormat: 'CSV',
  selectors: ['Year', 'Season', 'Cloud_Rate', 'Count']
});

// 3. 画个图看看趋势
var chartCloud = ui.Chart.feature.byFeature(statsFC, 'Year', ['Cloud_Rate'])
    .setChartType('ColumnChart')
    .setOptions({
        title: 'Sentinel-2 Historical Cloud Rate (L1C Basis)',
        vAxis: {title: 'Avg Cloud Cover (%)', viewWindow: {min: 0, max: 80}},
        colors: ['#d73027'] // 红色
    });
print(chartCloud);

var chartCount = ui.Chart.feature.byFeature(statsFC, 'Year', ['Count'])
    .setChartType('LineChart')
    .setOptions({
        title: 'Image Availability Count (L1C Basis)',
        vAxis: {title: 'Number of Images'},
        colors: ['#4575b4'] // 蓝色
    });
print(chartCount);