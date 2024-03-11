var iniDate = '2023-08-01';
var endDate = '2023-11-30';

// cloud %
var cloudPerc = 15

Map.centerObject(geometry, 7)

// Import Collections //

// sentinel-2 L1C
var MSI = ee.ImageCollection('COPERNICUS/S2');
// sentinel-2 L2A atmospherically corrected
var MSIA = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")

// landsat-8 surface reflactance product (for masking purposes)
var SRP = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');

// toms / omi
var ozone = ee.ImageCollection('TOMS/MERGED');

var pi = ee.Image(3.141592);

// water mask
var startMonth = 5;
var endMonth = 9;
var startYear = 2020;
var endYear = 2023;

var forMask = SRP.filterBounds(geometry).select('B6').filterMetadata('CLOUD_COVER', "less_than", cloudPerc).filter(ee.Filter.calendarRange(startMonth, endMonth, 'month')).filter(ee.Filter.calendarRange(startYear, endYear, 'year'));
var mask = ee.Image(forMask.select('B6').median().lt(300)) 
mask = mask.updateMask(mask)

// filter sentinel 2 collection
var FC = MSI.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than", cloudPerc);
print(FC, 'Sentinel 2 Collection')
var FCA = MSIA.filterDate(iniDate, endDate).filterBounds(geometry).filterMetadata('CLOUDY_PIXEL_PERCENTAGE', "less_than", cloudPerc);

// single image
var FCList=FC.toList(999);
var s2=ee.Image(ee.List(FCList).get(0)); //note index 0 is the first image
print('Single Image',s2)


// Collection Processing //

// atmospheric correction
var Rrs_coll = FC.map(s2Correction);

// chlorophyll-a
var chlor_a_coll = Rrs_coll.map(chlorophyll);

// turbidity
var NDTI_coll = Rrs_coll.map(turbidity);
var Turbidityn = FCA.map(turbidityntu);

// sd
var sd_coll = Rrs_coll.map(secchi);


// Mapping functions //

function s2Correction(img){

// msi bands 
var bands = ['B1','B2','B3','B4','B5','B6','B7', 'B8', 'B8A', 'B11', 'B12'];

// rescale
var rescale = img.select(bands).divide(10000).multiply(mask)

// tile footprint
var footprint = rescale.geometry()

// dem
var DEM = ee.Image('USGS/SRTMGL1_003').clip(footprint); 

// ozone
var DU = ee.Image(ozone.filterDate(iniDate,endDate).filterBounds(footprint).mean());

//Julian Day
var imgDate = ee.Date(img.get('system:time_start'));
var FOY = ee.Date.fromYMD(imgDate.get('year'),1,1);
var JD = imgDate.difference(FOY,'day').int().add(1);

// earth-sun distance
var myCos = ((ee.Image(0.0172).multiply(ee.Image(JD).subtract(ee.Image(2)))).cos()).pow(2)
var cosd = myCos.multiply(pi.divide(ee.Image(180))).cos();
var d = ee.Image(1).subtract(ee.Image(0.01673)).multiply(cosd).clip(footprint)

// sun azimuth
var SunAz = ee.Image.constant(img.get('MEAN_SOLAR_AZIMUTH_ANGLE')).clip(footprint);

// sun zenith
var SunZe = ee.Image.constant(img.get('MEAN_SOLAR_ZENITH_ANGLE')).clip(footprint);
var cosdSunZe = SunZe.multiply(pi.divide(ee.Image(180))).cos(); // in degrees
var sindSunZe = SunZe.multiply(pi.divide(ee.Image(180))).sin(); // in degrees

// sat zenith
var SatZe = ee.Image.constant(img.get('MEAN_INCIDENCE_ZENITH_ANGLE_B5')).clip(footprint);
var cosdSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).cos();
var sindSatZe = (SatZe).multiply(pi.divide(ee.Image(180))).sin();
  
// sat azimuth
var SatAz = ee.Image.constant(img.get('MEAN_INCIDENCE_AZIMUTH_ANGLE_B5')).clip(footprint);

// relative azimuth
var RelAz = SatAz.subtract(SunAz);
var cosdRelAz = RelAz.multiply(pi.divide(ee.Image(180))).cos();

// Pressure
var P = (ee.Image(101325).multiply(ee.Image(1).subtract(ee.Image(0.0000225577).multiply(DEM)).pow(5.25588)).multiply(0.01)).multiply(mask);
var Po = ee.Image(1013.25);

// esun
var ESUN = ee.Image(ee.Array([ee.Image(img.get('SOLAR_IRRADIANCE_B1')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B2')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B3')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B4')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B5')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B6')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B7')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B8')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B8A')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B11')),
                  ee.Image(img.get('SOLAR_IRRADIANCE_B2'))]
                  )).toArray().toArray(1);

ESUN = ESUN.multiply(ee.Image(1))

var ESUNImg = ESUN.arrayProject([0]).arrayFlatten([bands]);

// create empty array for the images
var imgArr = rescale.select(bands).toArray().toArray(1);

// pTOA to Ltoa
var Ltoa = imgArr.multiply(ESUN).multiply(cosdSunZe).divide(pi.multiply(d.pow(2)));

// band centers
var bandCenter = ee.Image(443).divide(1000).addBands(ee.Image(490).divide(1000))
                                        .addBands(ee.Image(560).divide(1000))
                                        .addBands(ee.Image(665).divide(1000))
                                        .addBands(ee.Image(705).divide(1000))
                                        .addBands(ee.Image(740).divide(1000))
                                        .addBands(ee.Image(842).divide(1000))
                                        .addBands(ee.Image(865).divide(1000))
                                        .addBands(ee.Image(1610).divide(1000))
                                        .addBands(ee.Image(2190).divide(1000))
                                        .toArray().toArray(1);

// ozone coefficients
var koz = ee.Image(0.0039).addBands(ee.Image(0.0213))
                        .addBands(ee.Image(0.1052))
                        .addBands(ee.Image(0.0505))
                        .addBands(ee.Image(0.0205))
                        .addBands(ee.Image(0.0112))
                        .addBands(ee.Image(0.0075))
                        .addBands(ee.Image(0.0021))
                        .addBands(ee.Image(0.0019))                          
                        .addBands(ee.Image(0))
                        .addBands(ee.Image(0))
                        .toArray().toArray(1);

// Calculate ozone optical thickness
var Toz = koz.multiply(DU).divide(ee.Image(1000));

// Calculate TOA radiance in the absense of ozone
var Lt = Ltoa.multiply(((Toz)).multiply((ee.Image(1).divide(cosdSunZe)).add(ee.Image(1).divide(cosdSatZe))).exp());

// Rayleigh optical thickness
var Tr = (P.divide(Po)).multiply(ee.Image(0.008569).multiply(bandCenter.pow(-4))).multiply((ee.Image(1).add(ee.Image(0.0113).multiply(bandCenter.pow(-2))).add(ee.Image(0.00013).multiply(bandCenter.pow(-4)))));

// Specular reflection (s- and p- polarization states)
var theta_V = ee.Image(0.0000000001);
var sin_theta_j = sindSunZe.divide(ee.Image(1.333));

var theta_j = sin_theta_j.asin().multiply(ee.Image(180).divide(pi));

var theta_SZ = SunZe;

var R_theta_SZ_s = (((theta_SZ.multiply(pi.divide(ee.Image(180)))).subtract(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)).divide((((theta_SZ.multiply(pi.divide(ee.Image(180)))).add(theta_j.multiply(pi.divide(ee.Image(180))))).sin().pow(2)));

var R_theta_SZ_p = (((theta_SZ.multiply(pi.divide(180))).subtract(theta_j.multiply(pi.divide(180)))).tan().pow(2)).divide((((theta_SZ.multiply(pi.divide(180))).add(theta_j.multiply(pi.divide(180)))).tan().pow(2)));

var R_theta_V_p = ee.Image(0.0000000001);

var R_theta_SZ = ee.Image(0.5).multiply(R_theta_SZ_s.add(R_theta_SZ_p));

var R_theta_V = ee.Image(0.5).multiply(R_theta_V_s.add(R_theta_V_p));
  
// Sun-sensor geometry
var theta_neg = ((cosdSunZe.multiply(ee.Image(-1))).multiply(cosdSatZe)).subtract((sindSunZe).multiply(sindSatZe).multiply(cosdRelAz));

var theta_neg_inv = theta_neg.acos().multiply(ee.Image(180).divide(pi));

var theta_pos = (cosdSunZe.multiply(cosdSatZe)).subtract(sindSunZe.multiply(sindSatZe).multiply(cosdRelAz));

var theta_pos_inv = theta_pos.acos().multiply(ee.Image(180).divide(pi));

var cosd_tni = theta_neg_inv.multiply(pi.divide(180)).cos(); // in degrees

var cosd_tpi = theta_pos_inv.multiply(pi.divide(180)).cos(); // in degrees

var Pr_neg = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tni.pow(2))));

var Pr_pos = ee.Image(0.75).multiply((ee.Image(1).add(cosd_tpi.pow(2))));
  
// Rayleigh scattering phase function
var Pr = Pr_neg.add((R_theta_SZ.add(R_theta_V)).multiply(Pr_pos));

// rayleigh radiance contribution
var denom = ee.Image(4).multiply(pi).multiply(cosdSatZe);
var Lr = (ESUN.multiply(Tr)).multiply(Pr.divide(denom));

// rayleigh corrected radiance
var Lrc = Lt.subtract(Lr);
var LrcImg = Lrc.arrayProject([0]).arrayFlatten([bands]);

// Aerosol Correction //

// Bands in nm
var bands_nm = ee.Image(443).addBands(ee.Image(490))
                            .addBands(ee.Image(560))
                            .addBands(ee.Image(665))
                            .addBands(ee.Image(705))
                            .addBands(ee.Image(783))
                            .addBands(ee.Image(842))
                            .addBands(ee.Image(865))                            
                            .addBands(ee.Image(0))
                            .addBands(ee.Image(0))
                            .toArray().toArray(1);

// Lam in SWIR bands
var Lam_10 = LrcImg.select('B11');
var Lam_11 = LrcImg.select('B12');

// Calculate aerosol type
var eps = ((((Lam_11).divide(ESUNImg.select('B12'))).log()).subtract(((Lam_10).divide(ESUNImg.select('B11'))).log())).divide(ee.Image(2190).subtract(ee.Image(1610)));

// Calculate multiple scattering of aerosols for each band
var Lam = (Lam_11).multiply(((ESUN).divide(ESUNImg.select('B12')))).multiply((eps.multiply(ee.Image(-1))).multiply((bands_nm.divide(ee.Image(2190)))).exp());

// diffuse transmittance
var trans = Tr.multiply(ee.Image(-1)).divide(ee.Image(2)).multiply(ee.Image(1).divide(cosdSatZe)).exp();

// Compute water-leaving radiance
var Lw = Lrc.subtract(Lam).divide(trans);

// water-leaving reflectance
var pw = (Lw.multiply(pi).multiply(d.pow(2)).divide(ESUN.multiply(cosdSunZe)));

// remote sensing reflectance
var Rrs_coll = (pw.divide(pi).arrayProject([0]).arrayFlatten([bands]).slice(0,9));

return(Rrs_coll.set('system:time_start',img.get('system:time_start')));

}
function chlorophyll(img){
  var NDCI_coll = (img.select('B5').subtract(img.select('B4'))).divide(img.select('B5').add(img.select('B4')));
  var chlor_a_coll = ee.Image(14.039).add(ee.Image(86.115).multiply(NDCI_coll)).add(ee.Image(194.325).multiply(NDCI_coll.pow(ee.Image(2))));
  return(chlor_a_coll.updateMask(chlor_a_coll.lt(100)).set('system:time_start',img.get('system:time_start')))
}
function turbidity(img){
  var NDTI_coll = img.normalizedDifference(['B4','B3']);
  return(NDTI_coll.updateMask(NDTI_coll.lt(100)).set('system:time_start',img.get('system:time_start')))
}
function secchi(img){
  var blueRed_coll = (img.select('B2').divide(img.select('B4'))).log()
  var lnMOSD_coll = (ee.Image(1.4856).multiply(blueRed_coll)).add(ee.Image(0.2734)); // R2 = 0.8748 with Anthony's in-situ data
  var MOSD_coll = ee.Image(10).pow(lnMOSD_coll);
  var sd_coll = (ee.Image(1.1777).multiply(MOSD_coll)).add(ee.Image(1.0813));
  return(sd_coll.updateMask(sd_coll.lt(10)).set('system:time_start',img.get('system:time_start')))
}
function trophicState(img){
  var tsi_coll =  ee.Image(30.6).add(ee.Image(9.81).multiply(img.log()));
  return(tsi_coll.updateMask(tsi_coll.lt(200)).set('system:time_start',img.get('system:time_start')))
}
function turbidityntu(img){
  var turb_ntu =  ee.Image(294.79).multiply(img.select('B5').multiply(img.select('B5').divide(img.select('B2')))).add(ee.Image(0.9061));
  return(turb_ntu.updateMask(mask).divide(ee.Image(10000)).set('system:time_start',img.get('system:time_start')))
}
// Zhan, Y., Delegido, J., Erena, M., Soria, J.M., Ruiz-Verdú, A., Urrego, P., Sòria-Perpinyà, X., Vicente, E., Moreno, J., 2022. Mar Menor lagoon (SE Spain) chlorophyll-a and turbidity estimation with Sentinel-2. Limnetica 41, 1. https://doi.org/10.23818/limn.41.18


// Time Series //

// Chlorophyll-a time series
var chlorTimeSeries = ui.Chart.image.seriesByRegion(
    chlor_a_coll, geometry2, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean Chlorphyll-a',
          vAxis: {title: 'Chlor-a [micrograms/L]'},
          lineWidth: 1,
          pointSize: 4,
          legend: 'Chlorphyll-a',
});

// turbidity time series
var turbTimeSeries = ui.Chart.image.seriesByRegion(
    Turbidityn, geometry2, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean turbidity',
          vAxis: {title: 'turbidity [NTU/DTU]'},
          lineWidth: 1,
          pointSize: 4,
          legend: 'turbidity',
});

// SD time series
var sdTimeSeries = ui.Chart.image.seriesByRegion(
    sd_coll, geometry2, ee.Reducer.mean())
        .setChartType('ScatterChart')
        .setOptions({
          title: 'Mean Secchi Depth',
          vAxis: {title: 'Zsd [m]'},
          lineWidth: 1,
          pointSize: 4,
});


print(chlorTimeSeries)
print(turbTimeSeries)
print(sdTimeSeries)

// Map Layers //
Map.addLayer(mask, {}, 'mask', false)
Map.addLayer(chlor_a_coll.mean(), {min: 0, max: 40, palette: ['darkblue','blue','cyan','limegreen','yellow', 'orange', 'orangered', 'darkred']}, 'Mean chlor-a', false);
Map.addLayer(NDTI_coll.mean(), {min: -1, max: 1, palette: ['darkblue','blue','cyan','limegreen','yellow', 'orange', 'orangered', 'darkred']}, 'Mean turbidity', false);
Map.addLayer(sd_coll.mean(), {min: 0, max: 2, palette: ['800000', 'FF9700', '7BFF7B', '0080FF', '000080']}, 'Mean Zsd', false);
Map.addLayer(FC.mean(),{bands:['B4','B3','B2'],min:0,max:5000},'RGB', false)
Map.addLayer(Turbidityn.mean(), {min: 0, max: 100, palette: ['darkblue','blue','cyan','limegreen','yellow', 'orange', 'orangered', 'darkred']}, 'turbidity ntu', false);





var FCAT = FCA.map(turbidityntuu);
function turbidityntuu(img){
  var turb_ntu =  ee.Image(194.79).multiply(img.select('B5').multiply(img.select('B5').divide(img.select('B2')))).add(ee.Image(0.9061));
  return(img.addBands(turb_ntu.updateMask(mask))  ) // .divide(ee.Image(10000))
}

var S2_comp =  FCAT.map(function(image){
  var start = ee.Date(image.get('system:time_start'));
  var end = ee.Date(image.get('system:time_end'));
  var label = start.format('YYYY-MM-dd');//.cat(' - ').cat(end.format('YYYY-MM-dd'));
  
  return image.visualize({
    forceRgbOutput: true,
    //palette: palettes.BrBG[9], //original palette was "000000", "fdbb84"
    bands: ['constant', 'B3', 'B2'],
    min: 440,
    max: 2190
  }).clip(geometry).set({label: label}); 
});

var videoArgs = {
 dimensions: 600,
 region: geometry,
 framesPerSecond: 5,
 format:'gif',
 maxPixels: 26214400,
  //crs: 'EPSG:3857',
};
var text = require('users/gena/packages:text');

var annotations = [
  {
    position: 'right', offset: '5%', margin: '20%', property: 'label', scale: 150
  }
]
  
S2_comp = S2_comp.map(function(image) {
  return text.annotateImage(image, {}, geometry, annotations) // bounds, 
});

print(S2_comp.getVideoThumbURL(videoArgs));

