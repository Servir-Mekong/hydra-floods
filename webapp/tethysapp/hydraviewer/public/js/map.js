/* global variables to be tossed around like hot potatoes */
var map,
    selected_date,
    browse_layer,
    precip_layer,
    historical_layer,
    sentinel1_layer,
    admin_layer,
    flood_layer,
    drawing_polygon,
    $layers_element;

// init function for the page
$(function() {

  selected_date = $('#date_selection').val();

  $('.js-range-slider').ionRangeSlider({
        skin: "round",
	type: "double",
        grid: true,
	from: 0,
	to: 11,
        values: [
            "Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
        ]
    });


  // set slider vars with the correct elements
  var browseSlider = $('#browse-opacity').slider();
  var precipSlider = $('#precip-opacity').slider();
  var historicalSlider = $('#historical-opacity').slider();
  var sentinel1Slider = $('#sentinel1-opacity').slider();
  var floodSlider = $('#flood-opacity').slider();
  var floodSlider1 = $('#flood1-opacity').slider();


  // check if this init is a load from the use case redirects
  // only need to change map center for one use case (only works on production server)
  if (window.location.href == 'http://tethys-servir.adpc.net/apps/hydraviewer/mapviewer/?sDate=2016-07-14&sensor_txt=sentinel1')  {
    var centerPt = [21,94]
  }
  else{
    centerPt = [16.8,95.7]
  }
  console.log(centerPt)

  // init map
  map = L.map('map',{
    center: centerPt,
    zoom: 8,
    minZoom:2,
    maxZoom: 16,
    maxBounds: [
     [-120, -220],
     [120, 220]
   ],
 });

// Initialise the FeatureGroup to store editable layers
var editableLayers = new L.FeatureGroup();
map.addLayer(editableLayers);

var drawPluginOptions = {
  draw: {
    polygon: {
      allowIntersection: false, // Restricts shapes to simple polygons
      drawError: {
        color: '#e1e100', // Color the shape will turn when intersects
        message: '<strong>Oh snap!<strong> you can\'t draw that!' // Message that will show when intersect
      },
      shapeOptions: {
        color: '#97009c'
      }
    },
    // disable toolbar item by setting it to false
    polyline: false,
    circle: false, // Turns off this drawing tool
    circlemarker: false,
    rectangle: true,
    marker: false,
    },
  edit: {
    featureGroup: editableLayers, //REQUIRED!!
    remove: false
  }
};

// Initialise the draw control and pass it the FeatureGroup of editable layers
var drawControl = new L.Control.Draw(drawPluginOptions);
map.addControl(drawControl);

var editableLayers = new L.FeatureGroup();
map.addLayer(editableLayers);

map.on('draw:created', function(e) {
  editableLayers.clearLayers();
  var type = e.layerType,
    layer = e.layer;

  if (type === 'marker') {
    layer.bindPopup('A popup!');
  }
  userPolygon = layer.toGeoJSON();
  drawing_polygon = JSON.stringify(userPolygon.geometry.coordinates[0]);


  editableLayers.addLayer(layer);
});

  var positron = L.tileLayer('http://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png', {
        attribution: '©OpenStreetMap, ©CartoDB'
      }).addTo(map);

  $layers_element = $('#layers');
  var $update_element = $('#update_button');
  var viirs_product = "VIIRS_SNPP_CorrectedReflectance_BandsM11-I2-I1"

  //browse_layer = addGibsLayer(browse_layer,viirs_product,selected_date)
  browse_layer = addGibsLayer(browse_layer,viirs_product,'2019-05-24')


  //sentinel1_layer = addMapLayer(sentinel1_layer,$layers_element.attr('data-sentinel1-url'))
  flood_layer = addMapLayer(flood_layer,$layers_element.attr('data-flood-url'))
  historical_layer = addMapLayer(historical_layer,$layers_element.attr('data-historical-url'))
  precip_layer = addMapLayer(precip_layer,$layers_element.attr('data-precip-url'))
  admin_layer = addMapLayer(admin_layer,$layers_element.attr('data-admin-url'))

  precip_layer.setOpacity(0)
  precipSlider.slider('disable')

  $('#date_selection').change(function(){
    selected_date = $('#date_selection').val();
    var prec = $('#product_selection').val();
    var cmap = $('#cmap_selection').val()
    var accum = prec.split('|')[0]
    var xhr = ajax_update_database('get_precipmap',{'sDate':selected_date,'accum':accum,'cmap':cmap},"layers");
    xhr.done(function(data) {
        if("success" in data) {
          precip_layer.setUrl(data.url)
        }else{
          precip_layer.setUrl('')
          alert('Opps, there was a problem processing the request. Please see the following error: '+data.error);
        }
    });

    var prod = $('#browse_selection').val();
    var id = prod.split('|')[1]
    var template =
      '//gibs-{s}.earthdata.nasa.gov/wmts/epsg3857/best/' +
      id + '/default/' + selected_date + '/{tileMatrixSet}/{z}/{y}/{x}.jpg'
    browse_layer.setUrl(template)

    var sensor_val = $('#sensor_selection').val();
    console.log(sensor_val)

    if (sensor_val != 'none'){
      var xhr = ajax_update_database('get_surfacewatermap',{'sDate':selected_date,'sensor_txt':sensor_val},"layers");
      xhr.done(function(data) {
          if("success" in data) {
            flood_layer.setUrl(data.url)
          }else{
            flood_layer.setUrl('')
            alert(data.error);
          }
      });
    }
  });

  $('#cmap_selection').change(function(){
    var prod = $('#product_selection').val();
    var cmap = $('#cmap_selection').val()
    var accum = prod.split('|')[0]
    var xhr = ajax_update_database('get_precipmap',{'sDate':selected_date,'accum':accum,'cmap':cmap},"layers");
    xhr.done(function(data) {
        if("success" in data) {
          precip_layer.setUrl(data.url)
          $("#precip-cb").attr("src", cb_url+"?rnd="+Math.random())
        }else{
          alert('Opps, there was a problem processing the request. Please see the following error: '+data.error);
        }
    });
  });

  $('#product_selection').change(function(){
    var prod = $('#product_selection').val();
    var cmap = $('#cmap_selection').val()
    var accum = prod.split('|')[0]
    var xhr = ajax_update_database('get_precipmap',{'sDate':selected_date,'accum':accum,'cmap':cmap},"layers");
    xhr.done(function(data) {
        if("success" in data) {
          precip_layer.setUrl(data.url)
          $("#precip-cb").attr("src", cb_url+"?rnd="+Math.random())
        }else{
          alert('Opps, there was a problem processing the request. Please see the following error: '+data.error);
        }
    });
  });

  $('#sensor_selection').change(function(){
    var sensor_val = $('#sensor_selection').val();
    var xhr = ajax_update_database('get_surfacewatermap',{'sDate':selected_date,'sensor_txt':sensor_val},"layers");
    xhr.done(function(data) {
        if("success" in data) {
          flood_layer.setUrl(data.url)
        }else{
          alert(data.error);
        }
    });

  });

  $("#sensor_selection option[value='atms']").attr('disabled','disabled');

  $("#update-button").on("click",function(){
    var startYear = $('#start_year_selection_historical').val();
    var endYear = $('#end_year_selection_historical').val();
    var slider = $("#month_range").data("ionRangeSlider");

   // Get values
    var startMonth = slider.result.from + 1;
    var endMonth= slider.result.to + 1;
    var method = 'discrete';

    if (startMonth == endMonth) { endMonth += 1 }

    var xhr = ajax_update_database('update_historical',{'startYear':startMonth,'endYear':endYear,'startMonth': startMonth,'endMonth': endMonth, 'method': method},"layers");
        xhr.done(function(data) {
        if("success" in data) {
          historical_layer.setUrl(data.url)
        }else{
          alert(data.error);
        }
    });
  });


  $("#btn_download").on("click",function(){
 if(drawing_polygon === undefined){
  alert("Please draw a polygon");
  }else{
     // var selected_date = $('#selected_date').val();
     var sensor_val = $('#sensor_selection').val();
     var xhr = ajax_update_database('download_surfacewatermap',{'sDate':selected_date,'sensor_txt':sensor_val,'poly_coordinates': drawing_polygon},"json");
    xhr.done(function(data) {
        if("success" in data) {
          //alert('Download URL: \n'+ data.url)
          window.open(data.url, '_blank');
        }else{
          alert('Opps, there was a problem processing the request. Please see the following error: '+data.error);
        }
    });
}
  });

  $('#browse_selection').change(function(){
    var prod = $('#browse_selection').val();
    var id = prod.split('|')[1]
    var template =
      '//gibs-{s}.earthdata.nasa.gov/wmts/epsg3857/best/' +
      id + '/default/' + selected_date + '/{tileMatrixSet}/{z}/{y}/{x}.jpg'
    browse_layer.setUrl(template)
  });

  $('#browse-opacity').change(function(){
    var opac = parseFloat($('input[id="browse-opacity"]').slider('getValue'))
    browse_layer.setOpacity(opac)
  });

  $('#precip-opacity').change(function(){
    var opac = parseFloat($('input[id="precip-opacity"]').slider('getValue'))
    precip_layer.setOpacity(opac)
  });

  $('#historical-opacity').change(function(){
    var opac = parseFloat($('input[id="historical-opacity"]').slider('getValue'))
    historical_layer.setOpacity(opac)
  });

  $('#flood1-opacity').change(function(){
    var opac = parseFloat($('input[id="flood1-opacity"]').slider('getValue'))
    flood_layer.setOpacity(opac)
  });

  $("#browse-check").on("click",function(){
    if(this.checked){
      browseSlider.slider('enable')
      var opac = parseFloat($('input[id="browse-opacity"]').slider('getValue'))
      browse_layer.setOpacity(opac)
    }
    else{
      browseSlider.slider('disable')
      browse_layer.setOpacity(0)
    }
  });

  $("#precip-check").on("click",function(){
    if(this.checked){
      precipSlider.slider('enable')
      var opac = parseFloat($('input[id="precip-opacity"]').slider('getValue'))
      precip_layer.setOpacity(opac)
    }
    else{
      precipSlider.slider('disable')
      precip_layer.setOpacity(0)
    }
  });

  $("#historical-check").on("click",function(){
    if(this.checked){
      historicalSlider.slider('enable')
      var opac = parseFloat($('input[id="historical-opacity"]').slider('getValue'))
      historical_layer.setOpacity(opac)
    }
    else{
      historicalSlider.slider('disable')
      historical_layer.setOpacity(0)
    }
  });

$("#flood-check").on("click",function(){
    if(this.checked){
      floodSlider1.slider('enable')
      var opac = parseFloat($('input[id="flood1-opacity"]').slider('getValue'))
      flood_layer.setOpacity(opac)
    }
    else{
      floodSlider1.slider('disable')
      flood_layer.setOpacity(0)
    }
  });

$("#download_flood-check").on("click",function(){
    if(this.checked ){
      $("#btn_download").removeAttr('disabled');
    }
    else{
      $("#btn_download").attr('disabled','disabled');
    }
  });


  $(".legend-info-button").click(function () {
    $(".legend-tabs").toggle();
    $("#legend-content").toggle();
    if ($("#legend-content").is(":visible") == true) {
      $("#legend-collapse").css("display","inline-block");
      $("#legend-expand").css("display","none");
    }
    else {
      $("#legend-collapse").css("display","none");
      $("#legend-expand").css("display","inline-block");
    }
  });

// end of init function
});

function openLegendTab(event, name) {
  tabcontent = document.getElementsByClassName("legend-tab-content");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  }
  tablinks = document.getElementsByClassName("legend-tab");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }
  document.getElementById('legend-tab-' + name).style.display = "block";
  event.currentTarget.className += " active";
}
$(document).ready(function() {
  document.getElementById("legend-tab-default").click();
});

// function to add and update tile layer to map
function addMapLayer(layer,url){
  layer = L.tileLayer(url,{attribution:
    '<a href="https://earthengine.google.com" target="_">' +
    'Google Earth Engine</a>;'}).addTo(map);
  return layer
}

function addGibsLayer(layer,product,date){
  var template =
    '//gibs-{s}.earthdata.nasa.gov/wmts/epsg3857/best/' +
    '{layer}/default/{time}/{tileMatrixSet}/{z}/{y}/{x}.jpg';

  layer = L.tileLayer(template, {
    layer: product,
    tileMatrixSet: 'GoogleMapsCompatible_Level9',
    maxZoom: 9,
    time: date,
    tileSize: 256,
    subdomains: 'abc',
    noWrap: true,
    continuousWorld: true,
    // Prevent Leaflet from retrieving non-existent tiles on the
    // borders.
    bounds: [
      [-85.0511287776, -179.999999975],
      [85.0511287776, 179.999999975]
    ],
    attribution:
      '<a href="https://wiki.earthdata.nasa.gov/display/GIBS" target="_">' +
      'NASA EOSDIS GIBS</a>;'
  });

  map.addLayer(layer);

  return layer
}

/*
* Workaround for 1px lines appearing in some browsers due to fractional transforms
* and resulting anti-aliasing.
* https://github.com/Leaflet/Leaflet/issues/3575
*/
// (function () {
//   var originalInitTile = L.GridLayer.prototype._initTile;
//   L.GridLayer.include({
//     _initTile: function (tile) {
//       originalInitTile.call(this, tile);
//
//       var tileSize = this.getTileSize();
//
//       tile.style.width = tileSize.x + 1 + 'px';
//       tile.style.height = tileSize.y + 1 + 'px';
//     }
//   });
// })();
