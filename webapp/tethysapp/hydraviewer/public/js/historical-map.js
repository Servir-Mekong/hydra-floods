var map;
var water_layer;
var water_source;
var climo;

$(function() {

  climo = 1
  slideMons = {1:'Jan.',2:'Feb.',3:'Mar.',4:'Apr.',5:'May',6:'Jun.'
         ,7:'Jul.',8:'Aug.',9:'Sep.',10:'Oct.',11:'Nov.',12:'Dec.'}
  slideFormat = {
  	formatter: function(value) {
  		return slideMons[value];
  	}
  }


  // With JQuery
  var slider = $('#slider1').slider(slideFormat);

  // Get the Open Layers map object from the Tethys MapView
  var map = TETHYS_MAP_VIEW.getMap();
  var $layers_element = $('#layers');

  var base_map = new ol.layer.Tile({
            crossOrigin: 'anonymous',
            source: new ol.source.XYZ({
                // attributions: [attribution],
                url: 'https://services.arcgisonline.com/ArcGIS/rest/services/Canvas/' +
                'World_Light_Gray_Base/MapServer/tile/{z}/{y}/{x}'
            })
        });

  map.getLayers().item(0).setVisible(false);
  map.addLayer(base_map);

  water_url = $layers_element.attr('data-water-url');

  water_source = new ol.source.XYZ({url:water_url});
  water_layer = new ol.layer.Tile(
    {
      source:water_source
    }
  );

  map.addLayer(water_layer)

  $("#climo-enabled").click(function() {
  	if(this.checked) {
  		// With JQuery
  		slider.slider("enable");
      document.getElementById("climo-selector").innerHTML = "&emsp;Yes";
      climo = 1
  	}
  	else {
  		// With JQuery
  		slider.slider("disable");
      document.getElementById("climo-selector").innerHTML = "&emsp;No";
      climo = 0
	  }
  });

  $('[name="update-button"]').on('click',function() {

    $('#spinner').show();
    var start = $('#date_picker1').val()
    var end = $('#date_picker2').val()
    var mon = $('input[id="slider1"]').slider('getValue')
    var algo = $('#algorithm_selection').val()

    console.log(start,end,mon,algo,climo)
    var xhr = ajax_update_database('update_historical',{'sDate':start,'eDate':end,'month':mon,'algo':algo,'climo':climo});
    xhr.done(function(data) {
        if("success" in data) {
          console.log(data)
          water_source = new ol.source.XYZ({url:data.url});
          water_layer.setSource(water_source)
          $('#spinner').hide();
        }else{
            alert('Opps, there was a problem processing the request. Please see the following error: '+data.error);
            $('#spinner').hide();
        }
    });
  });


  return
});
