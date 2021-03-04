# Vorobox

## Importing GEOJSON 

_Vorobox_ can read input files which use a modified form of the _topojson_ format which captures topological relationships between data points. As there are many utilities available for converting _topojson_ files to and from _geojson_ it is possible to prepare a _geojson_ file so it can be used as input for _vorobox_.

Consider a the data describing national boundaries, coasylines and lakes from [Natural Earth](). Using Mike Bostock’s excellent [ndjson](https://github.com/mbostock/ndjson-cli) library this file can be filtered to extract the data describing the required countries. In this case let's select those of _Western Europe_ and _Northern Europe_. 
```bash
ndjson-split 'd.features' < ne_50m_admin_0_countries_lakes.json > natural_earth.ndjson
ndjson-filter 'd.properties.SUBREGION === "Western Europe" || d.properties.SUBREGION === "Northern Europe"' < natural_earth.ndjson > nw_europe.ndjson
```

The properties list associated with each feature now needs trimming to retain only a few components; namely "POP_EST", "WOE_ID", "NAME", "SOVEREIGNT", and here again _ndjson_ can help

```bash
ndjson-map 'd.properties = {nation: d.properties.SOVEREIGNT, id:d.properties.WOE_ID, numberPoints: Math.floor(0.001 * d.properties.POP_EST), name: d.properties.name}, d' < nw_europe.ndjson > nations.ndjson
```
Yielding data like this 
```json
{"type":"Feature","properties":{"nation":"Finland","id":12577865,"numberPoints":27},"geometry":{"type":"MultiPolygon","coordinates":[[[[20.611328125000057,60.04067382812502],[20.603417968750023,60.01694335937498],[20.521777343750017,60.01166992187498],[20.4875,60.03276367187502],[20.411230468750034,60.030126953125006],[20.39794921875,60.04067382812502],[20.42958984375005,60.061718749999955],[20.49013671875005,60.07490234374998],[20.569140625000074,60.069628906250045],[20.611328125000057,60.04067382812502]]],[[[19.662304687500068,60.18715820312502],[19.667480468750057,60.16474609374998],[19.629199218750074,60.17036132812498],[19.59980468750004,60.162695312500006],[19.579882812500045,60.135058593750074],[19.536523437500023,60.144970703124955],[19.51904296875,60.18457031250003],[19.551367187500006,60.24384765625001],[19.62880859375008,60.24609375000003],[19.662304687500068,60.18715820312502]]],[[[19.989550781250017,60.351171875000034],[20.020214843750068,60.35087890624999],[20.03388671875001,60.35932617187501],[20.08740234374997,60.353417968749966],[20.167871093750023,60.31469726562497],[20.184082031250057,60.29375],[20.239550781250074,60.28300781250002],[20.258886718750006,60.26127929687502],[20.19472656250008,60.19355468749998],[20.15507812499999,60.192285156249994],[20.12548828125,60.20087890624998],[20.07324218750003,60.19345703124998],[20.042578125,60.18066406250003],[20.032324218750034,60.152490234374966],[20.033984375000017,60.093554687500045],[19.7998046875,60.081738281249955],[19.74599609375005,60.098974609375006],[19.67226562500005,60.233007812500034],[19.686914062500023,60.26762695312499],[19.736523437499983,60.282373046874966],[19.77900390625004,60.28554687499999],[19.785253906250063,60.21337890625003],[19.84765624999997,60.22055664062498],[19.867187500000057,60.26811523437499],[19.871582031250057,60.30161132812498],[19.854687499999983,60.318505859374966],[19.812304687500074,60.331591796875045],[19.78779296875004,60.35405273437502],[19.823046875000045,60.39018554687499],[19.88828125,60.40581054687499],[19.94453125000001,60.357519531250006],[19.989550781250017,60.351171875000034]]]]}}
```

Now the data can be converted back to _geojson_, again using _ndjson_. 
```bash
ndjson-reduce 'p.features.push(d), p' '{type: "FeatureCollection", features: []}' < nations.ndjson > nations.geojson
```

But these data are in coordinates of latitude and longitude; it would be better to use an equal area projection 
```bash 
geoproject 'd3.geoAlbers()
    .rotate([-20.0, 0.0])
    .center([0.0, 52.0])
    .parallels([35.0, 65.0]).reflectY(true)' < nations.geojson > nw_albers.geojson
```
The resulting _geojson_ file can now be converted into topojson; which holds the topological information necessary for _vorobox_. This can be converted to ndjson (again) 

```bash
geo2topo -q 0 nw_europe.geojson > nw_europe.topojson 
ndjson-split 'd.objects.nw_albers.geometries' < nw_europe.topojson > nw_europe_1.ndjson
ndjson-split '[{"arcs" : d.arcs}]' < nw_europe.topojson >> nw_europe_1.ndjson
```

Producing a data stream with leading entries as follows: with a mixture of polygons and multipolygons but not grouped by nation.  
```json
{"type":"MultiPolygon","arcs":[[[0]],[[1]],[[2]]],"properties":{"nation":"Finland","id":12577865,"number":27153}}
{"type":"Polygon","arcs":[[3,4,5,6,7,8]],"properties":{"nation":"Austria","id":23424750,"number":8754413}}
{"type":"Polygon","arcs":[[9,10,11,12,13,14,15,16]],"properties":{"nation":"Belgium","id":23424757,"number":11491346}}
{"type":"Polygon","arcs":[[17,-6,18,-4,19,20,21]],"properties":{"nation":"Switzerland","id":23424957,"number":8236303}}
{"type":"MultiPolygon","arcs":[[[22]],[[23]],[[24]],[[25]],[[26,-8,27,-22,28,29,-12,30,31,32]],[[33]]],"properties":{"nation":"Germany","id":23424829,"number":80594017}}
{"type":"MultiPolygon","arcs":[[[34]],[[35]],[[36]],[[37]],[[38]],[[39]],[[40]],[[41]],[[42]],[[43]],[[44]],[[-33,45]]],"properties":{"nation":"Denmark","id":23424796,"number":5605948}}
{"type":"MultiPolygon","arcs":[[[46]],[[47]],[[48]],[[49,50]]],"properties":{"nation":"Estonia","id":23424805,"number":1251581}}
{"type":"MultiPolygon","arcs":[[[51]],[[52]],[[53]],[[54]],[[55]],[[56]],[[57]],[[58,59,60]]],"properties":{"nation":"Finland","id":23424812,"number":5491218}}
```
This can be remedied using _ndjson-reduce_ and then piped into a ruby script _mkvbx_ which splits multipolygons into polygons sharing common properties and produces a _vorobox_ file

```bash
ndjson-reduce < nw_europe_1.ndjson | mkvbx -m nation -s 1 > nw_europe_1.vxjson
```

where mkpoly is a simple ruby script to unpack multipolygons so all are now polygons. and to collect the resulting list into a _FeatureCollection_ to produce a complete geojson file.

And the resulting output 
