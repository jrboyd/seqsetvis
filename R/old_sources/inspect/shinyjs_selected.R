selected_jscode = "

shinyjs.getRange
gd.on('plotly_selected', (eventData) => {
  var xRange = eventData.range.x;
  var yRange = eventData.range.y;
  
  Plotly.relayout('graph', 'title',
                  `x range: [${xRange.map(formatter).join(', ')}]<br>
                    y range: [${yRange.map(formatter).join(', ')}]`
  );
});
"
selected_jscode_orig = "
var gd = document.getElementById('graph');
var d3 = Plotly.d3;

var formatter = d3.format('.2f');

Plotly.plot(gd, [{
mode: 'markers',
x: Array.apply(null, Array(100)).map(() => Math.random()),
y: Array.apply(null, Array(100)).map(() => Math.random()),
}], {
dragmode: 'select'
});

gd.on('plotly_selected', (eventData) => {
var xRange = eventData.range.x;
var yRange = eventData.range.y;

Plotly.relayout('graph', 'title',
`x range: [${xRange.map(formatter).join(', ')}]<br>
y range: [${yRange.map(formatter).join(', ')}]`
);
});
"


selected_html ="
  <head>
  <script src='https://cdn.plot.ly/plotly-latest.min.js'></script>
  </head>
  <body>
  <div id='graph'></div>
  </body>
"
