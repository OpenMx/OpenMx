// See http://scotch.io/tutorials/javascript/build-a-restful-api-using-node-and-express-4

var express    = require('express'); 		// call express
var app        = express(); 				// define our app using express
var bodyParser = require('body-parser');

// configure app to use bodyParser()
// this will let us get the data from a POST
app.use(bodyParser.urlencoded({ extended: true }));
app.use(bodyParser.json());

var port = process.env.PORT || 1337; 		// set our port

// ROUTES FOR OUR API
// =============================================================================
var router = express.Router(); 				// get an instance of the express Router

// test route to make sure everything is working (accessed at GET http://localhost:8080/api)
router.get('/', function(req, res) {
	res.json({ message: 'hooray! welcome to our api!' });	
});

var gModel = {};

router.post('/model/:name', function(req, res) {
    if (!(req.params.name in gModel)) {
	model = req.body
	model.evaluation = 0
	gModel[req.params.name] = model
	console.log('stored model ' + req.params.name);
	res.json({ok: true});
    } else {
	res.json({ok: false});
    }
})

router.get('/model', function(req, res) {
    var models = new Array();
    for (var foo in gModel) {
	models[models.length] = foo;
    }
    res.json({models: models});
})

router.get('/model/:name', function(req, res) {
    console.log('handing over model ' + req.params.name);
    res.json({model: gModel[req.params.name].model})
})

// -----------------

router.put('/model/:name/param', function(req, res) {
    model = gModel[req.params.name];
    model.param = req.body.param;
    model.evaluation += 1;
    model.fit = new Array();
    model.fitsum = 0;
    console.log('set ' + req.params.name + ' param to ' + model.param);
    gModel[req.params.name] = model;
    res.json({ok: true})
})

router.get('/model/:name/param', function(req, res) {
    model = gModel[req.params.name];
    console.log('hand over ' + req.params.name + ' param ' + model.param);
    res.json({evaluation: model.evaluation, at: model.param });
})

router.post('/model/:name/fit', function(req, res) {
    agent = req.body.agent;
    model = gModel[req.params.name];
    if (model.evaluation != req.body.evaluation || model.fit[agent]) {
	console.log("ignored fit reported for old evaluation " + req.body.evaluation);
	res.json({ok: false});
    } else {
	model.fit[agent] = true;
	model.fitsum += req.body.fit[0];
	console.log("got fit " + req.body.fit[0] + " from " + agent);
	gModel[req.params.name] = model;
	res.json({ok: true});
    }
})

router.get('/model/:name/fit', function(req, res) {
    model = gModel[req.params.name];
    var count = 0;
    for (var fit1 in model.fit) {
	if (fit1) count += 1;
    }
    if (count == 0) {
	res.json({evaluation: model.evaluation})
	return;
    }
    var mfit = model.fitsum / count;
    if (isNaN(mfit)) {
	res.json({evaluation: model.evaluation, nan: true });
    } else {
	res.json({evaluation: model.evaluation, fit: mfit });
    }
})

app.use('/api', router);

app.listen(port);
console.log('MIDDLE coordinator on port ' + port);
