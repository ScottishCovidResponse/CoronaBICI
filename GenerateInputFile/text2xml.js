const fs = require('fs');

// ----------------------------------------------------------------------------
// Importing existing javascript functionality from BICI GUI.

eval(fs.readFileSync('js/const.js').toString());
eval(fs.readFileSync('js/var.js').toString());
eval(fs.readFileSync('js/io.js').toString());
eval(fs.readFileSync('js/model.js').toString());
eval(fs.readFileSync('js/prior.js').toString());
eval(fs.readFileSync('js/import.js').toString());
eval(fs.readFileSync('js/init.js').toString());
eval(fs.readFileSync('js/convdata.js').toString());
eval(fs.readFileSync('js/inference.js').toString());

// ----------------------------------------------------------------------------
// Stubs of functions mocked out to remove dependencies on GUI code.

function alertimp(msg)                                     // Error message for imported file
{
    console.log("Error on line " + (jsto+1) + ": " + msg);
    process.exit(1)
}

function alertimp2(msg)                                     // Error message for imported file
{
    console.log("error: " + msg);
    process.exit(1)
}

function alertp(msg)
{
    console.log("error: " + msg);
    process.exit(1)
}

function textwidth(text)
{
    return 0;
}

function initcompmult()
{}

function startloading()
{}

function transplotinit()
{}

function buttoninit()
{}

function hex(c) {                                         // Converts decimal to hexidecimal
	var hex = (Math.floor(c)).toString(16);
	return hex.length == 1 ? "0" + hex : hex;
}

function alert(msg)
{
    console.log(msg);
}

startspawn = function() {}

ty = ""

// ----------------------------------------------------------------------------
// Processing command line arguments.

const usage = "usage: node " + process.argv[1] + " model_file output_file";

if (process.argv.length != 4) {
    console.log(usage);
    process.exit(1);
}

var model_file = process.argv[2];
var output_file = process.argv[3];

// ----------------------------------------------------------------------------
// Setting up minimal model parameters to allow for model to be imported from file.

initvar();

addclass('Age');
addclass('Time');

addcompartment(0, 'A0+', '', 0, 0, 0, 0, 0);
addcompartment(1, 'All', '', 0, 0, 0, 0, 0);

text = fs.readFileSync(model_file).toString()
importfile(text);

// ----------------------------------------------------------------------------
// Converting model to XML.

startinference(0, 1, output_file)