/**
 * PopG Web.
 * 
 * This is a JavaScript implementation of the PopG genetic simulation program:
 * 
 *      http://evolution.gs.washington.edu/popg/
 *
 *      Copyright 1993-2016. University of Washington and Joseph Felsenstein. 
 *      All rights reserved. Permission is granted to reproduce, perform, and 
 *      modify this program. Permission is granted to distribute or provide 
 *      access to this program provided that this copyright notice is not 
 *      removed, this program is not integrated with or called by any product
 *      or service that generates revenue, and that your distribution of this
 *      program is free. Any modified versions of this program that are 
 *      distributed or accessible shall indicate that they are based on this 
 *      program. Educational institutions are granted permission to distribute
 *      this program to their students and staff for a fee to recover distribution
 *      costs. Permission requests for any other distribution of this program 
 *      should be directed to license (at) u.washington.edu.
 *
 * Authors: Conor K. Camplisson, Kevin M. Schilling 
 * 
 */

// configure script
var NO_DRIFT_COLOR='#007bff';           // color for "No gene drift" ideal trace 
var TRACE_COLOR='#212529';              // color for simulated population traces
var AXIS_FONT_COLOR=TRACE_COLOR;        // axis font color
var AXIS_FONT_SIZE=22;                  // axis font size
var AXIS_FONT_FAMILY = 'Roboto';        // axis fonts
var AXIS_X_TITLE="Generation";          // X-axis label
var AXIS_Y_TITLE="P(A)";                // Y-axis label 
var SEED = NaN;                         // default seed for pseudo-random number generation
var FRAME_TIME = 10;
var FRAME_LIMIT = 100;

// execute when the DOM is fully loaded
$(function() {

    // configure website footer
    var d = new Date();
    $("#footer_wrapper").append(
        '<div class="footer"><center>Copyright &copy; ' + d.getFullYear()
        + ' | <a href="terms.html">Terms of Use</a>'
        + ' | <a target="_blank" href="https://github.com/basecon/popg-web">View on GitHub</a>'
    );

    // initialize config parameters
    var config = config_from_url();

    // update URL and control panel with latest paramters
    update_url(config);
    update_ctrl_panel(config);

    // update random number seed
    SEED = config.randSeed;

    // update URL on control panel form change
    $("#sim_form :input").change(function() {

        // disable continue button because params have changed        
        $('#continue_button').prop('disabled', true);

        // update config using changed form value
        config = config_from_url();
        config[$(this).attr('id')] = $(this).val();        

        // update URL    
        update_url(config);

    });    
    
    // configure default parameter buttons
    $('#defaults_button').click(function() {
        
        // disable continue button because params have changed        
        $('#continue_button').prop('disabled', true);
        
        // set default parameters        
        restore_defaults(); 

    });

    // hide loading icon and reveal control panel once plot os loaded
    $('.load_icon_container').hide();
    $('.load_hide').show();    

    // draw initial plot
    var run = new Run();

    // update plot dimensions when window is resized
    $(window).resize(function() {
        run.plot_result();
    });

    // run a new simulation on button click
    $("#sim_button").click(function() {
        run = new Run();
    });
    
    // continue existing simulation on button click
    $('#continue_button').click(function() {
        run.continue();
    });

});


/**
 * Run object definition.
 */
function Run() {

    // initialize this run
    this.init();
}
Run.prototype.init = function() {

    // temporarily disable plot buttons
    disable_buttons();

    // initialize config object with latest parameters
    this.config = config_from_url();
    SEED = this.config.randSeed;

    this.genArray = [];
    this.popArray = [];
    this.gensSoFar = 0;
    this.newRun = true;

    // run simulation
    this.result = this.calc_popg();
    
    // create plot
    this.plot_result()

};
Run.prototype.continue = function() {

    // run simulation
    this.result = this.calc_popg();
    
    // create plot
    this.plot_result();

};
Run.prototype.calc_popg = function() {   

    /*
    loop through all generations
        loop through all populations
            if migration. do that to edit p
            if not lost or fixed
                get mean fitness
                get freq of Aa and AA
                if anything but the starting piont (first gen)
                    calculate the next gen
                        use binomial
                if the starting point, calculate what average would be with infinite population
            if lost or fixed
                propagate this
            store population
        store generation
    */

    var config = this.config;
    var beginGen = this.gensSoFar;
    var endGen = beginGen + config.genRun;
    this.gensSoFar = this.gensSoFar + config.genRun; // Technically not true until the end of all gens in this run, shoudl move this down

    var currentPopSize = config.popSize * 2; // TODO look into pop gen (diploids, # ppl vs # chromosomes)
    
    if (this.newRun == true){
        // Initialize the big generation array with the starting values for each population
        var i;
        for (i = 0; i <= config.numPop; i++){
            this.popArray.push(config.initFreq);
        }
        // If the first run, add the starting values. If not, just continue with them
        this.genArray.push(this.popArray);
    }
    this.newRun = false;
   
    // Loop through all generations
    var pbar;
    var p;
    var q;
    var w;
    var pp1;
    var pp2;
    var nx;
    var ny;
    var numFixedPops;
    var numLostPops;
    var generation;
    for (generation = beginGen; generation < endGen; generation++){
        pbar = 0;
        var nextPopArray = [];
        this.popArray = this.genArray[generation];
        
        // Count number of populations no longer active
        for (var j = 1; j <= config.numPop; j++){ // Starting at 1 because our idealized population wont do this?
            var end = this.popArray[j];
            pbar += end;
            if (end <= 0.0) {
                numLostPops += 1;
            }            
            if (end >= 1.0) {
                numFixedPops += 1;
            }
        }
        pbar /= config.numPop;
        
        // Loop through all populations
        var population;
        for(population = 0; population <= config.numPop; population++){
            p = this.popArray[population];
            // All but population 0 get migrants
            if (population > 0){
                p = p * (1.0 - config.migRate) + config.migRate * pbar;
            }
            p = (1 - config.mutAa) * p + config.mutaA * (1-p);
            
            // if genotype is not fixed, calculate the new frequency
            if ((p > 0.0) && (p < 1.0)) {
                q = 1-p;
                w = (p * p * config.fitGenAA) + (2.0 * p * q * config.fitGenAa) + (q * q * config.fitGenaa); //get mean fitness
                pp1 = (p * p * config.fitGenAA) / w; //get frequency of AA after selection
                pp2 = (2.0 * p * q * config.fitGenAa) / w; //get frequency of Aa after selection                
                
                if (population > 0) { // calculate the next generation for real populations
                    nx = binomial(config.popSize, pp1); // how many AA survive

                    if (pp1 < 1.0 && nx < config.popSize) {// TODO: why no nested parens here, but yes above?
                        ny = binomial((config.popSize - nx), (pp2 / (1.0- pp1)));// how many Aa survive
                    }
                    else {
                        ny = 0;
                    }
                    // compute new number of A's and make it a frequency
                    nextPopArray.push(((nx * 2.0) + ny) / currentPopSize);
                }
                else { // calculate what the true average would be with an infinite population, for the "no-drift" idealized population with infinite size
                    nextPopArray.push(pp1 + (pp2 / 2.0));
                }
            }
            else {
                if (p <= 0.0) {
                    p = 0.0;
                }
                else {
                    p = 1.0;
                }
                nextPopArray.push(p);
            }
            
        } // end population
        this.genArray.push(nextPopArray);
        
    } // end generation
    
    // transpose data for plotting
    this.plot_genArray = this.genArray[0].map((col, i) => this.genArray.map(row => row[i]));

    // success
    return this.plot_genArray;

};
Run.prototype.plot_result = function() {

    // determine whether this is a new run
    var start_idx = 0;
    var first_run = this.result[0].length - 1 == this.config.genRun;
    if (first_run == false) {
        start_idx = this.result[0].length - this.config.genRun;
    }

    // construct animation frames, each frame is a complete set of traces
    var frames = [];
    var frame_threshold = Math.floor(this.result[0].length / FRAME_LIMIT);
    var frame_counter = 0;
    for (var frame_idx = start_idx; frame_idx < this.result[0].length; frame_idx++) {
        if (frame_counter < frame_threshold) {
            frame_counter += 1;
            continue;
        }
        frame_counter = 0;

        // construct all traces for this frame
        var traces = [];
        var x = xrange(0, frame_idx);
        for (var trace_idx = 0; trace_idx < this.result.length; trace_idx++) {

            // create a new trace
            var trace = {
                mode: 'lines',
                name: "",
                x: x,
                y: this.result[trace_idx].slice(0, frame_idx + 1),
                type: 'scatter',
                line: {
                  color: TRACE_COLOR,
                },                
                showlegend: false
            };

            // re-color and enable legend for no drift trace            
            if (trace_idx == 0) {
                trace.line.color = NO_DRIFT_COLOR;
                trace.name = "No gene drift";
                trace.showlegend = true;
            }            

            // add this trace            
            traces.push(trace);

        }
        
        // move first trace to the end
        traces.push(traces.shift());        
        
        // add this frame
        frames.push(traces);

    }

    // construct layouts with per-frame annotations
    var layouts = [];
    for (var i = 0; i < frames.length; i++) {

        // count number of populations that have fixed or lost the allele
        var num_fixed = 0;
        var num_lost = 0;
        for (var j = 0; j < frames[i].length - 1; j++) {
            var y_series = frames[i][j].y;
            var last_value = y_series.slice(-1)[0];
            if (last_value == 1) {
                num_fixed += 1;
            } else if (last_value == 0) {
                num_lost += 1
            }
        }

        // configure plot
        var layout = {
            yaxis: {
                range: [-0.05555555555555555, 1.0555555555555556],
                title: AXIS_Y_TITLE,
                titlefont: {
                    family: AXIS_FONT_FAMILY,
                    size: AXIS_FONT_SIZE,
                    color: AXIS_FONT_COLOR
                },
            },
            xaxis: {
                range: [0, this.result[0].length - 1],
                title: AXIS_X_TITLE,
                titlefont: {
                    family: AXIS_FONT_FAMILY,
                    size: AXIS_FONT_SIZE,
                    color: AXIS_FONT_COLOR
                },
            },
            hovermode: 'closest',
            autosize: true,
            margin: {
                l: 60,
                r: 10,
                b: 50,
                t: 30,
                pad: 0,
            },


            annotations: [
                {
                    text: 'Fixed: ' + parseInt(num_fixed),
                    align: 'left',
                    showarrow: false,
                    xref: 'paper',
                    yref: 'paper',
                    x: 0.98,
                    y: 0.965,
                    xshift: 100,
                },
                {
                    text: ' Lost: ' + parseInt(num_lost),
                    align: 'left',
                    showarrow: false,
                    xref: 'paper',
                    yref: 'paper',
                    x: 0.98,
                    y: 0.04,
                    xshift: 100,
                },
            ],
            legend: {
                xref: 'paper',
                yref: 'paper',
                x: 1.0,
                y: 0.9,
                yshift: 100,
            }
        };
        layouts.push(layout);
    }
    
    var config = {
        responsive: true,
        displaylogo: false,
        modeBarButtonsToRemove: [
            // 'toImage',
            'zoom2d', 
            'pan2d',
            'zoomIn2d',
            'zoomOut2d', 
            'autoScale2d',
            // 'resetScale2d',
            'toggleSpikelines',
            'hoverClosestCartesian',
            'hoverCompareCartesian',
        ],

        
    };

    // plot initial frame if new run, else update axes for continued sim
    if (first_run) {
        Plotly.newPlot('plot_div', frames[0], layouts[0], config);
    } else {
        Plotly.relayout('plot_div', layouts[0]);
        disable_buttons();
    }

    // add animated frames to plot
    for (var i = 0; i < frames.length; i++) {
        Plotly.animate('plot_div', {
                data: frames[i],
                layout: layouts[i]
            }, {
                frame: {
                    duration: FRAME_TIME,
                    redraw: true
                }
            }
        );
    }

    // re-enable control buttons after animation
    $('#plot_div').on('plotly_animated', function() {
        enable_buttons();
    });

};


/**
 * Default PopG algorithm parameters
 */
function get_default_config() {
    
    // set default values
    var config = {
        popSize: 100,
        fitGenAA: 1.0,
        fitGenAa: 1.0,
        fitGenaa: 1.0,
        mutAa: 0.0,
        mutaA: 0.0,
        migRate: 0.0,
        initFreq: 0.5,
        genRun: 100,
        numPop: 10,
        randSeed: NaN,     
    };
    
    // success
    return config;
    
}


/**
 * Helper functions to enable/disable simulation control buttons.
 */
function disable_buttons() {
    $('#continue_button').prop('disabled', true);
    $('#sim_button').prop('disabled', true);
}
function enable_buttons() {
    $('#continue_button').prop('disabled', false);
    $('#sim_button').prop('disabled', false);
}


/**
 * Update URL parameters using the current config object.
 */
function update_url(config) {
    var param_names = Object.keys(config);
    for (var i = 0; i < param_names.length; i++) {
        setUrlParam(param_names[i], config[param_names[i]]);
    }
}


/**
 * Update control panel values using the current config object.
 */
function update_ctrl_panel(config) {
    var param_names = Object.keys(config);
    for (var i = 0; i < param_names.length; i++) {
        $('#' + param_names[i]).val(config[param_names[i]]);
    }
}


/**
 * Load default parameters and update URL and control panel.
 */
function restore_defaults() {

    // load default params, update URL and control panel
    var config = get_default_config();
    update_url(config);
    update_ctrl_panel(config);

}


/**
 * Load default parameters and override with provided URL values, if any.
 */
function config_from_url() {
    
    // load default parameters and override with URL params as needed
    var config = get_default_config();
    var param_names = Object.keys(config);
    for (var i = 0; i < param_names.length; i++) {
        var url_value = getUrlParam(param_names[i]);
        if (url_value) {
            config[param_names[i]] = parseFloat(url_value);
        }
    }
    
    // success
    return config;
    
}


/**
 * Generates a seeded pseudo-random number.
 */
function random() {

    // generate a seeded pseudo-random number
    var rand;
    if (isNaN(SEED)) {
        rand = Math.random();
    } else {
        // adapted from: http://indiegamr.com/generate-repeatable-random-numbers-in-js/
        SEED = (SEED * 9301 + 49297) % 233280;
        rand = SEED / 233280;        
    }

    // success
    return rand;   

}


/**
 * Returns the number of successes for n trials,
 * with the probability of success in any given trial as p
 */
function binomial(n, p) {
    var successes = 0;
    var trial;
    for (trial = 0; trial < n; trial++){
        if (random() < p){
            successes++;
        }
    }
    return successes;  
}


/**
 * url query parameter manipulation from https://stackoverflow.com/a/40648561
 */
function setUrlParam(key, value) {

    var baseUrl = [location.protocol, '//', location.host, location.pathname].join(''),
        urlQueryString = document.location.search,
        newParam = key + '=' + value,
        params = '?' + newParam;

    // If the "search" string exists, then build params from it
    if (urlQueryString) {
        updateRegex = new RegExp('([\?&])' + key + '[^&]*');
        removeRegex = new RegExp('([\?&])' + key + '=[^&;]+[&;]?');
        if( typeof value == 'undefined' || value == null || (value == '' && typeof value != 'number') ) { // Remove param if value is empty
            params = urlQueryString.replace(removeRegex, "$1");
            params = params.replace( /[&;]$/, "" );
        } else if (urlQueryString.match(updateRegex) !== null) { // If param exists already, update it
            params = urlQueryString.replace(updateRegex, "$1" + newParam);
        } else { // Otherwise, add it to end of query string
            params = urlQueryString + '&' + newParam;
        }
    }
    window.history.replaceState({}, "", baseUrl + params);
};


/**
 * get url query param from https://html-online.com/articles/get-url-parameters-javascript/
 */
function getUrlVars() {
    var vars = {};
    var parts = window.location.href.replace(/[?&]+([^=&]+)=([^&]*)/gi, function(m,key,value) {
        vars[key] = value;
    });
    return vars;
}
function getUrlParam(parameter, defaultvalue){
    var urlparameter = defaultvalue;
    if(window.location.href.indexOf(parameter) > -1){
        urlparameter = getUrlVars()[parameter];
        }
    return urlparameter;
}


/**
 * Generates an array of integers from min to max.
 */
function xrange(min, max) {

    // generate array    
    var x = [];
    for (var i = min; i <= max; i++) {
        x.push(i);
    }

    // success
    return x;
    
}
