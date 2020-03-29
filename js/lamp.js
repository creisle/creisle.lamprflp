/*jslint browser:false, sub: true*/
/*globals Snap,alert,window,document,createLampReaction*/

if (!('trim' in String.prototype)) {
    String.prototype.trim = function() {
        return this.replace(/^\s*|\s*$/g,'');
    };
}


let run = 0; //keeps track of how many times the routine has been run for labelled the figures

function startLampReaction(){
    run += 1;
    const fastaContent = document.getElementById('fasta').value;
    //create the ren_list
    const layout = [];

    for (let i = 1; i <= 8; i++){
        const label = 'lane' + i;
        const {childNodes: [, el]} = document.getElementById(label);
        const item = el.options[el.selectedIndex].innerHTML;
        layout.push(item);
    }

    try {
        const {lanes, undigestedProduct, reaction} = createLampReaction(fastaContent,  layout);
        drawGel(reaction.targetName, layout, lanes, undigestedProduct);
    } catch (err) {
        alert(err);
    }
}

function clearReactionResults(){
    let elem = document.getElementsByTagName('svg');

    while (elem.length > 0){
        elem[0].parentNode.removeChild(elem[0]);
    }
    elem = document.getElementsByTagName('h2');

    while (elem.length > 0){
        elem[0].parentNode.removeChild(elem[0]);
    }
    run = 0;
}

/*****************************************************************
*
* function: draw_gel()
* purpose: generates the GUI output of the digestion reactions
* input:
*     lanes: the fragments in each lane
*     layout: from the user input, the list of what will be in each "lane" of the gel
*     rxn: the reaction information, used for labelling the gel/GUI output
*     uncut: list of uncut products so we can find max values
* output:
*     outputs the predicted "gel"/GUI output to the webpage
*
*****************************************************************/
function drawGel(name, layout, lanes, uncut){
    var w = window.innerWidth;
    var h = window.innerHeight;
    var max = Math.log(Math.max.apply(Math, uncut)) / Math.LN10 + 0.25; //maximum log value

    var min = Math.log(100) / Math.LN10;

    for (let i = 0; i < lanes.length; i++){
        var temp = Math.log(Math.min.apply(Math, lanes[i])) / Math.LN10;

        if (temp < min) {
            min = temp;
        }
    }
    min -= 0.5;
    var increment = h / (max - min); //intervals
    var paper = new Snap(w, h);
    paper.rect(0, 0, w - 15, h).attr({'stroke-width': 2, fill: 'white', stroke: 'black'}); //create border for snap space

    var lane_width = w / (lanes.length + 2.5); //two for spacing
    var offset = (lane_width * 2) / (lanes.length); //spacing between lanes

    for (let i = 0; i < lanes.length; i++){
        var x = offset + lane_width * i;

        for (var j = 0; j < lanes[i].length; j++){
            var y = h - (Math.log(lanes[i][j]) / Math.LN10 - min) * increment;
            //draw a rectangle
            var opaq = (lanes[i][j] / 1000) * 1;

            if (opaq > 1) {
                opaq = 1;
            }
            paper.rect(x, y, lane_width, 1).attr({'stroke-opacity': opaq, 'fill-opacity': opaq}); //
        }
        paper.text(x + lane_width / 2, (max - min) * increment - 15, layout[i]).attr('text-anchor', 'middle');
        offset += lane_width / lanes.length;
    }

    //draw the dragable line...
    //http://svg.dabbles.info/snaptut-drag.html
    var measure = paper.group();
    var ybp100 = h - (Math.log(100) / Math.LN10 - min) * increment;

    measure.add(paper.rect(10,ybp100, lane_width * lanes.length + offset - 10, 2).attr({fill: 'blue', stroke: 'blue', 'stroke-opacity': 0, 'fill-opacity': 0.5, 'stroke-width': 4}));
    var label = paper.text(offset + lane_width * lanes.length, ybp100 + 5, '100.00 bp');
    measure.add(label);

    /*dragging functions: move, start, stop
    reference for dragging function: http://svg.dabbles.info/snaptut-drag.html
    note that this will have to be edited to use on a group that does not contain a text element*/
    var move = function(dx, dy){
        this.attr({ transform: this.data('origTransform') + (this.data('origTransform')
            ? 'T'
            : 't'
        ) + [0, dy]});

        //code to deal with updating our bp label on the text element of our group
        var original_y = this[0].node.attributes[1].value;
        var shift = this.matrix.f;
        var yy = parseFloat(original_y) + parseFloat(shift);
        yy = Math.pow(10, (h - yy) / increment + min); // (h - y)/increment = Math.log(lanes[i][j])/Math.LN10 - min;
        this[1].node.textContent = yy.toFixed(2) + ' bp';
    };

    var start = function(){
        this.data('origTransform', this.transform().local);
    };

    measure.drag(move, start);

    var figure_title = document.createElement('h2');
    figure_title.innerHTML = 'Figure ' + run + '. Restriction Analysis of ' + name + ' LAMP products';
    document.body.appendChild(figure_title);
}

