/*jslint browser:false, sub: true*/
/*globals Snap,alert,window,document,mina,renLib*/

/******************************************************************************************************************
 *
 * Program: script for LAMP RFLP web tool
 * Author: Caralyn Reisle
 * Last Update: 2014-March-25
 *
 * purpose: given lamp primers and a target sequence, predict fragment patterns for a given set of restriction enzymes
 *
 ******************************************************************************************************************/

//primer BIP = b1c-linker-b2
//primer FIP = f1c-linker-f2

//modules from notomi
// bplus = b1c-b2.....b1
// bminus = b1c....b2c-b1
// fminus = f1c....f2c-f1
// fplus = f1c-f2.....f1

//loop primers k Nagamine et al Molecular and Cellular Probes (2002) 16, 223-229 doi:10.1006/mcpr.2002.0415
//lb b1....-->.....b2 //on bminus only
//lf f2....<-......f1 //on fminus only

if (!('trim' in String.prototype)) {  
  String.prototype.trim = function() {
    return this.replace(/^\s*|\s*$/g,""); 
  };   
}

//var primer = "taccga";
var cent_ladder = [1517, 1200, 1000, 900, 800, 700, 600, 517, 500, 400, 300, 200, 100];
var fifty_ladder = [1350, 916, 766, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50];

var has_loop = false;

var run = 0; //keeps track of how many times the routine has been run for labelled the figures

/*****************************************************************
*
* function: strt()
* purpose: this is where the script begins
* input: none
* output: returns -1 if an error occurred. otherwise returns 0
*
*****************************************************************/
function strt(){
    run += 1;
    has_loop = false;
    var rxn = read_fasta();
    if (rxn==null){
	return -1;
    }
    //create the ren_list
    var layout = [];
    for(var i=1; i<=8; i++){
	var label = "lane"+i;
	var el = document.getElementById(label).childNodes[1];
	var item = el.options[el.selectedIndex].innerHTML;
	layout.push(item);
    }
    
    var result = idiot_proof_input(rxn.target, rxn);
    if (result===null){
	return -1;
    }
    var blocks = create_blocks(result[0], result[1], result[2]);
    var products = amplify_blocks(blocks);
    if (has_loop==true){
	products = loop_amp(products, rxn);
    }
    
    var rens = [get_ren("RsaI"), get_ren("BslI"), get_ren("HpaII"), get_ren("HinfI")];
    
    var lanes = find_cut_sites(products, layout);
    var uncut = generate_uncut(products);
    draw_gel(rxn, layout, lanes, uncut);
    return 0;
};

/*****************************************************************
*
* function: clear_results()
* purpose: removes all the current snap "pages" from the web page
* input: none
* output: none
*
*****************************************************************/
function clear_results(){
    elem = document.getElementsByTagName("svg");
    while(elem.length > 0){
        elem[0].parentNode.removeChild(elem[0]);
    }
    elem = document.getElementsByTagName("h2");
    while(elem.length > 0){
        elem[0].parentNode.removeChild(elem[0]);
    }
    run = 0;
};

/*****************************************************************
*
* function: loop_amp()
* purpose: amplifies the products using the loop primers, adds the
* 	new products to the current list of products
* input:
* 	products: the list of products generated from the base lamp primers b1, b2, f2 and f2
* 	rxn: the primers and sequence, i.e. the user input stored in a rxn object
* output: returns the new list of products
*
*****************************************************************/
function loop_amp(products, rxn){
    var result = [];
    result.concat(products);
    if(rxn.lf==null||rxn.lf.length!=0){
	for(var i=0; i<products.length; i++){
	    var temp = align_primer(products[i], rxn.lf);
	    var len = products[i].length;
	    if(temp.reverse==true){
		    var p = products[i].substring(0, temp.index-1)+rxn.lf;
		    result.push(p);
	    }else{
		    var st = temp.index+rxn.lf.length;
		    if (st<len-1) {
		        var p = rev_comp(rxn.lf)+products[i].substring(st);
		        result.push(p);
		    }
	    }
	}
	
    }
    if(rxn.lb==null||rxn.lb.length!=0){
	    for(var i=0; i<products.length; i++){
	        var temp = align_primer(products[i], rxn.lb);
	        if(temp.reverse==true){
	    	    var p = products[i].substring(0, temp.index-1)+rxn.lb;
	    	    result.push(p);
	        }else{
	    	    var st = temp.index+rxn.lb.length;
	    	    if (st<len-1){
	    	        var p = rev_comp(rxn.lb)+products[i].substring(st);
	    	        result.push(p);
	    	    }
	        }
	    }
    }
    return result;
};

/*****************************************************************
*
* function: read_fasta()
* purpose: splits the user input based on regular expressions to recognize fasta formatted text
* 	puts the split components into the rxn object. checks that the minimum elements have been added
* input: none
* output: the rxn object that stores the primers and target sequence information
*
*****************************************************************/
function read_fasta(){
    var fst = document.getElementById("fasta").value;
    var re = /\s*>\s*.*[\n\r]\s*[\w\s]*/gi;
    var matches = fst.match(re);
    var rxn = {b2: "", b1: "", f1: "", f2: "", lf: "", lb: "", target: "", target_name: "", linker: "tttt"};
    for(var i=0 ; i<matches.length; i++){
	    //get the tag name
	    var name = /\s*>\s*.*[\n\r]\s*/gi
	    var id = name.exec(matches[i]);
	    var split_point = id.index + id[0].length;
	    var tag = id[0].trim(); tag = tag.substring(1).toLowerCase();
	    var seq = matches[i].substring(split_point).trim().replace(/\s+/g, '').toLowerCase();
	    
	    switch(tag){
	        case "b1":
	    	    rxn.b1 = seq;
	    	    break;
	        case "b2":
	    	    rxn.b2 = seq;
	    	    break
	        case "f1":
	    	    rxn.f1 = seq;
	    	    break;
	        case "f2":
	    	    rxn.f2 = seq;
	    	    break;
	        case "lf":
	    	    rxn.lf = seq;
	    	    has_loop = true;
	    	    break;
	        case "lb":
	    	    rxn.lb = seq;
	    	    has_loop = true;
	    	    break;
	        case "linker":
	    	    rxn.linker = seq;
	    	    break;
	        default:
	    	    rxn.target = seq;
	    	    rxn.target_name = tag;
	    	    break;
	    }
    }
    if (rxn.b2==""||rxn.b1==""||rxn.f1==""||rxn.f2==""||rxn.target=="") {
	    console.log("error. base components required for reaction. b1, b2, f1, f2 primers and a single target");
	    return null;
    }
    return rxn;
}

/*****************************************************************
*
* function: idiot_proof_input()
* purpose: aligns the primers and then makes sure they were entered in the
* 	correct order/orientation and that the rxn components will permit a lamp reaction.
* 	if the best matches of the primers do not correspond to the correct positioning for a lamp reaction
* 	then this will return null
* input:
* 	str: the target sequence
* 	primers: the rxn object that stores the rxn information, primers etc
* 	
* output:
* 	null if the rxn will not permit a lamp of the target
* 	else an array containing {str, matches, primers}
* 		str:  being the rxn target (may have been updated)
* 		matches being the computed alignment positions (the first index in the target where each primer would anneal)
* 		primers: the list of primers (may have been updated)
*
*****************************************************************/
function idiot_proof_input(str, primers){
    var b2 = align_primer(str, primers.b2);
    var f2 = align_primer(str, primers.f2);
    
    if(b2.index<f2.index){
        if(b2.reverse==true){
            console.log("input string was the other strand from what was needed. will take the compliment and re-compute");
            str = rev_comp(str, false); //just gives the compliment, not reversed
            b2 = align_primer(str, primers.b2);
            f2 = align_primer(str, primers.f2);
        }else{
            console.log("should be the expected format");
        }
    }else{
        if(b2.reverse==true){
            console.log("need reverse compliment of the input strand. generating...");
            str = rev_comp(str, true); //just gives the compliment, not reversed
            b2 = align_primer(str, primers.b2);
            f2 = align_primer(str, primers.f2);
        }else{
            console.log("just need to reverse the input strand");
            str = str.reverse();
            b2 = align_primer(str, primers.b2);
            f2 = align_primer(str, primers.f2);
        }
    }
    
    //after this all should be b2.....f2c
    
    var temp = str.substring(b2.index+primers.b2.length, f2.index);
    if(temp.length<primers.b1.length){
        console.log("error in the primer input. please check the provided sequences");
        alert("error in the primer input. primers do not align to the target in the expected format. please check the provided sequences");
        return null;
    }
    var b1 = align_primer(temp, primers.b1);
    b1.index += b2.index+primers.b2.length;
    temp = str.substring(b1.index+primers.b1.length, f2.index);
    if(temp.length<primers.b1.length){
        console.log("error in the primer input. please check the provided sequences");
        alert("error in the primer input. primers do not align to the target in the expected format. please check the provided sequences");
        return null;
    }
    var f1 = align_primer(temp, primers.f1);
    f1.index += b1.index+primers.b1.length;
    //all the indicies must be right due to our substring restriction. now check the directions
    //already know that b2 is not reversed and order is b2...b1....f1....f2
    if(b1.reverse==false){ //should be true, revcomp the primer
        primers.b1 = rev_comp(primers.b1, true);
        b1.reverse = true;
    }
    if(f1.reverse==true){ //should be false
        primers.f1 = rev_comp(primers.f1, true);
        f1.reverse = true;
    }
    if(f2.reverse==false){ //should be true
        primers.f2 = rev_comp(primers.f2, true);
        f2.reverse = true;
    }
    
    var matches = {b1: b1, b2: b2, f1: f1, f2: f2};
    return [str, matches, primers];
    
};

/*****************************************************************
*
* function: create_blocks()
* purpose: creates the lamp reaction building blocks using the previously found
* 	alignment positions of the primers and the pre-defined template take from
* 	the notomoi et al (2000) paper
* input:
* 	str: the target sequence
* 	primers: the rxn object that stores the rxn information, primers etc
* 	matches: the computed alignment index positions on the target sequence for each primer
* output:
* 	the building blocks of the lamp reaction
*
*****************************************************************/
function create_blocks(str, matches, primers){
    var blocks = {};
    blocks.bplus = primers.b1+primers.linker+primers.b2+str.substring(matches.b2.index+primers.b2.length, matches.b1.index+primers.b1.length);
    blocks.bminus = rev_comp(blocks.bplus, true);
    blocks.fplus = str.substring(matches.f1.index, matches.f2.index-1)+primers.f2+primers.linker+primers.f1;
    blocks.fminus = rev_comp(blocks.fplus, true);
    blocks.plus = str.substring(matches.b1.index+primers.b1.length, matches.f1.index-1);
    blocks.minus = rev_comp(blocks.plus, true);
    return blocks;
};

/*****************************************************************
*
* function: amplify_blocks()
* purpose: concatanates the building bloacks in a pre-defined pattern to predict the lamp products
* input:
* 	b: the lamp reaction building blocks
* output:
* 	the products of the lamp reaction
*
*****************************************************************/
function amplify_blocks(b){
    products = [];
    var temp = b.bplus+b.plus+b.fplus; //6
    products.push(temp);
    temp += b.minus+b.bminus; //7
    products.push(temp);
    temp += b.plus+b.fminus+b.minus+b.bminus; //9
    products.push(temp);
    temp += b.plus+b.fplus+b.minus+b.bplus+b.plus+b.fminus+b.minus+b.bminus; //13
    products.push(temp);
    temp += b.plus+b.fplus+b.minus+b.bminus+b.plus+b.fminus+b.minus+b.bplus+b.plus+b.fplus+b.minus+b.bplus+b.plus+b.fminus+b.minus+b.bminus; //15
    products.push(temp);
    temp = b.fminus+b.minus+b.bplus+b.plus+b.fplus+b.minus+b.bplus+b.plus+b.fminus+b.minus+b.bminus; //17
    products.push(temp);
    temp = b.fminus+b.minus+b.bminus; //20
    products.push(temp);
    temp += b.plus+b.fplus; //10
    products.push(temp);
    temp = b.bplus+b.plus+b.fminus+b.minus+b.bminus; //16
    products.push(temp);
    temp = b.bplus+b.plus+b.fplus+b.minus+b.bplus+b.plus+b.fminus+b.minus+b.bminus; // 18
    products.push(temp);
    temp += b.plus+b.fplus+b.minus+b.bminus+b.plus+b.fminus+b.minus+b.bminus; //19
    products.push(temp);
    
    return products;
};

/*****************************************************************
*
* function: generate_uncut()
* purpose: when a restriction enzyme does not have any cut sites on the products,
* 	this gives the lengths of the uncut fragments since that is what we will see in the gel if the product is undigested
* input:
* 	products: the lamp reaction products
* output:
* 	the lengths of the fragments of these products
*
*****************************************************************/
function generate_uncut(products){
    var uncut = [];
    for(var i=0; i<products.length; i++){
	    uncut.push(products[i].length);
    }
    return uncut;
}

/*****************************************************************
*
* function: find_cut_sites()
* purpose: looks for cut site matches for a list of restriction enzymes and generates the fragment lengths where found
* input:
* 	products: the lamp reaction products
* 	layout: from the user input, the list of what will be in each "lane" of the gel
* output:
* 	the fragments that will be in each lane in order of the output gel
* 	note that this can contain ladders and uncut products as well as digested products
*
*****************************************************************/
function find_cut_sites(products, layout){
    var uncut = generate_uncut(products);
    var lanes = []
    for(var i=0; i<layout.length; i++){
        var fragments;
        var temp;
	    switch(layout[i]){
	        case "100 bp ladder":
	    	    var temp = cent_ladder.concat(cent_ladder);
	    	    temp = temp.concat(temp);
	    	    lanes.push(temp);
	    	    layout[i] = "L: 100bp";
	    	    break;
	        case "50 bp ladder":
	    	    var temp = fifty_ladder.concat(fifty_ladder);
	    	    temp = temp.concat(temp);
	    	    lanes.push(temp);
	    	    layout[i] = "L: 50bp";
	    	    break;
	        case "None":
	    	    lanes.push(uncut);
	    	    break;
	        default:
	    	    var ren = get_ren(layout[i]);
	    	    if((temp=ren_palindrome(ren))===null){
	    	        fragments = cut_products(products, ren);
	    	    }else{
	    	        fragments = cut_products(products, ren);
	    	        fragments.concat(cut_products(products, temp));
	    	    }
	    	    if(fragments.length==0){
	    	        lanes.push(uncut);
	    	    }else{
	    	        lanes.push(fragments);
	        	}
	    	    break;
	    }
    }
    return lanes;
};

/*****************************************************************
*
* function: draw_gel()
* purpose: generates the GUI output of the digestion reactions
* input:
* 	lanes: the fragments in each lane
* 	layout: from the user input, the list of what will be in each "lane" of the gel
* 	rxn: the reaction information, used for labelling the gel/GUI output
* 	uncut: list of uncut products so we can find max values
* output:
* 	outputs the predicted "gel"/GUI output to the webpage
*
*****************************************************************/
function draw_gel(rxn, layout, lanes, uncut){
    var w = window.innerWidth;
    var h = window.innerHeight;
    var max = Math.log(Math.max.apply(Math, uncut))/Math.LN10+0.25; //maximum log value
    
    var min = Math.log(100)/Math.LN10;
    for(var i=0; i<lanes.length; i++){
	    var temp = Math.log(Math.min.apply(Math, lanes[i]))/Math.LN10;
	    if (temp<min) {
	        min = temp;
	    }
    }
    min -= 0.5;
    var increment = h/(max-min); //intervals
    var paper = new Snap(w, h);
    paper.rect(0, 0, w-15, h).attr({"stroke-width": 2, fill: "white", stroke: "black"}); //create border for snap space
    
    var lane_width = w/(lanes.length+2.5); //two for spacing
    var offset = (lane_width*2)/(lanes.length); //spacing between lanes
    for(var i=0; i<lanes.length; i++){
	    var x = offset+lane_width*i;
	    for(var j=0; j<lanes[i].length; j++){
	        var y = h - (Math.log(lanes[i][j])/Math.LN10 - min)*increment;
	        //draw a rectangle
	        var opaq = (lanes[i][j]/1000)*1;
	        if (opaq>1) {
		        opaq = 1;
	        }
	        paper.rect(x, y, lane_width, 1).attr({"stroke-opacity": opaq, "fill-opacity": opaq}); //
	    }
	    paper.text(x+lane_width/2, (max-min)*increment-15, layout[i]).attr("text-anchor", "middle");
	    offset += lane_width/lanes.length;
    }
    
    //draw the dragable line...
    //http://svg.dabbles.info/snaptut-drag.html
    var measure = paper.group();
    var ybp100 = h -(Math.log(100)/Math.LN10 - min)*increment;
    
    measure.add(paper.rect(10,ybp100, lane_width*lanes.length+offset-10, 2).attr({fill: "blue", stroke: "blue", "stroke-opacity": 0, "fill-opacity": 0.5, "stroke-width": 4}));
    var label = paper.text(offset+lane_width*lanes.length, ybp100+5, "100.00 bp");
    measure.add(label);
    
    /*dragging functions: move, start, stop
    reference for dragging function: http://svg.dabbles.info/snaptut-drag.html
    note that this will have to be edited to use on a group that does not contain a text element*/
    var move = function(dx, dy, x, y){
	    this.attr({ transform: this.data('origTransform') + (this.data('origTransform') ? "T" : "t") + [0, dy]});
	
	    //code to deal with updating our bp label on the text element of our group
	    var original_y = this[0].node.attributes[1].value;
	    var shift = this.matrix.f;
	    var yy = parseFloat(original_y)+ parseFloat(shift);
	    yy = Math.pow(10, (h-yy)/increment+min); // (h - y)/increment = Math.log(lanes[i][j])/Math.LN10 - min;
	    this[1].node.textContent = yy.toFixed(2) +" bp";
    }
    var start = function(){
	    console.log("started dragging");
	    this.data('origTransform', this.transform().local);
    }
    var stop = function(){
	    console.log("finished dragging");
    }
	
    measure.drag(move, start, stop);
    
    var figure_title = document.createElement('h2');
    figure_title.innerHTML = "Figure "+run+". Restriction Analysis of "+rxn.target_name+" LAMP products";
    document.body.appendChild(figure_title);
}

/*****************************************************************
*
* function: ren_palindrome()
* purpose: checks if the REN we are looking at has a reverse palindrome so we
* 	know if we need to check the reverse compliment of the cut site or not
* input:
* 	ren: the restriction enzym
* output:
* 	returns null if it is not a reverse comp palindrome, otherwise returns a new ren with the new revcomp cutsite
*
*****************************************************************/
//find if its a palindrome, if not return a reverse comp version
function ren_palindrome(ren){
    var site = rev_comp(ren.cutsite);
    var halfway = Math.ceil(ren.cutsite.length/2.0);
    if(site.localeCompare(ren.cutsite)&&ren.cutindex==halfway){
        return null;
    }
    var result = {name: ren.name, cutindex: ren.cutsite.length - ren.cutindex, cutsite: site};
    return result;
};

/*****************************************************************
*
* function: cut_products()
* purpose: splits the products by ren cut sites
* input:
* 	ren: the restriction enzym
* 	products: the products from the lamp reaction
* output:
* 	returns null if there are no matching cut sites found, else returns the fragments
*
*****************************************************************/
function cut_products(products, ren){
    var re = generate_regex(ren);
    var fragments = [];
    for(var i=0; i<products.length; i++){
        var temp = products[i];
        var result;
        var indicies = [];
        while ((result = re.exec(temp)) !== null){
            var index = result.index;
            indicies.push(index);
        }
        var prev = 0;
	    if(indicies.length==0){
		    fragments.push(products[i].length);
	    }
        for(var j=0; j<indicies.length; j++){
            var f = indicies[j]+ren.cutindex-prev;
            fragments.push(f);
            prev = indicies[j];
        }
    }
    if(fragments.length==0){
	    return null;
    }
    //make a new array without the duplicates
    //fragments = remove_duplicates(fragments);
    return fragments;
};

/*****************************************************************
*
* function: remove_duplicates()
* purpose: removes duplicate elements from an array of integers
* input:
* 	arr: the array we wish to remove duplicates from
* output:
* 	returns the sorted array of distinct elements
*
*****************************************************************/
function remove_duplicates(arr){
    var result = [];
    arr.sort();
    var prev = arr[0];
    result.push(prev);
    for(var i=1; i<arr.length; i++){
        if(arr[i]!==prev){
            prev = arr[i];
            result.push(prev);
        }
    }
    return result;
};

/*****************************************************************
*
* function: generate_regex()
* purpose: generats a regular expression for a restriction enzyme cutsite
* input:
* 	ren: the restriction enzyme
* output:
* 	returns the regular expression for this ren
*
*****************************************************************/
function generate_regex(ren){
    var temp = ren.cutsite.toLowerCase();
    var reg = [];
    for(var i=0; i<temp.length; i++){
        var str = "";
        switch(temp.charAt(i)){
            case 'a':
            case 't':
            case 'c':
            case 'g':
                str = temp.charAt(i);
                break;
            case 'm':
                str = '[ac]'; break;
            case 'k':
                str = '[gt]'; break;
            case 'r':
                str = '[ag]'; break;
            case 'y':
                str = '[ct]'; break;
            case 'w':
                str = '[at]'; break;
            case 's':
                str = '[gc]'; break;
            case 'h':
                str = '[act]'; break;
            case 'v':
                str = '[acg]'; break;
            case 'n':
                str = '[atcgmkrywshvnbd]'; break;
            case 'b':
                str = '[cgt]'; break;
            case 'd':
                str = '[agt]'; break;
            default:
                console.log("error, unrecognized character");
                break;
        }
        reg.push(str);
    }
    var str = reg.join("");
    var re = new RegExp(str, "gi");
    return re;
};

/*****************************************************************
*
* function: align_primer()
* purpose: find the best alignment position for an primer on a target
* input:
* 	str; the target string
* 	prm: the primer string (prm.length<str.length)
* output:
* 	returns a match object with the index position and whether or not
* 	the best match was found with the original primer string or the reverse complement
*
*****************************************************************/
function align_primer(str, prm){
    var revcomp = rev_comp(prm, true);
    var minscore = prm.length; //min number of differences btwn primer and sequence
    var reverse = false; //min number of differences btwn revcomp and sequence
    var index = 0;
    for(var i=0; i<str.length-prm.length; i++){
        var score = 0;
        var revscore = 0;
        for(var j=0; j<prm.length; j++){
            if(str.charAt(i+j)!=prm.charAt(j)){
                score++;
            }
            if(str.charAt(i+j)!=revcomp.charAt(j)){
                revscore++;
            }
            if(score>minscore&&revscore>minscore){
                break;
            }
        }
        if(score<=revscore&&score<minscore){
            minscore = score; reverse = false; index = i;
        }else if(revscore<score&&revscore<minscore){
            minscore = revscore; reverse = true; index = i;
        }
    }
    return {index: index, reverse: reverse};
};

/*****************************************************************
*
* function: rev_comp()
* purpose: generate the reverse compliment or just compliment of a sequence
* input:
* 	str: the target string
* 	rev: boolean, if true: will generate revcomp, if false: will just generate the compliment
* output:
* 	returns the compliment or reverse compliment of the input str depending on the boolean rev parameter
*
*****************************************************************/
function rev_comp(str, rev){
    if (str==null) {
	    return null
    }
    str = str.toLowerCase();
    var result = [];
    for(var k=0; k<str.length; k++){
        var i = k;
        if(rev==true){
            i = str.length - k -1;
        }
        switch(str.charAt(i)){
            case 'c':
                result.push('g');
                break;
            case 'g':
                result.push('c');
                break;
            case 't':
                result.push('a');
                break;
            case 'a':
                result.push('t');
                break;
	        case 'r': //a, g --> t, c
		    	result.push('y');
	            break;
	        case 'y': //c, t --> g, a
		    	result.push('r');
	            break;
	        case 'k': //g, t --> c, a
		    	result.push('m');
	            break;
	        case 'm': //a, c --> t, g
		    	result.push('k');
	            break;
	        case 'b': //c, g, t --> g, c, a
		    	result.push('v');
	            break;	
	        case 'd': //a, g, t --> t, c, a
		    	result.push('h');
	            break;
	        case 'h': //a, c, t --> t, g, a
		    	result.push('d');
	            break;
	        case 'v': //a, c, g --> t, g, c
		    	result.push('b');
	            break;
            default: //s, w, n
                result.push(str.charAt(i));
                break;
        }
    }
    return result.join("");
};

/*****************************************************************
*
* function: get_ren()
* purpose: get the information for a given Restriction enzyme (ren)
* input:
* 	name: the name of the ren
* output:
* 	returns the values associated with that ren, i.e. cutsite, cutindex, etc
*
*****************************************************************/
function get_ren(name){
	for(var i=0; i<renlib.length; i++){
		if(renlib[i].name.toLowerCase()==name.toLowerCase()){
			return renlib[i];
		}
	}
	return null;
};


  
