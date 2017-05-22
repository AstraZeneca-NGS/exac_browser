
var make_cnvs_svg = function(_cnvs, _transcript, scale_type, skip_utrs, container, cnv_svg)
{
// plot CNVs as separate function
// implemented this way so that just the CNVs svg can be redrawn without fudging with any of the other bits!

    var coding_coordinate_params = get_coding_coordinate_params(_transcript, skip_utrs);
    var chart_width;
    if (scale_type == 'overview') {
        chart_width = gene_chart_width;
    } else {
        chart_width = coding_coordinate_params.size*2;
    }

    var exon_x_scale = d3.scale.linear()
        .domain([0, coding_coordinate_params.size])
        .range([0, chart_width]);


    var cnv_filter_status = $('#filtered_checkbox').is(":checked");
    var m = get_max_cnv(_cnvs, cnv_filter_status);

    var svg = d3.select(container).append("svg")
        .attr("width", chart_width + cnv_chart_margin.left + cnv_chart_margin.right)
        .attr("height",100 + 40)
        .attr('id', cnv_svg)
        .append("g")
        .attr('id', 'track')
        .attr("transform", "translate(" + cnv_chart_margin.left + "," + 0 + ")");

    var cnv_scale = d3.scale.linear()
         .domain([(-1*(m+1)),(m+1)])
         .range([130,30]);


    var yAxis = d3.svg.axis()
        .scale(cnv_scale)
        .orient("left")
        .tickFormat(function(d,i){return Math.abs(d)})
        .ticks(5);

    var del_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        var x = get_cnv_pop(_cnvs, d.start, d.stop, 'del', cnv_filter_status)
                    .replace(/,/g, '<br/>')
                    .replace(/:/g,': ')
                    .replace(/\//g,' / ')
                    .replace('|', '<br/>');
        return x;
    });
    svg.call(del_tip);

    var dup_tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        var x =  get_cnv_pop(_cnvs, d.start, d.stop, 'dup', cnv_filter_status)
                    .replace(/,/g, '<br/>')
                    .replace(/:/g,': ')
                    .replace(/\//g,' / ')
                    .replace('|', '<br/>');
        return x;
    });
    svg.call(dup_tip);

    //deletions
        svg.selectAll("bar")
        .data(_transcript.exons)
            .enter()
            .append("rect")
            .call(yAxis)
            .attr('pointer-events', 'all')
            .on('mouseover', function(d) {
                del_tip.show(d);
            })
            .on('mouseout', function(d) {
                del_tip.hide(d);
            })
            .attr('class', 'track_bar')
            .style("fill", "#cd2932")
            .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(_transcript, d.start, skip_utrs)); })
            .attr("y", function(d, i) { return cnv_scale(0);})
            .attr("width", function(d, i) {
                if (get_coding_coordinate(_transcript, d.start, true) == -100) {
                    return exon_x_scale(175);
                }
                return exon_x_scale(d.stop-d.start+1);
            })
            .attr("height", function(d, i) {
                if (d.feature_type == 'CDS') {
                    return cnv_scale(0) - cnv_scale(get_cnv(_cnvs, d.start, d.stop, 'del', cnv_filter_status));
                }
            });

    // duplications
    svg.selectAll("bar")
            .data(_transcript.exons)
            .enter()
            .append("rect")
            .attr('pointer-events', 'all')
            .on('mouseover', function(d) {
                dup_tip.show(d);
            })
            .on('mouseout', function(d) {
                dup_tip.hide(d);
            })
            .attr('class', 'track_bar')
            .style("fill", "#a96500")
            .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(_transcript, d.start, skip_utrs)); })
            .attr("y", function(d, i) { return cnv_scale(get_cnv(_cnvs, d.start, d.stop, 'dup', cnv_filter_status));})

            .attr("width", function(d, i) {
                if (get_coding_coordinate(_transcript, d.start, true) == -100) {
                    return exon_x_scale(175);
                }
                return exon_x_scale(d.stop-d.start+1);
            })
            .attr("height", function(d, i) {
                if (d.feature_type == 'CDS') {
                    <!-- console.log(d); -->
                    <!-- console.log(get_cnv(_cnvs, d.start, d.stop, 'dup', cnv_filter_status)); -->
                    <!-- console.log(cnv_scale(0)); -->
                    <!-- console.log(get_cnv(_cnvs, d.start, d.stop, 'dup', cnv_filter_status)); -->
                    <!-- console.log(cnv_scale(0) - cnv_scale(get_cnv(_cnvs, d.start, d.stop, 'dup', cnv_filter_status))); -->

                    return cnv_scale(0) - cnv_scale(get_cnv(_cnvs, d.start, d.stop, 'dup', cnv_filter_status));
                }
            });


    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    svg.append("text")
        .attr("x", (chart_width / 2))
        .attr("y", 0 + (cnv_chart_margin.top / 10)*8)
        .attr("text-anchor", "middle")
        .style("font-size", "14px")
        .style("font-weight", "bold")
        .text("CNV Counts")

    svg.append("text")
        .attr("x", (chart_width / 2) + 120)
        .attr("y", 0 + (cnv_chart_margin.top / 10)*8)
        .attr("text-anchor", "middle")
        .attr("fill", "#428bca")
        .style("font-size", "12px")
        .style("cursor", "pointer")
        .style("font-weight", "bold")
        .style("font-family", "FontAwesome")
        .html("<a>(view individual CNVs &#xf08e;)</a>")
        .on("click", function() {
            window.open("http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=chr{{ gene.chrom }}%3A{{ gene.start - 1 }}-{{ gene.stop - 1 }}&hgt.customText=http://personal.broadinstitute.org/ruderfer/exac/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed", "_blank");
        });



    //plotting lines underneath bars
        svg.selectAll("bar")
            .data(_transcript.exons)
            .enter()
            .append("rect")
            .attr('class', 'track_bar')
            .style("fill", "white")
            .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(_transcript, d.start, skip_utrs)); })
            .attr("y", function(d, i) {
                if (d.feature_type == 'CDS') {
                    return cnv_scale(0);
                }
            })
            .attr("width", function(d, i) {
                if (get_coding_coordinate(_transcript, d.start, true) == -100) {
                    return exon_x_scale(175);
                }
                return exon_x_scale(d.stop-d.start+1);
            })
            .attr("height", function(d, i) {
                if (d.feature_type == 'CDS') {
                    return 1;
                }
            });
}


function gene_chart(data, new_data, variant_data, _transcript, _cnvs, container, cnv_svg, genome, project_name) {
    var coords = 'pos_coding_noutr';
    var coding_coordinate_params = get_coding_coordinate_params(_transcript, true);
    var chart_width = gene_chart_width;
    var metric = 'mean';

    var exon_x_scale = d3.scale.linear()
        .domain([0, coding_coordinate_params.size])
        .range([0, chart_width]);

    var max_cov = (metric == 'mean' || metric == 'median') ? d3.max(data, function(d) { return d[metric]; }) : 1;
    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([gene_chart_height, 0]);

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    var svg = d3.select(container).append("svg")
        .attr("width", chart_width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .attr('id', 'inner_svg')
        .attr('class', 'hidden-xs')
        .append("g")
        .attr('id', 'inner_graph')
        .attr("transform", "translate(" + gene_chart_margin.left + "," + gene_chart_margin.top + ")");

    var area = d3.svg.area()
        .x( function(d) {
            return exon_x_scale(d[coords]);
        }).y0( function(d) {
            return gene_chart_height;
        }).y1( function(d) {
            return (metric in d) ? y(d[metric]) : gene_chart_height;
        });

    svg.append("path")
        .datum(new_data)
        .style("fill", "steelblue")
        .attr('class', 'area')
        .attr("d", area);

    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis);

    svg.selectAll("circle")
        .data(data)
        .enter()
        .append("circle")
        .attr("cx", function(d) {
            if (d['pos_coding_noutr'] == undefined) {
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        })
        .attr("cy", gene_chart_height + 3)
        .attr('r', 2 );

    d3.select(container).append("br"); //make sure track_svg is below graph and not to the right of it

    // plot exons
    var svg_outer = d3.select(container).append("svg")
        .attr("width", chart_width + gene_chart_margin_lower.left + gene_chart_margin_lower.right)
        .attr("height", lower_gene_chart_height)
        .attr('id', 'track_svg')
        .append("g")
        .attr('id', 'track')
        .attr("transform", "translate(" + gene_chart_margin_lower.left + "," + 0 + ")");

    var exon_color = "lightsteelblue";
    svg_outer.append("line")
        .attr("y1", lower_gene_chart_height/2)
        .attr("y2", lower_gene_chart_height/2)
        .attr("x1", 0)
        .attr("x2", exon_x_scale(coding_coordinate_params.size))
        .attr('id', 'boundary_line')
        .attr("stroke-width", 5)
        .attr("stroke", exon_color);

    // plot exon rects
    svg_outer.selectAll("bar")
        .data(_transcript.exons)
        .enter()
        .append("rect")
        .attr('class', 'track_bar')
        .style("fill", exon_color)
        .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(_transcript, d.start, true)); })
        .attr("y", function(d, i) {
            if (d.feature_type == 'CDS') {
                return 0;
            } else {
                return lower_gene_chart_height/4;
            }
        })
        .attr("width", function(d, i) {
            if (get_coding_coordinate(_transcript, d.start, true) == -100) {
                return exon_x_scale(175);
            }
            return exon_x_scale(d.stop-d.start+1);
        })
        .attr("height", function(d, i) {
            if (d.feature_type == 'CDS') {
                return lower_gene_chart_height;
            } else {
                return lower_gene_chart_height/2;
            }
        });


    var a_s = _transcript.strand == "-"? -1 : 1; //arrow direction
    var a_x = -5;  //arrow position on x-axis
    var a_y = lower_gene_chart_height/2.0; //arrow position on y-axis
    var a_w = 2; //arrow width
    var points = [[a_x+a_s*6, a_y], [a_x+a_s*1, a_y+a_w*3], [a_x+a_s*1, a_y+a_w], [a_x-a_s*9, a_y+a_w],
        [a_x-a_s*9, a_y-a_w], [a_x+a_s*1, a_y-a_w], [a_x+a_s*1, a_y-a_w*3]];
    svg_outer.append("polygon")
            .attr("points", points.join(" "))
            .attr("fill", "steelblue")
            .attr("stroke", "black");

    var bounds = get_af_bounds(variant_data);
    var min_af = bounds[0];
    var max_af = bounds[1];
    var variant_size_scale = d3.scale.log()
        .domain([min_af, max_af])
        //Circle/Ellipse
        .range([2, lower_gene_chart_height/3]);
        //Rectangle
//        .range([lower_gene_chart_height, 2]);

    // show variant category on hover
    var tip = d3.tip().attr('class', 'd3-tip').html(function(d) {
        if (d.category) {
            var csq = d.major_consequence.replace('_variant', '')
                    .replace('_', ' ')
                    .replace('utr', 'UTR')
                    .replace('3 prime', "3'")
                    .replace('5 prime', "5'")
                    .replace('nc ', "non-coding ");
            var output = csq + '<br/>' + d.chrom + ':' + d.pos + ' ' + d.ref + '&#8594;' + d.alt;
            if (d.major_consequence == 'missense_variant' || d.major_consequence == 'synonymous_variant') {
                output += '<br/>' + d.HGVSp;
            }
            if (d.allele_freq) {
                output += '<br/>Frequency: ' + d.allele_freq * 100 + '%';
            }
            return output;
        } else {
            return 'None';
        }
    });
    svg.call(tip);

    // plot variants
    svg_outer.selectAll("bar")
        .data(variant_data)
        .enter()
        .append("a")
        .attr('class', 'track_variant_link')
        .attr("xlink:href", function(d, i) { return "/" + genome + "/" + project_name + "/variant/" + d.chrom + "-" + d.pos + "-" + d.ref + "-" + d.alt; })
        .attr("data-toggle", "tooltip")
        .attr('filter_status', function(d) {
            return d.filter;
        })
        .attr('category', function(d) {
            return d.category;
        })
        .attr('major_consequence', function(d) {
            return d.major_consequence;
        })
        .attr('indel', function(d) {
            return d.indel;
        })
        .attr('variant_id', function(d) {
            return d.variant_id;
        })
        .on('mouseover', function(d) {
            $('#variant_' + d.variant_id).find('td').addClass('table_hover');
            tip.show(d);
        })
        .on('mouseout', function(d) {
            $('#variant_' + d.variant_id).find('td').removeClass('table_hover');
            tip.hide(d);
        })
        //Circle
//        .append("circle")
        //Ellipse
        .append("ellipse")
        .attr("class", function(d) {
            return "track_variant " + d.category;
        })
        .style("opacity", 0.5)
        .attr("cx",  function(d) {
            if (d['pos_coding_noutr'] == undefined) {
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        })
        .attr("cy", lower_gene_chart_height/2)
        //Circle
//        .attr("r", function(d, i) { return variant_size_scale(d.allele_freq); })
        //Ellipse
        .attr("rx", 2)
        .attr("ry", function(d, i) {
            if (!d.allele_freq) {
                return 3;
            } else {
                return variant_size_scale(d.allele_freq);
            }
        })
        // Workaround for exporting d3 to SVG (other delcaration in style.css).
        .attr('stroke', variant_colors)
        .attr('fill', variant_colors);
        //Rectangle
//        .append("rect")
//        .attr("class", "track_variant")
//        .style("fill", "darkred")
//        .style("opacity", 0.5)
//        .attr("x", function(d, i) {
//            var tx_coord = d.transcript_coordinates[transcript];
//            if (tx_coord == 0) {
//                return -1000;
//            } else {
//                var variant_exon_number = d.vep_annotations[0]['EXON'].split('/')[0] - 1;
//                return exon_x_scale(tx_coord + variant_exon_number*padding);
//            }
//        })
//        .attr("y", function(d, i) { return lower_gene_chart_height/2 - variant_size_scale(d.allele_freq)/2; } )
//        .attr("width", 2)
//        .attr("height", function(d, i) { return variant_size_scale(d.allele_freq); })
//        .attr("rx", 6)
//        .attr("ry", 6);

           // call to separate function to make the CNVs plot.
            var detail = get_plot_detail();
            if (_cnvs && _cnvs.length > 0) make_cnvs_svg(_cnvs, _transcript, detail, true, container, cnv_svg);
}

function variant_colors(d) {
    if (d.category == 'lof_variant') {
        return '#cd2932';
    } else if (d.category == 'missense_variant') {
        return '#a96500';
    } else if (d.category == 'synonymous_variant') {
        return '#157e28';
    }
}

function change_coverage_chart(data, new_data, variant_data, _transcript, scale_type, metric, skip_utrs, container, cnv_svg, _cnvs) {
    var coords = skip_utrs ? 'pos_coding_noutr' : 'pos_coding';
    var max_cov = (metric == 'mean' || metric == 'median') ? d3.max(data, function(d) { return d[metric]; }) : 1;
    var coding_coordinate_params = get_coding_coordinate_params(_transcript, skip_utrs);
    var chart_width;
    if (scale_type == 'overview') {
        chart_width = gene_chart_width;
    } else {
        chart_width = coding_coordinate_params.size*2;
    }

    var exon_x_scale = d3.scale.linear()
        .domain([0, coding_coordinate_params.size])
        .range([0, chart_width]);

    var svg = d3.select(container).select('#inner_svg')
        .attr("width", chart_width + gene_chart_margin.left + gene_chart_margin.right)
        .attr("height", gene_chart_height + gene_chart_margin.top + gene_chart_margin.bottom)
        .select('#inner_graph');

    var y = d3.scale.linear()
        .domain([0, max_cov])
        .range([gene_chart_height, 0]);

    var area = d3.svg.area()
        .x( function(d) {
            return exon_x_scale(d[coords]);
        }).y0( function(d) {
            return gene_chart_height;
        }).y1( function(d) {
            return (metric in d) ? y(d[metric]) : gene_chart_height;
        });

    var path = svg.selectAll("path")
        .datum(new_data)
        .transition()
        .duration(500)
        .attr("d", area)
        .style("fill", "steelblue")
        .attr('class', 'area');

    svg.selectAll("circle")
        .data(data)
        .transition()
        .duration(500)
        .attr("cx", function(d) {
            if (d[coords] == undefined) {
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        });

    // plot exons
    var svg_outer = d3.select(container).select('#track_svg')
        .attr("width", chart_width + gene_chart_margin_lower.left + gene_chart_margin_lower.right)
        .attr("height", lower_gene_chart_height).select('#track');

    svg_outer.select('#boundary_line')
        .attr("x2", exon_x_scale(coding_coordinate_params.size));

    // plot exon rounded rects
    svg_outer.selectAll("rect")
        .data(_transcript.exons)
        .transition()
        .duration(500)
        .attr("x", function(d, i) { return exon_x_scale(get_coding_coordinate(_transcript, d.start, skip_utrs)); })
        .attr("width", function(d, i) {
            if (get_coding_coordinate(_transcript, d.start, skip_utrs) == -100) {
                return exon_x_scale(175);
            }
            return exon_x_scale(d.stop-d.start+1);
        })
        .attr("height", function(d, i) {
            if (d.feature_type == 'CDS') {
                return lower_gene_chart_height;
            } else {
                return lower_gene_chart_height/2;
            }
        });

    // plot variants
    svg_outer.selectAll("a")
        .data(variant_data)
        .transition()
        .duration(500)
        .selectAll('ellipse')
        .attr("cx", function(d) {
            if (d[coords] == undefined) {
                return -100;
            } else {
                return exon_x_scale(d[coords]);
            }
        });

    var yAxis = d3.svg.axis()
        .scale(y)
        .orient("left");

    svg.select(".y.axis")
        .transition()
        .duration(200)
        .call(yAxis);

    $("#" + cnv_svg).remove();
    if (_cnvs && _cnvs.length > 0)
        make_cnvs_svg(_cnvs, _transcript, scale_type, skip_utrs, container, cnv_svg);
}

function create_new_data(data, coords) {
    var data_object = {};
    $.each(data, function(i, d) {
        data_object[d[coords]] = d;
    });
    var new_data = [];
    var max_coord = d3.max(data, function(d) { return d[coords] });
    for (var i = d3.min(data, function(d) { return d[coords] }); i < max_coord; i++) {
        var x = {'has_coverage': false};
        x[coords] = i;
        //Check the previous base to see if this is the beginning of an exon
        if (i in data_object && !(i-1 in data_object)) {
            new_data.push(x);
        }
        //Check the previous base to see if this is the end of an exon
        if (!(i in data_object) && i-1 in data_object) {
            x[coords] = i-1;
            new_data.push(x);
        }
        if (i in data_object) {
            new_data.push(data_object[i]);
        } else {
            new_data.push(x);
        }
    }
    return new_data;
}

function coverage_sum(key, coverage_stats) {
    var total = 0;
    $.map(coverage_stats, function(entry) {
        total += entry[key];
    });
    return (total / coverage_stats.length).toPrecision(4);
}

function get_plot_detail() {
    var details = $('.display_coverage_metric_buttons.active')
    var detail;
    if(details.length === 0){
        detail = 'overview'
    }
    else{
        detail = details.attr('id').replace('display_coverage_', '').replace('_button', '');
    }
    return detail;
}


function refresh_links(plot_id, _cnvs) {
    var gene_plot_container = "gene_plot_container_" + plot_id;
    $("#coverage_plot_download_" + plot_id).attr('href', set_plot_image(gene_plot_container, 0));
    $("#exon_plot_download_" + plot_id).attr('href', set_plot_image(gene_plot_container, 1));
    if (_cnvs && _cnvs.length > 0)
        $("#cnv_plot_download_" + plot_id).attr('href', set_plot_image(gene_plot_container, 2));
}