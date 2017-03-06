#!/usr/bin/env node

var fs = require('fs')
var readline = require('readline')
var Promise = require('bluebird')
var exec = Promise.promisify (require('child_process').exec)

var converters = require('../wtfgenes/lib/converters')
var Ontology = require('../wtfgenes/lib/ontology')
var gaf2json = converters.gaf2json
var obo2json = converters.obo2json

var go_url_prefix = "http://geneontology.org/ontology/"
var gaf_url_prefix = "http://geneontology.org/gene-associations/"

var go_url_suffix = "go-basic.obo"
var gaf_metadata_url_suffix = "go_annotation_metadata.all.js"

var metadata_url = gaf_url_prefix + gaf_metadata_url_suffix

function skip_id (id) {
  return /^goa_uniprot_all/.test(id) || /_complex$/.test(id) || /_rna$/.test(id) || id === 'jcvi'
}

var example_terms = [['GO:0006298'],  // mismatch repair
                     ['GO:0000027'],  // ribosomal large subunit assembly
                     ['GO:0004803','GO:0006302']]  // transposase activity + double-strand break repair

var deploy_dir = 'web'
var gaf_subdir = 'gaf'
var go_subdir = 'go'

var download_dir = 'download'

var go_path = go_subdir + '/' + go_url_suffix.replace('.obo','.json')

// set up deploy directory
var init_promise = fs.existsSync(deploy_dir)
    ? Promise.resolve(true)
    : (exec("cp -r wtfgenes/web " + deploy_dir)
       .then (() => exec("cd " + deploy_dir + "; make"))
       .then (() => exec("mkdir -p " + deploy_dir + '/' + gaf_subdir))
       .then (() => exec("mkdir -p " + deploy_dir + '/' + go_subdir))
       .then (() => exec("mkdir -p " + download_dir))
       .then (() => console.log("initialized " + deploy_dir + " and " + download_dir)))

// wrappers to download-and-cache
var readFile = Promise.promisify (fs.readFile)
function download_filename (url) {
  var filename = url.replace(/^.*\//,'')
  var local_path = download_dir + '/' + filename
  var local_unzipped = local_path.replace (/\.gz$/, '')
  var curl_promise = fs.existsSync (local_path)
      ? Promise.resolve(true)
      : (exec ("cd " + download_dir + "; curl -O " + url)
         .then (() => {
           if (url.match(/\.gz$/))
             return exec ("cd " + download_dir + "; gunzip --keep " + filename)
           else
             return true
         }))
  return curl_promise
    .then (() => local_unzipped)
}

function download_data (url) {
  return download_filename (url)
    .then ((filename) => readFile (filename))
    .then ((buffer) => buffer.toString())
}

// download GO
var term_name = {}, in_term_closure = {}
init_promise
  .then (() => download_data (go_url_prefix + go_url_suffix))
  .then (function (obo) {
    console.log ("processing " + go_url_suffix)
    var go_json = obo2json ({ obo: obo,
			      compress: true,
			      includeTermInfo: true })
    fs.writeFileSync (deploy_dir + '/' + go_path, JSON.stringify (go_json))
    var go = new Ontology (go_json), trans = go.transitiveClosure()
    go_json.termParents.forEach (function (term_parents, n) {
      var name = term_parents[0]
      term_name[name] = go_json.termInfo[n]
      var in_closure = {}
      trans[n].forEach ((c) => { in_closure[go_json.termParents[c][0]] = true })
      in_term_closure[name] = in_closure
    })
  })
  .then (() => download_data (metadata_url)
         // download GAF metadata
         .then (function (html) {
           eval (html)  // defines global_go_annotation_metadata
           return global_go_annotation_metadata
         }).then (function (metadata) {
           var resources = metadata.resources
           return Promise.all
           (resources
            .filter (resource => !skip_id(resource.id))
            .map (function (resource) {
              // download each GAF file
              var gaf_filename = resource.gaf_filename
              return download_data (gaf_url_prefix + gaf_filename)
                .then ((gaf_file) => quick_gaf2json (gaf_filename, gaf_file))
                .then (function (gaf_json) {
                  // examples
                  var examples = example_terms.map (function (terms) {
                    var found_term = {}
                    var genes = gaf_json.idAliasTerm.filter (function (id_alias_term) {
                      var has_term = false
                      id_alias_term[2].forEach (function (gterm) {
                        terms.forEach (function (term) {
                          if (in_term_closure[gterm][term])
                            has_term = found_term[term] = true
                        })
                      })
                      return has_term
                    }).map (function (id_alias_term) { return id_alias_term[0] })
                    var name = terms.filter ((term) => found_term[term])
                        .map ((term) => term_name[term]).join(" + ")
                    return { name, genes }
                  }).filter (function (example) { return example.genes.length > 0 })
                  // write
                  var gaf_path = gaf_subdir + '/' + resource.id + '.json'
                  fs.writeFileSync (deploy_dir + '/' + gaf_path, JSON.stringify (gaf_json))
                  return {
                    name: resource.label,
                    ontologies: [
                      { name: "Gene Ontology (basic)",
                        ontology: go_path,
                        assocs: gaf_path,
                        examples: examples }
                    ]
                  }
                })
            }))
         }).then (function (organisms) {
           organisms = organisms.sort (function (a, b) { return a.name < b.name ? -1 : +1 })
           fs.writeFileSync (deploy_dir + '/datasets.json', JSON.stringify ({ organisms }))
           console.log ("done")
         })
        )

function quick_gaf2json (gaf_filename, gaf_file) {
  console.log ("processing " + gaf_filename)
  var json = gaf2json ({ gaf: gaf_file,
                         mergeDuplicates: true })
  return json
}

