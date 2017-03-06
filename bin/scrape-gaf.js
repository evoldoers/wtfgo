#!/usr/bin/env node

var fs = require('fs')
var readline = require('readline')
var Promise = require('bluebird')
var exec = Promise.promisify (require('child_process').exec)

var gaf2json = require('../lib/converters').gaf2json
var obo2json = require('../lib/converters').obo2json

var go_url_prefix = "http://geneontology.org/ontology/"
var gaf_url_prefix = "http://geneontology.org/gene-associations/"

var go_url_suffix = "go-basic.obo"
var gaf_metadata_url_suffix = "go_annotation_metadata.all.js"

var uniprot_mapping_url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz"

var metadata_url = gaf_url_prefix + gaf_metadata_url_suffix

function skip_id (id) {
  return /^goa_uniprot_all/.test(id) || /_complex$/.test(id)
}

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
       .then (() => exec("mkdir " + deploy_dir + '/' + gaf_subdir))
       .then (() => exec("mkdir " + deploy_dir + '/' + go_subdir))
       .then (() => exec("mkdir " + download_dir))
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
init_promise
  .then (() => download_data (go_url_prefix + go_url_suffix))
  .then (function (obo) {
    console.log ("processing " + go_url_suffix)
    var go_json = obo2json ({ obo: obo,
			      compress: true,
			      includeTermInfo: true })
    fs.writeFileSync (deploy_dir + '/' + go_path, JSON.stringify (go_json))
  })
// download UniProt mapping
  .then (() => download_filename (uniprot_mapping_url))
  .then ((uniprot_mapping_filename) => download_data (metadata_url)
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
                .then ((gaf_file) => get_aliases (gaf_filename, gaf_file, uniprot_mapping_filename))
                .then (function (gaf_json) {
                  // process
                  var gaf_path = gaf_subdir + '/' + resource.id + '.json'
                  fs.writeFileSync (deploy_dir + '/' + gaf_path, JSON.stringify (gaf_json))
                  return {
                    name: resource.label,
                    ontologies: [
                      { name: "Gene Ontology (basic)",
                        ontology: go_path,
                        assocs: gaf_path,
                        examples: [] }
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

function get_aliases (gaf_filename, gaf_file, uniprot_mapping_filename) {
  return new Promise (function (resolve, reject) {
    console.log ("processing " + gaf_filename)
    var gaf_json_no_aliases = gaf2json ({ gaf: gaf_file })
    // find all gene symbols
    var is_symbol = {}
    gaf_json_no_aliases.idAliasTerm.forEach (function (id_alias_term) {
      is_symbol[id_alias_term[0]] = true
      id_alias_term[1].forEach (function (alias) { is_symbol[alias] = true })
    })
    // scan through Uniprot ID map for aliases
    var alias = []
    var line_reader = readline.createInterface({
      input: require('fs').createReadStream(uniprot_mapping_filename)
    })
    var re = /^(.*?) /;
    line_reader.on ('line', function (line) {
      var match = re.exec (line)
      if (match && is_symbol[match[1]])
        alias.push (line)
    })
    line_reader.on ('close', function() {
      console.log ("finished " + gaf_filename + " (" + alias.length + " aliases)")
      var aliases = alias.join('')
      var gaf_json = gaf2json ({ aliases: aliases, gaf: gaf_file })
      resolve (gaf_json)
    })
  })
}

