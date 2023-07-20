from typing import Dict
from collections import OrderedDict

from jsonpath_ng import parse

import rdflib
from pyld import jsonld
from rdflib.namespace import Namespace, RDF, RDFS, XSD
import json
from rdflib import BNode, Literal

# Namespaces
DCTERMS = Namespace('http://purl.org/dc/terms/')
NSG = Namespace('https://neuroshapes.org/')
SKOS = Namespace('http://www.w3.org/2004/02/skos/core#')
NIFRID = Namespace('http://uri.neuinfo.org/nif/nifstd/readable/')
ILXTR = Namespace('http://uri.interlex.org/tgbugs/uris/readable/')
OWL = Namespace('http://www.w3.org/2002/07/owl#')
SCHEMA = Namespace('http://schema.org/')
mba = Namespace('http://api.brain-map.org/api/v2/data/Structure/')
PAXRAT = Namespace('http://uri.interlex.org/paxinos/uris/rat/labels/')
PROV = Namespace('http://www.w3.org/ns/prov#')
UBERON = Namespace('http://purl.obolibrary.org/obo/UBERON_')
BMO = Namespace("https://bbp.epfl.ch/ontologies/core/bmo/")
TYPE_FOR_FRAMING = NSG["FrameType"]

ROOT_REGION = 997
CHILDREN = "children"

REPRESENTED = "representedInAnnotation"
REGIONVOLUME = "regionVolume"
REGIONVOLUMERATIO = "regionVolumeRatioToWholeBrain"

def hierarchy_json_to_jsonld(json_content: Dict) -> Dict:
    jsonpath_expr = parse(f'$..{CHILDREN}.[*]') # missing the root (997) node
    matches = jsonpath_expr.find(json_content)
    match_map = {str(ROOT_REGION): json_content}
    for match in matches:
        match_map[str(match.value['id'])] = match.value

    hierarchy_graph = rdflib.Graph()
    # Todo: provide a non blank node URI
    hierarchy_uri = BNode()
    hierarchy_graph.add((hierarchy_uri, RDF.type, OWL.Ontology))

    # init prefix mappings
    hierarchy_graph.bind('dcterms', DCTERMS)
    hierarchy_graph.bind('nsg', NSG)
    hierarchy_graph.bind('defines', NSG["defines"])
    hierarchy_graph.bind('uberon', UBERON)
    hierarchy_graph.bind('skos', SKOS)
    hierarchy_graph.bind('label', RDFS["label"])
    hierarchy_graph.bind('subClassOf', RDFS["subClassOf"])
    hierarchy_graph.bind('notation', SKOS["notation"])
    hierarchy_graph.bind('prefLabel', SKOS["prefLabel"])
    hierarchy_graph.bind('altLabel', SKOS["altLabel"])
    hierarchy_graph.bind('NIFRID', NIFRID)
    hierarchy_graph.bind('owl', OWL)
    hierarchy_graph.bind('schema', SCHEMA)
    hierarchy_graph.bind('hasPart', SCHEMA["hasPart"])
    hierarchy_graph.bind('isPartOf', SCHEMA["isPartOf"])
    hierarchy_graph.bind('identifier', SCHEMA["identifier"])
    hierarchy_graph.bind('isDefinedBy', RDFS["isDefinedBy"])
    hierarchy_graph.bind('mba', mba)
    hierarchy_graph.bind('PROV', PROV)
    hierarchy_graph.bind('BMO', BMO)

    for _, v in match_map.items():
        add_class_to_graph(v, hierarchy_graph, hierarchy_uri)

    hierarchy_graph_str = hierarchy_graph.serialize(format="json-ld", auto_compact=True, indent=2)
    hierarchy_graph_jsonld = json.loads(hierarchy_graph_str)

    if '@context' in hierarchy_graph_jsonld:
        context = hierarchy_graph_jsonld["@context"]
        frame_json = {
            "@context": context,
            "@type": str(OWL.Ontology),
            "contains": [{
                "@type": str(TYPE_FOR_FRAMING),
                "@embed": True
            }]
        }
        hierarchy_graph_framed = jsonld.frame(json.loads(hierarchy_graph_str), frame_json)

        hierarchy_framed = graph_free_jsonld(hierarchy_graph_framed)
        hierarchy_framed["defines"][0]["@type"] = "owl:Class" # it's already the case
        #hierarchy_framed["defines"] = [hierarchy_framed["defines"]]
        if "nsg:tempDefines" in hierarchy_framed:
            defined_ids = filter(lambda x: x["@id"] != "mba:"+str(ROOT_REGION), hierarchy_framed["nsg:tempDefines"])
            hierarchy_framed["defines"].extend(defined_ids)
            hierarchy_framed.pop("nsg:tempDefines")
        return hierarchy_framed

    else:
        return {}


def add_class_to_graph(json_content, hierarchy_graph, hierarchy_uri):
    ss = mba[str(json_content["id"])]
    hierarchy_graph.add((ss, RDFS.subClassOf, NSG.BrainRegion))

    if json_content["id"] == ROOT_REGION:
        hierarchy_graph.add((ss, RDF.type, TYPE_FOR_FRAMING))
    else:
        hierarchy_graph.add((ss, RDF.type, OWL.Class))
    hierarchy_graph.add((hierarchy_uri, NSG.defines, ss))
    hierarchy_graph.add((hierarchy_uri, NSG.tempDefines, ss))

    hierarchy_graph.add((ss, mba.color_hex_triplet, Literal(json_content["color_hex_triplet"])))
    hierarchy_graph.add((ss, mba.hemisphere_id, Literal(json_content["hemisphere_id"])))
    hierarchy_graph.add((ss, BMO.hasLayerLocationPhenotype, Literal(json_content["layers"])))
    hierarchy_graph.add((ss, BMO.representedInAnnotation, Literal(json_content[REPRESENTED])))
    if json_content[REPRESENTED]:
        hierarchy_graph.add((ss, BMO.regionVolume, Literal(json_content[REGIONVOLUME])))
        hierarchy_graph.add((ss, BMO.regionVolumeRatioToWholeBrain, Literal(json_content[REGIONVOLUMERATIO])))
        hierarchy_graph.add((ss, BMO.adjacentTo, Literal(json_content["adjacentTo"])))
        hierarchy_graph.add((ss, BMO.continuousWith, Literal(json_content["continuousWith"])))

    hierarchy_graph.add((SCHEMA.isPartOf, RDF.type, OWL.ObjectProperty))
    hierarchy_graph.add((SCHEMA.hasPart, RDF.type, OWL.ObjectProperty))

    children = json_content[CHILDREN]

    for child in children:
        child_fragment = child["id"]
        hierarchy_graph.add((ss, SCHEMA.hasPart, mba[str(child_fragment)]))
        hierarchy_graph.add((mba[str(child_fragment)], SCHEMA.isPartOf, ss))

    hierarchy_graph.add((ss, SKOS.prefLabel, Literal(json_content["name"])))
    hierarchy_graph.add((ss, RDFS.label, Literal(json_content["name"])))
    hierarchy_graph.add((ss, RDFS.isDefinedBy, hierarchy_uri))

    hierarchy_graph.add((ss, SKOS.notation, Literal(json_content["acronym"])))
    hierarchy_graph.add((ss, SCHEMA.identifier, Literal(str(json_content["id"]), datatype=XSD.string)))
    if json_content["name"] != json_content["acronym"]:
        hierarchy_graph.add((ss, SKOS.altLabel, Literal(json_content["acronym"])))


def graph_free_jsonld(jsonld_doc, context=None):
    if "@graph" in jsonld_doc and len(jsonld_doc["@graph"]) > 0:
        graph_free_jsonld_doc = jsonld_doc["@graph"][0]
        if not context:
            context = jsonld_doc["@context"]
        graph_free_jsonld_doc['@context'] = context
        graph_free = OrderedDict(graph_free_jsonld_doc)
        graph_free['@context'] = context
        graph_free.move_to_end('@context', last=False)

        return graph_free
    else:
        return jsonld_doc
